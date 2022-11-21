#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    extract_b0;
    extract_b0 as aff_extract_b0;
    extract_b0 as syn_extract_b0;
    extract_b0 as extract_target_b0;
    compute_powder_average as aff_pa_dwi;
    compute_powder_average as syn_pa_dwi;
    compute_powder_average as compute_target_pdavg
} from '../processes/preprocess.nf'
include {
    scil_compute_dti_fa as dti_fa_eroded;
    scil_compute_dti_fa as dti_fa
} from '../processes/measure.nf'
include {
    registration_wkf as t1_to_reference_affine;
    registration_wkf as b0_to_reference_affine;
    registration_wkf as t1_to_reference_syn;
    registration_wkf as b0_to_reference_syn;
    registration_wkf as t1_to_b0_registration_wkf
} from "./preprocess.nf"
include {
    ants_transform as transform_t1_to_b0;
    ants_transform as transform_mask_to_b0;
    ants_transform as transform_t1_to_reference;
    ants_transform as transform_t1_mask_to_reference;
    ants_transform as transform_b0_to_reference;
    ants_transform as transform_dwi_mask_to_reference;
    ants_transform as transform_dwi_to_reference;
    ants_transform as transform_t1_syn;
    ants_transform as transform_t1_mask_syn;
    ants_transform as transform_b0_syn;
    ants_transform as transform_dwi_mask_syn;
    ants_transform as transform_dwi_syn;
    ants_transform as transform_t1_mask_to_b0
} from '../processes/register.nf'
include {
    prepend_sid as prepend_sid_template;
    prepend_sid as prepend_sid_template_mask;
    prepend_sid as prepend_sid_template_dilated_mask;
    clean_mask_borders;
    dilate_mask as dilate_t1_mask;
    dilate_mask as dilate_dwi_mask;
    erode_mask as erode_dwi_mask;
    apply_mask as mask_b0_dilated;
    apply_mask as mask_pa_dwi_dilated;
    apply_mask as mask_t1_dilated;
    apply_mask as mask_b0;
    apply_mask as mask_pa_dwi;
    apply_mask as mask_t1;
    apply_mask as mask_target_b0;
    apply_mask as mask_target_pdavg;
    apply_mask as mask_moving_t1;
    apply_mask as mask_template;
    bet_mask
} from '../processes/utils.nf'
include {
    resampling_reference;
    scilpy_resample_to_reference as resample_template;
    scilpy_resample_to_reference as resample_dilated_mask
} from '../processes/upsample.nf'
include {
    get_data_path;
    get_config_path
} from "../functions.nf"


params.use_quick = false
params.tissue_segmentation_root = "${get_data_path()}/maccaca_mulatta/tissue_segmentation"

params.t1_registration_extract_b0_config = file("${get_config_path()}/extract_mean_b0_base_config.py")
params.t1_to_b0_registration_config = file("${get_config_path()}/t1_to_b0_registration_affine_config.py")
params.ants_transform_mask_config = file("${get_config_path()}/ants_transform_mask_config.py")
params.ants_transform_base_config = file("${get_config_path()}/ants_transform_base_config.py")

params.t1_to_template_affine_config = file("${get_config_path()}/t1_to_template_affine_config.py")
params.b0_to_template_affine_config = file("${get_config_path()}/b0_to_template_affine_config.py")
params.t1_to_template_affine_quick_config = file("${get_config_path()}/t1_to_template_affine_quick_config.py")
params.b0_to_template_affine_quick_config = file("${get_config_path()}/b0_to_template_affine_quick_config.py")

params.t1_to_template_syn_config = file("${get_config_path()}/t1_to_template_syn_config.py")
params.b0_to_template_syn_config = file("${get_config_path()}/b0_to_template_syn_config.py")
params.t1_to_template_syn_quick_config = file("${get_config_path()}/t1_to_template_syn_quick_config.py")
params.b0_to_template_syn_quick_config = file("${get_config_path()}/b0_to_template_dyn_quick_config.py")


workflow t12b0_registration {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        dwi_mask_channel
        dwi_metadata_channel
        publish_mask
        publish_t1
        use_quick
    main:
        extract_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)

        template_channel = prepend_sid_template(
            t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] }
        )
        template_mask_channel = prepend_sid_template_mask(
            t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_mask_no_bv.nii.gz")] }
        )
        template_dilated_mask_channel = prepend_sid_template_dilated_mask(
            t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_mask_no_bv_dilated.nii.gz")] }
        )

        resampling_reference(dwi_channel.map{ it.subList(0, 2) }.join(t1_channel).join(template_channel).map{ [it[0], it[1..-1]] }, "preprocess")
        resample_template(
            template_channel
                .join(resampling_reference.out.reference)
                .join(template_mask_channel)
                .map{ it + [""] },
            "preprocess",
            "lin",
            false,
            false,
            "", ""
        )
        resample_dilated_mask(
            template_dilated_mask_channel
                .join(resampling_reference.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "nn",
            false,
            false,
            "", ""
        )
        mask_template(
            resample_template.out.image
                .join(resample_template.out.mask)
                .map{ it + [""] },
            "preprocess",
            false
        )

        t1_to_b0_affine(
            dwi_channel,
            t1_channel,
            dwi_mask_channel,
            t1_mask_channel,
            resample_template.out.image,
            resample_dilated_mask.out.image,
            dwi_metadata_channel,
            use_quick ? params.t1_to_template_affine_quick_config : params.t1_to_template_affine_config,
            use_quick ? params.b0_to_template_affine_quick_config : params.b0_to_template_affine_config,
            false,
            false
        )

        t1_to_b0_syn(
            t1_to_b0_affine.out.dwi,
            t1_to_b0_affine.out.t1,
            t1_to_b0_affine.out.dwi_mask,
            t1_to_b0_affine.out.t1_mask,
            mask_template.out.image,
            resample_template.out.mask,
            dwi_metadata_channel,
            use_quick ? params.t1_to_template_syn_quick_config : params.t1_to_template_syn_config,
            use_quick ? params.b0_to_template_syn_quick_config : params.b0_to_template_syn_config,
            false,
            false
        )

        template_to_b0_transform = t1_to_b0_affine.out.b0_inverse_transform
            .join(t1_to_b0_syn.out.b0_inverse_transform)
            .map{ [it[0], it[1] + it[3], it[2] + it[4]] }

        t1_to_b0_transform = template_to_b0_transform
            .join(t1_to_b0_syn.out.t1_transform)
            .map{ [it[0], it[1] + it[3], it[2] + it[4]] }
            .join(t1_to_b0_affine.out.t1_transform)
            .map{ [it[0], it[1] + it[3], it[2] + it[4]] }

        transform_t1_to_b0(
            t1_channel
                .join(extract_b0.out.b0)
                .join(t1_to_b0_transform)
                .view()
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_t1", "t1",
            params.ants_transform_base_config
        )

        transform_mask_to_b0(
            t1_to_b0_syn.out.t1_mask
                .join(extract_b0.out.b0)
                .join(template_to_b0_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_mask_config
        )

        clean_mask_borders(transform_mask_to_b0.out.image, 2, "preprocess")

    emit:
        t1 = transform_t1_to_b0.out.image
        mask = clean_mask_borders.out.mask
        transform = t1_to_b0_transform
        reference = extract_b0.out.b0
}

workflow t1_to_b0_affine {
    take:
        dwi_channel
        t1_channel
        dwi_mask_channel
        t1_mask_channel
        reference_channel
        reference_dilated_mask_channel
        dwi_metadata_channel
        t1_affine_config
        b0_affine_config
        publish_t1
        publish_b0
    main:
        dilate_t1_mask(t1_mask_channel, 8, "preprocess")
        dilate_dwi_mask(dwi_mask_channel, 8, "preprocess")
        erode_dwi_mask(dwi_mask_channel, 8, "preprocess")

        aff_extract_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)
        mask_b0_dilated(aff_extract_b0.out.b0.join(dilate_dwi_mask.out.mask).map{ it + [""] }, "preprocess", "false")
        aff_pa_dwi(dwi_channel.map{ it.subList(0, 3) }.map{ it + ["", ""] }, "preprocess", "false")
        mask_pa_dwi_dilated(aff_pa_dwi.out.image.join(dilate_dwi_mask.out.mask).map{ it + [""] }, "preprocess", "false")
        dti_fa_eroded(dwi_channel.join(erode_dwi_mask.out.mask), "preprocess", "preprocess", false)
        mask_t1_dilated(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")

        b0_moving_channel = mask_b0_dilated.out.image
            .join(mask_pa_dwi_dilated.out.image)
            .join(dti_fa_eroded.out.fa)
            .map{ [it[0], it[1..-1]] }
        t1_moving_channel = mask_t1_dilated.out.image
            .map{ [it[0], [it[1]]] }

        t1_to_reference_affine(
            reference_channel.map{ [it[0], [it[1]]] },
            t1_moving_channel,
            null,
            reference_dilated_mask_channel.join(dilate_t1_mask.out.mask).map{ [it[0], it[1..-1]] },
            null,
            null,
            "",
            false,
            "",
            "",
            t1_affine_config,
            params.ants_transform_mask_config
        )

        b0_to_reference_affine(
            reference_channel.map{ [it[0], [it[1]]] },
            b0_moving_channel,
            null,
            reference_dilated_mask_channel.join(dilate_dwi_mask.out.mask).map{ [it[0], it[1..-1]] },
            null,
            null,
            "",
            false,
            "",
            "",
            b0_affine_config,
            params.ants_transform_mask_config
        )

        t1_transform = t1_to_reference_affine.out.transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ "false" }] }
        b0_transform = b0_to_reference_affine.out.transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ "false" }] }

        transform_t1_to_reference(
            t1_channel
                .join(t1_to_reference_affine.out.reference)
                .join(t1_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_t1", "t1",
            params.ants_transform_base_config
        )

        transform_t1_mask_to_reference(
            t1_mask_channel
                .join(t1_to_reference_affine.out.reference)
                .join(t1_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_t1", "t1_mask",
            params.ants_transform_mask_config
        )

        transform_b0_to_reference(
            aff_extract_b0.out.b0
                .join(b0_to_reference_affine.out.reference)
                .join(b0_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_b0", "b0",
            params.ants_transform_base_config
        )

        transform_dwi_mask_to_reference(
            dwi_mask_channel
                .join(b0_to_reference_affine.out.reference)
                .join(b0_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_b0", "dwi_mask",
            params.ants_transform_mask_config
        )

        transform_dwi_to_reference(
            dwi_channel.map{ it[0..1] }
                .join(b0_to_reference_affine.out.reference)
                .join(b0_transform)
                .join(dwi_channel.map{ [it[0], it[3]] })
                .join(dwi_metadata_channel),
            "preprocess",
            "","$publish_b0", "dwi",
            params.ants_transform_base_config
        )

    emit:
        t1 = transform_t1_to_reference.out.image
        t1_mask = transform_t1_mask_to_reference.out.image
        b0 = transform_b0_to_reference.out.image
        dwi = transform_dwi_to_reference.out.image
            .join(dwi_channel.map{ [it[0], it[2]] })
            .join(transform_dwi_to_reference.out.bvec)
        dwi_mask = transform_dwi_mask_to_reference.out.image
        dwi_metadata = transform_dwi_to_reference.out.metadata
        t1_transform = t1_transform
        t1_inverse_transform = t1_to_reference_affine.out.inverse_transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ f -> f.name.contains("syn_inverse") ? "false": "true" }] }
        b0_transform = b0_transform
        b0_inverse_transform = b0_to_reference_affine.out.inverse_transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ f -> f.name.contains("syn_inverse") ? "false": "true" }] }
}

workflow t1_to_b0_syn {
    take:
        dwi_channel
        t1_channel
        dwi_mask_channel
        t1_mask_channel
        reference_channel
        reference_mask_channel
        dwi_metadata_channel
        t1_syn_config    
        b0_syn_config
        publish_t1
        publish_b0
    main:
        syn_extract_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)
        mask_b0(syn_extract_b0.out.b0.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        syn_pa_dwi(dwi_channel.map{ it.subList(0, 3) }.map{ it + ["", ""] }, "preprocess", "false")
        mask_pa_dwi(syn_pa_dwi.out.image.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess", false)
        mask_t1(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")

        b0_moving_channel = mask_b0.out.image
            .join(mask_pa_dwi.out.image)
            .join(dti_fa.out.fa)
            .map{ [it[0], it[1..-1]] }
        t1_moving_channel = mask_t1.out.image
            .map{ [it[0], [it[1]]] }
        reference_fixed_channel = reference_channel.map{ [it[0], [it[1]]] }

        t1_to_reference_syn(
            reference_fixed_channel,
            t1_moving_channel,
            t1_mask_channel,
            reference_mask_channel.join(t1_mask_channel).map{ [it[0], it[1..-1]] },
            null,
            null,
            "",
            false,
            "",
            "",
            t1_syn_config,
            params.ants_transform_mask_config
        )

        b0_to_reference_syn(
            reference_fixed_channel,
            b0_moving_channel,
            null,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            b0_syn_config,
            params.ants_transform_mask_config
        )

        t1_transform = t1_to_reference_syn.out.transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ "false" }] }
        b0_transform = b0_to_reference_syn.out.transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ "false" }] }

        transform_t1_syn(
            t1_channel
                .join(t1_to_reference_syn.out.reference)
                .join(t1_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_t1", "t1",
            params.ants_transform_base_config
        )

        transform_t1_mask_syn(
            t1_mask_channel
                .join(t1_to_reference_syn.out.reference)
                .join(t1_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_t1", "t1_mask",
            params.ants_transform_mask_config
        )

        transform_b0_syn(
            syn_extract_b0.out.b0
                .join(b0_to_reference_syn.out.reference)
                .join(b0_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_b0", "b0",
            params.ants_transform_base_config
        )

        transform_dwi_mask_syn(
            dwi_mask_channel
                .join(b0_to_reference_syn.out.reference)
                .join(b0_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_b0", "dwi_mask",
            params.ants_transform_mask_config
        )

        transform_dwi_syn(
            dwi_channel.map{ it[0..1] }
                .join(b0_to_reference_syn.out.reference)
                .join(b0_transform)
                .map{ it + [""] }
                .join(dwi_metadata_channel),
            "preprocess",
            "","$publish_b0", "dwi",
            params.ants_transform_base_config
        )

    emit:
        t1 = transform_t1_syn.out.image
        t1_mask = transform_t1_mask_syn.out.image
        b0 = transform_b0_syn.out.image
        dwi = transform_dwi_syn.out.image
            .join(dwi_channel.map{ [it[0], it[2], it[3]] })
        dwi_mask = transform_dwi_mask_syn.out.image
        dwi_metadata = transform_dwi_syn.out.metadata
        t1_transform = t1_transform
        t1_inverse_transform = t1_to_reference_syn.out.inverse_transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ f -> f.name.contains("syn_inverse") ? "false": "true" }] }
        b0_transform = b0_transform
        b0_inverse_transform = b0_to_reference_syn.out.inverse_transform
            .map{ [it[0], [it[1]].flatten()] }
            .map{ [it[0], it[1], it[1].collect{ f -> f.name.contains("syn_inverse") ? "false": "true" }] }
}

workflow t1_mask_to_b0 {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        publish_mask
    main:
        extract_target_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)

        bet_mask(extract_target_b0.out.b0, "preprocess", "false")
        dwi_mask_channel = bet_mask.out.mask

        mask_target_b0(extract_target_b0.out.b0.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")

        compute_target_pdavg(dwi_channel.map{ it.subList(0, 3) }.map{ it + ["", ""] }, "preprocess", "false")
        mask_target_pdavg(compute_target_pdavg.out.image.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")

        mask_moving_t1(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")

        target_channel = mask_target_b0.out.image
            .join(mask_target_pdavg.out.image)
            .map{ [it[0], it[1..-1]] }
        moving_channel = mask_moving_t1.out.image
            .map{ [it[0], [it[1]]] }

        t1_to_b0_registration_wkf(
            target_channel,
            moving_channel,
            null,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            params.t1_to_b0_registration_config,
            ""
        )
        transform_t1_mask_to_b0(
            t1_mask_channel
                .join(extract_target_b0.out.b0)
                .join(t1_to_b0_registration_wkf.out.transform)
                .map{ it + ["", "", ""] },
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_mask_config
        )
    emit:
        mask = transform_t1_mask_to_b0.out.image
}