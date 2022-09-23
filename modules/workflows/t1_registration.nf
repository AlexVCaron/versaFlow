#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    extract_b0;
    extract_b0 as extract_target_b0;
    compute_powder_average;
    compute_powder_average as compute_target_pdavg
} from '../processes/preprocess.nf'
include { 
    scil_compute_dti_fa;
    scil_compute_dti_fa as compute_target_fa
} from '../processes/measure.nf'
include {
    registration_wkf as t1_to_template_registration_wkf;
    registration_wkf as b0_to_template_registration_wkf;
    registration_wkf as t1_to_b0_registration_wkf
} from "./preprocess.nf"
include {
    ants_transform as ants_transform_t1_to_b0;
    ants_transform as ants_transform_t1_mask_to_b0;
    ants_transform as transform_t1_mask_to_b0
} from '../processes/register.nf'
include { merge_channels_non_blocking; is_data } from '../functions.nf'
include {
    apply_mask as apply_mask_to_b0_for_reg;
    apply_mask as apply_mask_to_t1_for_reg;
    apply_mask as apply_mask_to_pdavg_for_reg;
    apply_mask as mask_target_b0;
    apply_mask as mask_target_pdavg;
    apply_mask as mask_moving_t1;
    prepend_sid as prepend_sid_template;
    bet_mask
} from '../processes/utils.nf'
include {
    get_data_path; get_config_path
} from "../functions.nf"

params.tissue_segmentation_root = "${get_data_path()}/maccaca_mulatta/tissue_segmentation"

params.ants_transform_base_config = file("${get_config_path()}/ants_transform_base_config.py")
params.ants_transform_mask_config = file("${get_config_path()}/ants_transform_mask_config.py")
params.t1_registration_extract_b0_config = file("${get_config_path()}/extract_mean_b0_base_config.py")
params.t1_to_template_registration_config = file("${get_config_path()}/t1_to_template_registration_config.py")
params.b0_to_template_registration_config = file("${get_config_path()}/b0_to_template_registration_config.py")
params.t1_to_template_registration_quick_config = file("${get_config_path()}/t1_to_template_registration_quick_config.py")
params.b0_to_template_registration_quick_config = file("${get_config_path()}/b0_to_template_registration_quick_config.py")
params.t1_to_b0_registration_config = file("${get_config_path()}/t1_to_b0_registration_affine.py")

workflow t12b0_registration {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        dwi_mask_channel
        publish_mask
        publish_t1
        use_quick
    main:
        extract_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)
        apply_mask_to_b0_for_reg(extract_b0.out.b0.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        apply_mask_to_t1_for_reg(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")
        compute_powder_average(dwi_channel.map{ it.subList(0, 3) }.map{ it + ["", ""] }, "preprocess", "false")
        apply_mask_to_pdavg_for_reg(compute_powder_average.out.image.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        scil_compute_dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess", false)

        target_channel = apply_mask_to_b0_for_reg.out.image
            .join(apply_mask_to_pdavg_for_reg.out.image)
            .join(scil_compute_dti_fa.out.fa)
            .map{ [it[0], it[1..-1]] }
        moving_channel = apply_mask_to_t1_for_reg.out.image
            .map{ [it[0], [it[1]]] }
        template_channel = prepend_sid_template(
            t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] }
        ).map{ [it[0], [it[1]]] }

        t1_to_template_registration_wkf(
            template_channel,
            moving_channel,
            t1_mask_channel,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            use_quick ? params.t1_to_template_registration_quick_config : params.t1_to_template_registration_config,
            params.ants_transform_mask_config
        )

        b0_to_template_registration_wkf(
            template_channel,
            target_channel,
            null,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            use_quick ? params.b0_to_template_registration_quick_config : params.b0_to_template_registration_config,
            params.ants_transform_mask_config
        )

        inverse_transform = b0_to_template_registration_wkf.out.inverse_transform
            .map{ [it[0], it[1], it[1].collect{ f -> f.name.contains("syn_inverse") ? "false": "true" }] }

        ants_transform_t1_to_b0(
            t1_to_template_registration_wkf.out.registration
                .join(extract_b0.out.b0)
                .join(inverse_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_base_config
        )
        ants_transform_t1_mask_to_b0(
            t1_to_template_registration_wkf.out.image
                .join(extract_b0.out.b0)
                .join(inverse_transform)
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_mask_config
        )

        transformation = t1_to_template_registration_wkf.out.transform
            .map{ [it[0], it[1], it[1].collect{ "false" }] }
            .join(inverse_transform)
            .map{ [it[0], it[1] + it[3], it[2] + it[4]] }

    emit:
        t1 = ants_transform_t1_to_b0.out.image
        mask = ants_transform_t1_mask_to_b0.out.image
        transform = transformation
        reference = extract_b0.out.b0
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

        compute_target_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess", false)
        mask_moving_t1(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")

        target_channel = mask_target_b0.out.image
            .join(mask_target_pdavg.out.image)
            .join(compute_target_fa.out.fa)
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
                .map{ it + ["", ""] },
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_mask_config
        )
    emit:
        mask = transform_t1_mask_to_b0.out.image
}