#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { registration_wkf as nmt_registration_wkf; registration_wkf as wm_seg_registration_wkf } from "./preprocess.nf"
include { atropos } from '../processes/segment.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include { 
    prepend_sid as prepend_sid_template;
    prepend_sid as prepend_sid_template_mask;
    prepend_sid as prepend_sid_segmentation;
    prepend_sid as prepend_sid_template_fa;
    prepend_sid as prepend_sid_wm_atlas;
    pvf_to_mask
} from '../processes/utils.nf'
include {
    resampling_reference;
    resampling_reference as resampling_reference_fa;
    scilpy_resample_to_reference as resample_template;
    scilpy_resample_to_reference as resample_segmentation;
    scilpy_resample_to_reference as resample_template_fa;
    scilpy_resample_to_reference as resample_wm_atlas
} from '../processes/upsample.nf'
include { get_config_path; get_data_path } from '../functions.nf'   

params.tissue_segmentation_root = "${get_data_path()}/maccaca_mulatta/tissue_segmentation"
params.wm_segmentation_root = "${get_data_path()}/maccaca_mulatta/wm_segmentation"

params.segmentation_registration_config = file("${get_config_path()}/segmentation_registration_config.py")
params.ants_transform_segmentation_config = file("${get_config_path()}/ants_transform_segmentation_config.py")


workflow segment_nmt_wkf {
    take:
        t1_channel
        mask_channel
    main:
        template_channel = prepend_sid_template(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] })
        template_mask_channel = prepend_sid_template_mask(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_mask.nii.gz")] })
        segmentation_channel = prepend_sid_segmentation(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation.nii.gz")] })

        resampling_reference(t1_channel.join(template_channel).map{ [it[0], it[1..-1]] }, "segmentation")
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
        resample_segmentation(
            segmentation_channel
                .join(resampling_reference.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "nn",
            false,
            false,
            "", ""
        )

        nmt_registration_wkf(
            t1_channel.map{ [it[0], [it[1]]] },
            resample_template.out.image.map{ [it[0], [it[1]]] },
            resample_segmentation.out.image.map{ [it[0], [it[1]]] },
            mask_channel.join(resample_template.out.mask).map{ [it[0], [it[1]]] },
            null,
            null,
            "segmentation",
            false,
            "", "",
            params.segmentation_registration_config,
            params.ants_transform_segmentation_config
        )

        atropos(t1_channel.join(mask_channel).join(nmt_registration_wkf.out.image), "segment")
        pvf_to_mask(atropos.out.vol_fractions.map{ [it[0]] + it[1].reverse() }.join(mask_channel), "segment", "segmentation")
    emit:
        segmentation = atropos.out.segmentation
        volume_fractions = atropos.out.vol_fractions
        tissue_masks = pvf_to_mask.out.wm_mask.join(pvf_to_mask.out.gm_mask).join(pvf_to_mask.out.csf_mask).map{ [ it[0], it[1..-1] ]}
        safe_wm_mask = pvf_to_mask.out.safe_wm_mask
}

workflow segment_wm_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        template_fa_channel = prepend_sid_template_fa(dwi_channel.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_fa.nii.gz")] })
        wm_atlas_channel = prepend_sid_wm_atlas(dwi_channel.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_atlas.nii.gz")] })

        resampling_reference_fa(dwi_channel.map{ it.subList(0, 2) }.join(template_fa_channel).map{ [it[0], it[1..-1]] }, "segmentation")
        resample_template_fa(
            template_fa_channel
                .join(resampling_reference_fa.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "lin",
            false,
            false,
            "", ""
        )
        resample_wm_atlas(
            wm_atlas_channel
                .join(resampling_reference_fa.out.reference)
                .map{ it + ["", ""] },
            "preprocess",
            "nn",
            false,
            false,
            "", ""
        )

        scil_compute_dti_fa(dwi_channel.join(mask_channel), "segmentation", "segmentation", false)

        wm_seg_registration_wkf(
            scil_compute_dti_fa.out.fa.map{ [it[0], [it[1]]] },
            resample_template_fa.out.image.map{ [it[0], [it[1]]] },
            resample_wm_atlas.out.image.map{ [it[0], [it[1]]] },
            null,
            null,
            null,
            "segmentation",
            true,
            "", "",
            params.segmentation_registration_config,
            params.ants_transform_segmentation_config
        )
    emit:
        segmentation = wm_seg_registration_wkf.out.image
}
