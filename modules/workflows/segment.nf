#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.tissue_segmentation_root = "$projectDir/.data/segmentation"
params.wm_segmentation_root = "$projectDir/.data/wm_segmentation"

params.segmentation_registration_config = file("$projectDir/.config/segmentation_registration_config.py")

include { registration_wkf as nmt_registration_wkf; registration_wkf as wm_seg_registration_wkf } from "./preprocess.nf"
include { atropos } from '../processes/segment.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include { prepend_sid as prepend_sid_template; prepend_sid as prepend_sid_segmentation; prepend_sid as prepend_sid_template_fa; prepend_sid as prepend_sid_wm_atlas } from '../processes/utils.nf'

workflow segment_nmt_wkf {
    take:
        t1_channel
        mask_channel
    main:
        nmt_registration_wkf(
            t1_channel.map{ [it[0], [it[1]]] },
            prepend_sid_template(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] }).map{ [it[0], [it[1]]] },
            prepend_sid_segmentation(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation.nii.gz")] }).map{ [it[0], [it[1]]] },
            mask_channel.map{ [it[0], [it[1], file("${params.tissue_segmentation_root}/tissue_segmentation_mask.nii.gz")]] },
            null,
            null,
            "segmentation",
            false,
            "", "",
            params.segmentation_registration_config
        )
        atropos(t1_channel.join(mask_channel).join(nmt_registration_wkf.out.image), "segment")
    emit:
        segmentation = atropos.out.segmentation
        masks = atropos.out.masks
}

workflow segment_wm_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        scil_compute_dti_fa(dwi_channel.join(mask_channel), "segment", "segment")
        wm_seg_registration_wkf(
            scil_compute_dti_fa.out.fa.map{ [it[0], [it[1]]] },
            prepend_sid_template_fa(scil_compute_dti_fa.out.fa.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_fa.nii.gz")] }).map{ [it[0], [it[1]]] },
            prepend_sid_wm_atlas(scil_compute_dti_fa.out.fa.map{ [it[0], file("${params.wm_segmentation_root}/wm_segmentation_atlas.nii.gz")] }).map{ [it[0], [it[1]]] },
            null,
            null,
            null,
            "segmentation",
            true,
            "", "",
            params.segmentation_registration_config
        )
    emit:
        segmentation = wm_seg_registration_wkf.out.image
}
