#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.nmt_root = "$projectDir/.data/segmentation"
params.nmt_ver = "2.0_asym"
params.wm_seg_root = "$projectDir/.data/wm_segmentation"

params.segmentation_registration_config = file("$projectDir/.config/segmentation_registration_config.py")

include { registration_wkf as nmt_registration_wkf; registration_wkf as wm_seg_registration_wkf } from "./preprocess.nf"
include { atropos } from '../processes/segment.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'

workflow segment_nmt_wkf {
    take:
        t1_channel
        mask_channel
    main:
        nmt_registration_wkf(
            t1_channel.map{ [it[0], [it[1]]] },
            t1_channel.map{ [it[0], [file("${params.nmt_root}/NMT_v${params.nmt_ver}_SS.nii.gz")]] },
            t1_channel.map{ [it[0], [file("${params.nmt_root}/NMT_v${params.nmt_ver}_segmentation.nii.gz")]] },
            mask_channel.map{ [it[0], [it[1], file("${params.nmt_root}/NMT_v${params.nmt_ver}_brainmask.nii.gz")]] },
            null,
            null,
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
            scil_compute_dti_fa.out.fa.map{ [it[0], [file("${params.wm_seg_root}/fa.nii.gz")]] },
            scil_compute_dti_fa.out.fa.map{ [it[0], [file("${params.wm_seg_root}/wm_atlas.nii.gz")]] },
            null,
            null,
            null,
            params.segmentation_registration_config
        )
    emit:
        segmentation = wm_seg_registration_wkf.out.image
}