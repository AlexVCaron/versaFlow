#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.nmt_root = "/nmt"
params.nmt_ver = "2.0_asym"

params.segmentation_registration_config = file("$projectDir/.config/segmentation_registration_config.py")

include { registration_wkf } from "./preprocess.nf"
include { atropos } from '../processes/segment.nf'

workflow segment_nmt {
    take:
        t1_channel
        mask_channel
    main:
        registration_wkf(
            t1_channel.map{ [it[0], [it[1]]] },
            t1_channel.map{ [it[0], [file("${params.nmt_root}/NMT_v${params.nmt_ver}_SS.nii.gz")]] },
            t1_channel.map{ [it[0], [file("${params.nmt_root}/NMT_v${params.nmt_ver}_segmentation.nii.gz")]] },
            mask_channel.map{ it + [file("${params.nmt_root}/NMT_v${params.nmt_ver}_brainmask.nii.gz")] },
            null,
            null,
            params.segmentation_registration_config
        )
        atropos(t1_channel.join(mask_channel).join(registration_wkf.out.image), "preprocess")
    emit:
        segmentation = atropos.out.segmentation
        masks = atropos.out.masks
}
