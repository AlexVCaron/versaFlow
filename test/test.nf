#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_root = false

include { t12b0_registration } from "../modules/workflows/t1_registration.nf"

workflow {
    root = file(params.data_root)
    t1_channel = Channel.fromFilePairs("$root/**/*t1.nii.gz", size: 1, flat: true)
    t1_mask_channel = Channel.fromFilePairs("$root/**/*t1_mask.nii.gz", size: 1, flat: true)
    dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{bval,bvec,nii.gz}", size: 3, flat: true).map{ [it[0], it[3], it[1], it[2]] }
    
    t12b0_registration(dwi_channel, t1_channel, t1_mask_channel, null, null, true, true)
}