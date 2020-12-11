#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_root = false
params.masked_dwi = false
params.masked_t1 = true
params.rev_is_b0 = true

include { key_from_filename } from "../modules/functions.nf"
include { prepare_metadata as pmeta_dwi; prepare_metadata as pmeta_rev } from "../modules/processes/io.nf"

workflow load_dataset {
    main:
        if ( !params.data_root )
            error "You must supply an input data root using --data_root"
        root = file(params.data_root)
        dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true).map{ [it[0], it[3], it[1], it[2]] }
        dwi_meta_channel = key_from_filename(Channel.fromPath("$root/**/*dwi.json"), "_").map{ [it[0].substring(0, it[0].lastIndexOf("_"))] + it.subList(1, it.size()) }
        affine_channel = key_from_filename(Channel.fromPath("$root/**/*.affine"), ".")
        anat_channel = key_from_filename(Channel.fromPath("$root/**/*t1.nii.gz"), "_")
        rev_channel = null
        rev_meta_channel = null

        if ( params.rev_is_b0 )
            rev_channel = key_from_filename(Channel.fromPath("$root/**/*rev.nii.gz"), "_")
        else {
            rev_channel = Channel.fromFilePairs("$root/**/*rev.{nii.gz,bval,bvec}", size: -1, flat: true).map{
                (0..3).collect{ i -> i >= it.size() ? null : it[i] }
            }.map{
                [it[0], it[3], it[1], it[2]]
            }.join(dwi_channel.map{
                [it[0]] + it.subList(2, it.size())
            }).map{
                it[2] ? it.subList(0, 4) + [it[-1]] : it.subList(0, 2) + [it[4], it[3], it[5]]
            }.map{
                it[3] ? it.subList(0, 4) : it.subList(0, 3) + [it[-1]]
            }
            rev_json_channel = key_from_filename(Channel.fromPath("$root/**/*rev.json"), "_").map{ [it[0].substring(0, it[0].lastIndexOf("_"))] + it.subList(1, it.size()) }
            in_meta_rev = rev_channel.map{ it.subList(0, 2) }.join(rev_json_channel, remainder: true).map{ it.size() > 2 ? it[-1] ? it : [it[0], it[1], ""] : it + [""] }
            rev_meta_channel = pmeta_rev(in_meta_rev.map{ it + ["true"] })
        }

        if ( params.masked_t1 )
            anat_channel = anat_channel.join(key_from_filename(Channel.fromPath("$root/**/*t1_mask.nii.gz"), "_"))
        if ( params.masked_dwi )
            dwi_channel = dwi_channel.join(key_from_filename(Channel.fromPath("$root/**/*dwi_mask.nii.gz"), "_"))

        dwi_meta_channel = pmeta_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_meta_channel, remainder: true).map{ it.size() > 2 ? it[-1] ? it : [it[0], it[1], ""] : it + [""] }.map{ it + ["false"] })
    emit:
        dwi = dwi_channel
        affine = affine_channel
        anat = anat_channel
        rev = rev_channel
        metadata = dwi_meta_channel
        rev_metadata = rev_meta_channel
}
