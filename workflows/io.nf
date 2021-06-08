#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_root = false
params.masked_dwi = false
params.msmt_odf = false
params.masked_t1 = true
params.rev_is_b0 = true

include { prepare_metadata as pmeta_dwi; prepare_metadata as pmeta_rev } from "../modules/processes/io.nf"

workflow load_dataset {
    main:
        if ( !params.data_root )
            error "You must supply an input data root using --data_root"
        root = file(params.data_root)
        dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true).map{ [it[0], it[3], it[1], it[2]] }
        dwi_meta_channel = Channel.fromFilePairs("$root/**/*dwi.json", size: 1, flat: true)
        anat_channel = Channel.fromFilePairs("$root/**/*t1.nii.gz", size: 1, flat: true)

        seg_channel = null
        if ( params.msmt_odf )
            seg_channel = Channel.fromFilePairs("$root/**/*{wm,gm,csf}_mask.nii.gz", size: 3, flat: true).map{ [it[0], [it[3], it[2], it[1]]] }
        rev_channel = null
        rev_meta_channel = null

        if ( params.has_reverse ) {
            if (params.rev_is_b0)
                rev_channel = Channel.fromFilePairs("$root/**/*rev.nii.gz", size: 1, flat: true)
            else {
                rev_channel = Channel.fromFilePairs("$root/**/*rev.{nii.gz,bval,bvec}", size: -1, flat: true).map {
                    (0..3).collect { i -> i >= it.size() ? null : it[i] }
                }.map {
                    [it[0], it[3], it[1], it[2]]
                }.join(dwi_channel.map {
                    [it[0]] + it.subList(2, it.size())
                }).map {
                    it[2] ? it.subList(0, 4) + [it[-1]] : it.subList(0, 2) + [it[4], it[3], it[5]]
                }.map {
                    it[3] ? it.subList(0, 4) : it.subList(0, 3) + [it[-1]]
                }
                rev_json_channel = Channel.fromFilePairs("$root/**/*rev.json", size: 1, flat: true)
                in_meta_rev = rev_channel.map { it.subList(0, 2) }.join(rev_json_channel, remainder: true).map { it.size() > 2 ? it[-1] ? it : [it[0], it[1], ""] : it + [""] }
                rev_meta_channel = pmeta_rev(in_meta_rev.map { it + ["true"] })
            }
        }

        dwi_mask_channel = null
        anat_mask_channel = null
        if ( params.masked_t1 )
            anat_mask_channel = Channel.fromFilePairs("$root/**/*t1_mask.nii.gz", size: 1, flat: true)
        if ( params.masked_dwi )
            dwi_mask_channel = Channel.fromFilePairs("$root/**/*dwi_mask.nii.gz", size: 1, flat: true)

        dwi_meta_channel = pmeta_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_meta_channel, remainder: true).map{ it.size() > 2 ? it[-1] ? it : [it[0], it[1], ""] : it + [""] }.map{ it + ["false"] })
    emit:
        dwi = dwi_channel
        dwi_mask = dwi_mask_channel
        anat = anat_channel
        anat_mask = anat_mask_channel
        rev = rev_channel
        seg = seg_channel
        metadata = dwi_meta_channel
        rev_metadata = rev_meta_channel
}
