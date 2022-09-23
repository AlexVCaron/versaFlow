#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_root = false

include { prepare_metadata as pmeta_dwi; prepare_metadata as pmeta_rev } from "../modules/processes/io.nf"
include { fill_missing_datapoints; exclude_missing_datapoints; separate_b0_from_dwi } from "../modules/functions.nf"
include {
    validate_affine as anat_validate_affine;
    validate_affine as dwi_validate_affine;
    validate_dwi_acquisition
} from "../modules/processes/validate.nf"

workflow load_dataset {
    main:
        if ( !params.data_root )
            error "You must supply an input data root using --data_root"
        root = file(params.data_root)

        // Load DWI and T1, those datapoints are all required for all subject/session
        dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true).map{ [it[0], it[3], it[1], it[2]] }
        anat_channel = Channel.fromFilePairs("$root/**/*t1.nii.gz", size: 1, flat: true)
        ref_id_channel = anat_channel.map{ [it[0]] }

        rev_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*rev.nii.gz", size: 1, flat: true),
            ref_id_channel,
            1, [""]
        )
        rev_bval_bvec = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*rev.{bval,bvec}", size: 2, flat: true),
            ref_id_channel,
            1, ["", ""]
        )
        rev_channel = rev_channel.join(rev_bval_bvec)

        // Load WM/GM/CSF segmentation if present
        pvf_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*{wm,gm,csf}_pvf.nii.gz", size: 3, flat: true),
            ref_id_channel,
            1, []
        )

        pvf_channel = pvf_channel.map{ [it[0], it.subList(1, it.size()).reverse()] }

        // Load per subject/session DWI json metadata specification and transform
        dwi_json_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*dwi.json", size: 1, flat: true),
            ref_id_channel,
            1, [""]
        )
        rev_json_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*rev.json", size: 1, flat: true),
            ref_id_channel,
            1, [""]
        )
        dwi_meta_channel = pmeta_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_json_channel).map{ it + ["false"] })
        rev_meta_channel = fill_missing_datapoints(
            pmeta_rev(exclude_missing_datapoints(rev_channel.map{ it.subList(0, 2) }.join(rev_json_channel), 1, "").map{ it + ["true"] }),
            ref_id_channel,
            1, [""]
        )

        // Load available masks (T1 and/or DWI)
        dwi_mask_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*dwi_mask.nii.gz", size: 1, flat: true),
            ref_id_channel,
            1, [""]
        )
        anat_mask_channel = fill_missing_datapoints(
            Channel.fromFilePairs("$root/**/*t1_mask.nii.gz", size: 1, flat: true),
            ref_id_channel,
            1, [""]
        )

    emit:
        dwi = dwi_channel
        dwi_mask = dwi_mask_channel
        anat = anat_channel
        anat_mask = anat_mask_channel
        rev = rev_channel
        pvf = pvf_channel
        metadata = dwi_meta_channel
        rev_metadata = rev_meta_channel
}

workflow anat_validation_wkf {
    take:
        anat_channel
        mask_channel
        other_channel
    main:
        compare_channel = mask_channel.mix(other_channel.transpose())
        compare_channel = exclude_missing_datapoints(compare_channel, 1, "")
        anat_validate_affine(anat_channel.combine(compare_channel, by: 0))

        errors = anat_validate_affine.out.bad_affine
        error_channel = anat_validate_affine.out.errors
    emit:
        error_count = errors.count()
        errors = error_channel
}

workflow dwi_validation_wkf {
    take:
        dwi_channel
        rev_channel
        mask_channel
    main:
        compare_channel = mask_channel.mix(rev_channel.map{ it[0..1] })
        compare_channel = exclude_missing_datapoints(compare_channel, 1, "")
        dwi_validate_affine(dwi_channel.map{ it[0..1] }.combine(compare_channel, by: 0))
        (rev_dwi_channel, rev_b0_channel) = separate_b0_from_dwi(rev_channel)
        validate_dwi_acquisition(dwi_channel.mix(rev_dwi_channel))

        errors = dwi_validate_affine.out.bad_affine.mix(validate_dwi_acquisition.out.error_files)
        error_channel = dwi_validate_affine.out.errors.mix(validate_dwi_acquisition.out.errors)
    emit:
        error_count = errors.count()
        errors = error_channel
}
