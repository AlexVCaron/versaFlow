#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { extract_b0 as extract_rep_b0 } from '../processes/preprocess.nf'
include { merge_repetitions; get_config_path } from '../functions.nf'
include { ants_register_dwi_repetition; ants_register_dwi_repetition as ants_register_rev_repetition; ants_register_t1_repetition } from '../processes/repetitions.nf'
include { cat_datasets as cat_repetitions } from '../processes/utils.nf'

params.reps_registration_extract_b0_config = file("${get_config_path()}/extract_mean_b0_base_config.py")
params.reps_registration_b0_registration_config = file("${get_config_path()}/reps_registration_b0_registration_config.py")
params.reps_registration_apply_registration_config = file("${get_config_path()}/reps_registration_apply_registration_config.py")
params.reps_registration_t1_registration_config = file("${get_config_path()}/reps_registration_t1_registration_config.py")
params.concatenate_base_config = file("${get_config_path()}/concatenate_base_config.py")


workflow register_dwi_repetitions_wkf {
    take:
        dwi_channel
        rev_channel
        metadata_channel
        rev_metadata_channel
        additional_publish_path
    main:
        dwi_channel = merge_repetitions(dwi_channel, true)
        if ( is_data(rev_channel) ) {
            rev_channel = merge_repetitions(rev_channel, true).transpose()
            rev_metadata_channel = merge_repetitions(rev_metadata_channel, false)
        }
        metadata_channel = merge_repetitions(metadata_channel, true)

        extract_rep_b0(dwi_channel.map{ [it[0], it[2][0], it[3][0]] }.join(metadata_channel.map{ [it[0]] + it[2..-1] }), "preprocess", params.reps_registration_extract_b0_config)
        main_b0 = extract_rep_b0.out.b0

        dwi_reg = dwi_channel.map{ [it[0]] + it[1..-1].collect{ i -> i[1..-1] } }.transpose()

        ants_register_dwi_repetition(
            main_b0.combine(dwi_reg, by: 0).combine(metadata_channel.map{ [it[0]] + it[2..-1] }, by: 0),
            "preprocess",
            additional_publish_path,
            params.reps_registration_extract_b0_config,
            params.reps_registration_b0_registration_config,
            params.reps_registration_apply_registration_config
        )
        if ( is_data(rev_channel) ) {
            ants_register_rev_repetition(
                main_b0.combine(rev_channel, by: 0).combine(rev_metadata_channel, by: 0),
                "preprocess",
                additional_publish_path,
                params.reps_registration_extract_b0_config,
                params.reps_registration_b0_registration_config,
                params.reps_registration_apply_registration_config
            )
        }
    emit:
        dwi = ants_register_dwi_repetition.out.dwi.concat(dwi_channel.map{ ["${it[0]}_${it[1][0]}"] + it[2..-1].collect{ i -> i[0] } })
        rev = rev_channel ? ants_register_rev_repetition.out.dwi : null
        metadata = ants_register_dwi_repetition.out.metadata.mix(metadata_channel.map{ ["${it[0]}_${it[1][0]}", it[2][0]] })
        rev_metadata = rev_channel ? ants_register_rev_repetition.out.metadata : null
}

workflow register_t1_repetitions_wkf {
    take:
        t1_channel
        mask_channel
    main:
        t1_channel = mask_channel ? t1_channel.join(mask_channel) : t1_channel.map{ it + [""] }
        t1_channel = merge_repetitions(t1_channel, true)
        template_t1 = t1_channel.map{ [it[0], it[2][0]] }
        t1_reg = t1_channel.map{ [it[0]] + it[1..-1].collect{ i -> i[1..-1] } }.transpose()
        ants_register_t1_repetition(template_t1.combine(t1_reg, by: 0), "preprocess", additional_publish_path, params.reps_registration_t1_registration_config)
    emit:
        t1 = ants_register_t1_repetition.out.t1.concat(t1_channel.map{ ["${it[0]}_${it[1][0]}", it[2][0]] })
        mask = mask_channel ? ants_register_t1_repetition.out.mask.concat(t1_channel.map{ ["${it[0]}_${it[1][0]}", it[3][0]] }) : null
}

workflow cat_dwi_repetitions_wkf {
    take:
        dwi_channel
        metadata_channel
        prefix
    main:
        dwi_channel = merge_repetitions(dwi_channel, false)
        metadata_channel = merge_repetitions(metadata_channel, false).map{ [it[0], it[1..-1].inject([]){ c, t -> c + t }] }
        cat_repetitions(dwi_channel.join(metadata_channel), 3, prefix, "preprocess", params.concatenate_base_config)
    emit:
        dwi = cat_repetitions.out.image.join(cat_repetitions.out.bval).join(cat_repetitions.out.bvec)
        metadata = cat_repetitions.out.metadata
}
