#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.b0_threshold = false

include { remove_alg_suffixes; add_suffix } from '../functions.nf'

process extract_b0 {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_b0") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(metadata)
        val(caller_name)
        val(publish)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}_b0.nii.gz"), emit: b0
        tuple val(sid), path("${dwi.simpleName}_b0*_metadata.*"), optional: true, emit: metadata
    script:
        def extra_args = params.b0_threshold ? "--ceil ${params.b0_threshold}" : ""
        """
        magic-monkey b0 extract --in $dwi --bvals $bval --out ${dwi.simpleName}_b0 --config $config $extra_args
        """
}

process squash_b0 {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(metadata)
        val(caller_name)
        val(publish)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__b0_squashed.nii.gz"), path("${dwi.simpleName}__b0_squashed.bval"), path("${dwi.simpleName}__b0_squashed.bvec"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__b0_squashed_metadata.*"), optional: true, emit: metadata
    script:
        def extra_args = params.b0_threshold ? "--ceil ${params.b0_threshold}" : ""
        """
        magic-monkey b0 squash --in $dwi --bvals $bval --bvecs $bvec --out ${dwi.simpleName}__b0_squashed --config $config $extra_args
        """
}
