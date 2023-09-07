#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.b0_threshold = false

include { remove_alg_suffixes; add_suffix } from '../functions.nf'

process compute_powder_average {
    label "FAST"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_b0") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), file(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${dwi.simpleName}_pd_avg.nii.gz"), emit: image
        tuple val(sid), path("${dwi.simpleName}_pd_avg_metadata.*"), optional: true, emit: metadata
    script:
        def after_script = ""
        if (!metadata.empty()) {
            after_script += "cp $metadata ${dwi.simpleName}_pd_avg_metadata.py"
        }
        """
        scil_compute_powder_average.py $dwi $bval ${dwi.simpleName}_pd_avg.nii.gz ${mask.empty() ? "" : "--mask $mask"}
        $after_script
        """

}

process extract_b0 {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_b0") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), file(metadata)
        val(caller_name)
        val(publish)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}_b0.nii.gz"), emit: b0
        tuple val(sid), path("${dwi.simpleName}_b0*_metadata.*"), optional: true, emit: metadata
    script:
        def extra_args = params.b0_threshold ? "--ceil ${params.b0_threshold}" : ""
        """
        mrhardi b0 extract --in $dwi --bvals $bval --out ${dwi.simpleName}_b0 --config $config $extra_args
        """
}

process squash_b0 {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(metadata)
        val(caller_name)
        val(publish)
        path(config)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${dwi.simpleName}__b0_squashed.nii.gz"), path("${dwi.simpleName}__b0_squashed.bval"), path("${dwi.simpleName}__b0_squashed.bvec"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__b0_squashed_metadata.*"), optional: true, emit: metadata
    script:
        def extra_args = params.b0_threshold ? "--ceil ${params.b0_threshold}" : ""
        """
        mrhardi b0 squash --in $dwi --bvals $bval --bvecs $bvec --out ${dwi.simpleName}__b0_squashed --config $config $extra_args
        """
}
