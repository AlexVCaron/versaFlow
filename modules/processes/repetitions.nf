#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_size_in_gb; remove_alg_suffixes } from '../functions.nf'


process ants_register_dwi_repetition {
    memory { 4f * (get_size_in_gb(target_b0) + get_size_in_gb(dwi)) }
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(target_b0), val(rep_idx), path(dwi), path(bval), path(bvec), path(metadata)
        val(caller_name)
        path(b0_config)
        path(reg_config)
        path(trans_config)
    output:
        tuple val("${sid}_${rep_idx}"), path("${dwi.simpleName}__rep_registered.nii.gz"), path("${dwi.simpleName}__rep_registered.bval"), path("${dwi.simpleName}__rep_registered.bvec"), emit: dwi
        tuple val("${sid}_${rep_idx}"), path("${dwi.simpleName}__rep_registered_metadata.*"), emit: metadata
    script:
    """
    export OMP_NUM_THREADS=$task.cpus
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OPENBLAS_NUM_THREADS=1
    magic-monkey b0 extract --in $dwi --bvals $bval --out rep_b0 --config $b0_config
    magic-monkey ants_registration --target ${target_b0} --moving rep_b0.nii.gz --out b0_registration --config $reg_config
    magic-monkey ants_transform --in $dwi --ref ${target_b0} --trans b0_registration0GenericAffine.mat --bvecs $bvec --out ${dwi.simpleName}__rep_registered --config $trans_config
    cp $bval ${dwi.simpleName}__rep_registered.bval
    """
}

process ants_register_t1_repetition {
    memory { 4f * (get_size_in_gb(ref_t1) + get_size_in_gb(t1)) }
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(ref_t1), val(rep_idx), path(t1), file(mask)
        val(caller_name)
        file(reg_config)
    output:
        tuple val("${sid}_${rep_idx}"), path("${t1.simpleName}__rep_registered.nii.gz"), emit: t1
        tuple val("${sid}_${rep_idx}"), path("${mask.simpleName}__rep_registered.nii.gz"), emit: mask, optional: true
    script:
        def command = ""
        if (!mask.empty()) {
            command += "magic-monkey ants_transform --in $mask --ref $ref_t1 --trans t1_registration0GenericAffine.mat --out ${mask.simpleName}__rep_registered"
        }
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        magic-monkey ants_registration --target $ref_t1 --moving $t1 --out t1_registration --config $reg_config
        magic-monkey ants_transform --in $t1 --ref $ref_t1 --trans t1_registration0GenericAffine.mat --out ${t1.simpleName}__rep_registered
        $command
        """
}