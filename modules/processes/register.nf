#!/usr/bin/env nextflow

import java.io.File

nextflow.enable.dsl=2


include { get_size_in_gb; swap_configurations } from '../functions.nf'

process ants_register {
    memory { 4f * (get_size_in_gb(moving) + get_size_in_gb(target)) }
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(moving), path(target), val(reference), file(mask), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${moving[0].simpleName}__registration_affine.mat"), path("${moving[0].simpleName}__registration_rigid.nii.gz"), optional: true, emit: affine
        tuple val(sid), path("${moving[0].simpleName}__registration_ref.nii.gz"), emit: reference
        tuple val(sid), path("${moving[0].simpleName}__registration_warped.nii.gz"), optional: true, emit: image
        tuple val(sid), path("${moving[0].simpleName}__registration*syn.nii.gz"), optional: true, emit: syn
        tuple val(sid), path("${moving[0].simpleName}__registration_warped_metadata.*"), optional: true, emit: metadata
    script:
        def mask_arg = ""
        if ( !mask.iterator().inject(false) { c, i -> c || i.empty() } ) {
            mask_arg = "--mask ${mask.iterator().collect{ it.name }.join(',')}"
        }

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        magic-monkey ants_registration --moving ${moving.join(",")} --target ${target.join(",")} --out ${moving[0].simpleName}__registration $mask_arg --config $config
        cp ${file(reference).name} ${moving[0].simpleName}__registration_ref.nii.gz
        cp ${moving[0].simpleName}__registration_warped.nii.gz ${moving[0].simpleName}__registration_rigid.nii.gz
        if [ -f ${moving[0].simpleName}__registration0GenericAffine.mat ]
        then
            mv ${moving[0].simpleName}__registration0GenericAffine.mat ${moving[0].simpleName}__registration_affine.mat
        fi
        if [ -f "${moving[0].simpleName}__registration1Warp.nii.gz" ]
        then
            mv ${moving[0].simpleName}__registration1Warp.nii.gz ${moving[0].simpleName}__registration_syn.nii.gz
            mv ${moving[0].simpleName}__registration1InverseWarp.nii.gz ${moving[0].simpleName}__registration_inv_syn.nii.gz
        fi
        """
}

process ants_correct_motion {
    memory { 4f * (get_size_in_gb(moving) + get_size_in_gb(target)) }
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(moving), path(target), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}__motion_correct_warped.nii.gz"), emit: image
        tuple val(sid), path("${sid}__motion_correct_warped_metadata.*"), optional: true, emit: metadata
    script:
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        magic-monkey ants_motion --moving ${moving.join(",")} --target ${target.join(",")} --out ${sid}__motion_correct --config $config
        """
}

process ants_transform {
    memory { 4f * get_size_in_gb([img, ref]) }
    label "res_single_cpu"

    publishDir "${params.output_root}/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(img), path(ref), path(affine), file(trans), file(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${img.simpleName}__transformed.nii.gz"), emit: image
        tuple val(sid), path("${img.simpleName}__transformed_metadata.*"), optional: true, emit: metadata
    script:
        args = "--in $img --ref $ref --mat $affine"
        if ( trans && !trans.empty() ) {
            args += " --trans ${trans}"
        }
        """
        magic-monkey ants_transform $args --out ${img.simpleName}__transformed.nii.gz --config $config
        """
}
