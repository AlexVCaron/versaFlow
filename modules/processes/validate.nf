#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.b0_threshold = false

process validate_affine {
    label "res_single_cpu"

    input:
        tuple val(sid), path(ref_image), path(cmp_image)
    output:
        tuple val(sid), env(ERROR_MSG), emit: errors
        tuple val(sid), path(ref_image), emit: ref
        tuple val(sid), path("${sid}_validation_correct_affine.txt"), emit: valid_affine, optional: true
        tuple val(sid), path("${sid}_validation_affine_components.txt"), emit: bad_affine, optional: true
    script:
        """
        ERROR_MSG="\$(mrhardi validate affine --ref $ref_image --in $cmp_image --stdout --out ${sid}_validation)"
        """
}

process validate_dwi_acquisition {
    label "res_single_cpu"

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), env(ERROR_MSG), emit: errors
        tuple val(sid), path(dwi), path(bval), path(bvec), emit: dwi
        tuple val(sid), path("${sid}_validation_errors.txt"), emit: error_files, optional: true
    script:
        def args = ""
        if (params.b0_threshold)
            args += " --b0-thr $params.b0_threshold"
        """
        ERROR_MSG="\$(mrhardi validate dwi --in $dwi --bvals $bval --bvecs $bvec --stdout $args --out ${sid}_validation)"
        """
}