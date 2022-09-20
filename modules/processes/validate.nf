#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.b0_threshold = false

process validate_affine {
    label "res_single_cpu"

    input:
        tuple val(sid), path(ref_image), path(cmp_image)
    output:
        tuple val(sid), env(ERROR_MSG), emit: errors
    script:
        """
        ERROR_MSG=$(mrHARDI compare affine --ref $ref_image --in $cmp_image --stdout)
        """
}

process validate_dwi_acquisition {
    label "res_single_cpu"

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), env(ERROR_MSG), emit: errors
    script:
        def args = ""
        if (params.b0_threshold)
            args += " --b0-thr $params.b0_threshold"
        """
        ERROR_MSG=$(mrHARDI validate_dwi --in $dwi --bvals $bval --bvecs $bvec $args)
        """
}