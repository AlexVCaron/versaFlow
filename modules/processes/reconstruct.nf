#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.frf_fa = 0.7
params.frf_min_fa = 0.5
params.frf_min_nvox = 300
params.frf_radii = false
params.frf_center = false


include { get_size_in_gb; uniformize_naming } from '../functions.nf'

process diamond {
    label params.on_hcp ? "res_full_node_override" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/diamond", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(input_dwi), path(mask), path(data)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_diamond*.nii.gz"), emit: diamond
    script:
        if ( "${mask}" != "" )
            args += " --mask $mask"

        """
        magic-monkey flip2ref --in $input_dwi --bvecs ${input_dwi.simpleName}.bvec --out flipped_bvecs
        cp $input_dwi flipped_bvecs.nii.gz
        cp ${input_dwi.simpleName}.bval flipped_bvecs.bval
        magic-monkey diamond --in flipped_bvecs.nii.gz --mask $mask --out ${sid}_diamond --config $config
        """
}

process mrtrix_dti {
    label "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(mask)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_dti_dti.nii.gz"), emit: dti
    script:
        args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" != "" )
            args += " --mask $mask"

        """
        magic-monkey dti $args --out ${sid}_dti --config $config
        """
}

process response {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(mask)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_response_*.txt"), emit: responses
    script:
        args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" != "" )
            args += " --mask $mask"

        """
        magic-monkey response $args --out ${sid}_response --config $config
        """
}

process csd {
    label "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(responses), path(dwi), path(bval), path(bvec), path(mask)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_csd_*.nii.gz"), emit: odfs
    script:
        args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" == "" )
            args += " --mask $mask"

        """
        magic-monkey csd $args --out ${sid}_csd --responses ${responses.join(',')} --config $config
        """
}

process scilpy_response {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(mask)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_response.txt"), emit: response
    script:
        args = ""
        if (params.frf_radii)
            args += " --roi_radii $params.frf_radii"
        if (params.frf_center)
            args += " --roi_center ${params.frf_center.join(" ")}"
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrconvert -datatype uint8 $mask mask4scil.nii.gz
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        magic-monkey shells --in $dwi --bvals $bval --bvecs flipped_bvecs.bvec --shells 1500 --keep leq --out dwi_leq_1500 --with_b0
        scil_compute_ssst_frf.py dwi_leq_1500.nii.gz dwi_leq_1500.bval dwi_leq_1500.bvec ${sid}_response.txt --mask mask4scil.nii.gz --fa $params.frf_fa --min_fa $params.frf_min_fa --min_nvox $params.frf_min_nvox $args
        """
}

process scilpy_msmt_response {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
    tuple val(sid), path(dwi), path(bval), path(bvec), path(mask), path(seg)
    val(caller_name)
    output:
    tuple val(sid), path("${sid}_wm_response.txt"), path("${sid}_gm_response.txt"), path("${sid}_csf_response.txt"), emit: response
    script:
    args = ""
    if (params.frf_radii)
        args += " --roi_radii $params.frf_radii"
    if (params.frf_center)
        args += " --roi_center ${params.frf_center.join(" ")}"
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrconvert -datatype uint8 $mask mask4scil.nii.gz
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        scil_compute_msmt_frf.py $dwi $bval flipped_bvecs.bvec ${sid}_wm_response.txt ${sid}_gm_response.txt ${sid}_csf_response.txt --mask mask4scil.nii.gz --mask_wm ${seg[0]} --mask_gm ${seg[1]} --mask_csf ${seg[2]} --fa_thr_wm $params.frf_fa --min_nvox $params.frf_min_nvox $args
        """
}

process scilpy_csd {
    label "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(response), path(mask)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_fodf.nii.gz"), emit: odfs
    script:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrconvert -datatype uint8 $mask mask4scil.nii.gz
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        scil_compute_ssst_fodf.py $dwi $bval flipped_bvecs.bvec $response ${sid}_fodf.nii.gz --mask mask4scil.nii.gz --force_b0_threshold --sh_order $params.sh_order --processes $task.cpus
        """
}

process scilpy_msmt_csd {
    label "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
    tuple val(sid), path(dwi), path(bval), path(bvec), path(wm_response), path(gm_response), path(csf_response), path(mask)
    val(caller_name)
    output:
    tuple val(sid), path("${sid}_wm_fodf.nii.gz"), path("${sid}_gm_fodf.nii.gz"), path("${sid}_csf_fodf.nii.gz"), emit: odfs
    tuple val(sid), path("${sid}_vf.nii.gz"), emit: vf
    script:
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrconvert -datatype uint8 $mask mask4scil.nii.gz
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        scil_compute_msmt_fodf.py $dwi $bval flipped_bvecs.bvec $wm_response $gm_response $csf_response --wm_out_fODF ${sid}_wm_fodf.nii.gz --gm_out_fODF ${sid}_gm_fodf.nii.gz --csf_out_fODF ${sid}_csf_fodf.nii.gz --vf ${sid}_vf.nii.gz --mask mask4scil.nii.gz --force_b0_threshold --sh_order $params.sh_order --processes $task.cpus
        """
}