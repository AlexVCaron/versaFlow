#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.frf_fa = 0.7
params.frf_min_fa = 0.5
params.frf_min_nvox = 300
params.frf_radii = false
params.frf_center = false
params.n_fascicles = 3
params.fascicle_model = "diamondNCcyl"
params.model_selection_with_tensor = false
params.estimate_restriction = false
params.estimate_hindered = false
params.restriction_tensor = false
params.free_water_tensor = false
params.normalized_fractions = true
params.strict_parameters = true
params.frf_on_dti_shell = false
params.max_dti_bvalue = 1300
params.random_seed = 1234

process diamond {
    label params.on_hcp ? "res_full_node_override" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/diamond", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(input_dwi), file(mask), path(data)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_diamond*.nii.gz"), emit: diamond
        tuple val(sid), path("${sid}_diamond.xml"), emit: xml_summary
    script:
        def args = ""
        if ( !mask.empty() )
            args += " --mask $mask"
        if ( params.model_selection_with_tensor )
            args += " --mose-tensor"
        if ( params.estimate_restriction )
            args += " --restricted"
        if ( params.restriction_tensor )
            args += " --res-tensor"
        if ( !params.normalized_fractions )
            args += " --nosum-fractions"
        if ( params.free_water_tensor )
            args += " --iso-tensor"
        if ( params.estimate_hindered )
            args += " --hindered"
        if ( !params.strict_parameters )
            args += " --lenient-params"

        """
        magic-monkey diamond --in $input_dwi --out ${sid}_diamond --n $params.n_fascicles --f $params.fascicle_model --p $task.cpus --config $config $args
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
        def args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" != "" )
            args += " --mask $mask"

        """
        export MRTRIX_RNG_SEED=$params.random_seed
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
        def args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" != "" )
            args += " --mask $mask"

        """
        export MRTRIX_RNG_SEED=$params.random_seed
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
        def args = "--in $dwi --bvals $bval --bvecs $bvec"
        if ( "${mask}" == "" )
            args += " --mask $mask"

        """
        export MRTRIX_RNG_SEED=$params.random_seed
        magic-monkey csd $args --out ${sid}_csd --responses ${responses.join(',')} --config $config
        """
}

process scilpy_response {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(mask), file(wm_mask)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_response.txt"), emit: response
    script:
        def args = ""
        def before_frf = ""
        if (params.frf_radii)
            args += " --roi_radii $params.frf_radii"
        if (params.frf_center)
            args += " --roi_center ${params.frf_center.join(" ")}"
        if (params.frf_on_dti_shell && params.max_dti_bvalue)
            before_frf += "magic-monkey shells --in $dwi --bvals $bval --bvecs $bvec --shells $params.max_dti_bvalue --keep leq --out dwi_frf_shells --with_b0\n"
        else {
            before_frf += "cp $dwi dwi_frf_shells.nii.gz\n"
            before_frf += "cp $bval dwi_frf_shells.bval\n"
            before_frf += "cp $bvec dwi_frf_shells.bvec\n"
        }
        if ( !wm_mask.empty() )
            args += " --mask_wm $wm_mask"

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py round $mask mask4scil.nii.gz --data_type uint8 -f
        $before_frf
        scil_compute_ssst_frf.py dwi_frf_shells.nii.gz dwi_frf_shells.bval dwi_frf_shells.bvec ${sid}_response.txt --mask mask4scil.nii.gz --fa $params.frf_fa --min_fa $params.frf_min_fa --min_nvox $params.frf_min_nvox $args
        """
}

process scilpy_msmt_response {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
    tuple val(sid), path(dwi), path(bval), path(bvec), path(mask), path(tissue_masks)
    val(caller_name)
    output:
    tuple val(sid), path("${sid}_wm_response.txt"), path("${sid}_gm_response.txt"), path("${sid}_csf_response.txt"), emit: response
    script:
        def args = ""
        if (params.frf_radii)
            args += " --roi_radii $params.frf_radii"
        if (params.frf_center)
            args += " --roi_center ${params.frf_center.join(" ")}"
        if (params.max_dti_bvalue)
            args += " --dti_bval_limit $params.max_dti_bvalue"
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py round $mask mask4scil.nii.gz --data_type uint8 -f
        scil_compute_msmt_frf.py $dwi $bval $bvec \
            ${sid}_wm_response.txt \
            ${sid}_gm_response.txt \
            ${sid}_csf_response.txt \
            --mask mask4scil.nii.gz \
            --mask_wm ${tissue_masks[0]} \
            --mask_gm ${tissue_masks[1]} \
            --mask_csf ${tissue_masks[2]} \
            --fa_thr_wm $params.frf_fa \
            --min_nvox $params.frf_min_nvox $args
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
        scil_image_math.py round $mask mask4scil.nii.gz --data_type uint8 -f
        scil_compute_ssst_fodf.py $dwi $bval $bvec $response ${sid}_fodf.nii.gz --mask mask4scil.nii.gz --force_b0_threshold --sh_order $params.sh_order --processes $task.cpus
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
    tuple val(sid), path("${sid}_vf.nii.gz"), path("${sid}_vf_rgb.nii.gz"), emit: vf
    script:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py round $mask mask4scil.nii.gz --data_type uint8 -f
        scil_compute_msmt_fodf.py $dwi $bval $bvec $wm_response $gm_response $csf_response \
            --wm_out_fODF ${sid}_wm_fodf.nii.gz \
            --gm_out_fODF ${sid}_gm_fodf.nii.gz \
            --csf_out_fODF ${sid}_csf_fodf.nii.gz \
            --vf ${sid}_vf.nii.gz \
            --vf_rgb ${sid}_vf_rgb.nii.gz \
            --mask mask4scil.nii.gz \
            --force_b0_threshold \
            --sh_order $params.sh_order \
            --processes $task.cpus
        """
}