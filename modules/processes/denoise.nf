#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.eddy_with_reverse = true
params.eddy_select_gpu = true
params.use_cuda = false
params.eddy_force_shelled = true
params.b0_threshold = false
params.b0_normalization_strategy = "linear"

include { get_size_in_gb; remove_alg_suffixes } from '../functions.nf'

process dwi_denoise {
    label params.on_hcp ? "res_full_node_override" : params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), file(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${dwi.simpleName}__dwidenoised.nii.gz"), emit: image
        tuple val(sid), path("${dwi.simpleName}__dwidenoised_metadata.*"), optional: true, emit: metadata
    script:
        after_denoise = "fslmaths -dt double dwidenoise.nii.gz -thr 0 ${dwi.simpleName}__dwidenoised.nii.gz -odt double\n"
        if ( !metadata.empty() )
            after_denoise += "cp $metadata ${dwi.simpleName}__dwidenoised_metadata.py"

        args = "-nthreads $task.cpus -datatype float64"
        if ( !mask.empty() )
            args += " -mask $mask"

        """
        dwidenoise $args $dwi dwidenoise.nii.gz
        $after_denoise
        """
}

process nlmeans_denoise {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${image.simpleName}__nlmeans_denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__nlmeans_denoised_metadata.*"), optional: true, emit: metadata
    script:
        def args = ""
        if ( !mask.empty() ) args += "--mask $mask"
        def after_script = ""
        if ( !metadata.empty() ) after_script += "cp $metadata ${image.simpleName}__nlmeans_denoised_metadata.py"
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        scil_run_nlmeans.py $image ${image.simpleName}__nlmeans_denoised.nii.gz 1 --processes $task.cpus -f $args
        $after_script
        """
}

process ants_gaussian_denoise {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}__ants_denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__ants_denoised_metadata.*"), optional: true, emit: metadata
    script:
        args = ""
        if ( !mask.empty() )
            args += "--mask-image $mask"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        DenoiseImage --input-image $image --noise-model Gaussian --output [${image.simpleName}__ants_denoised.nii.gz,${image.simpleName}__ants_denoised_noise_map.nii.gz] --verbose 1 $args
        """
}

process n4_denoise {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(anat), file(mask), file(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${image.simpleName}__n4denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__n4denoised_metadata.*"), optional: true, emit: metadata
    script:
        after_denoise = ""
        args = ""
        if ( anat.empty() )
            args += "--in $image"
        else
            args += "--in $anat --apply $image"

        if ( !metadata.empty() ) {
            after_denoise += "mv n4denoise_metadata.py ${image.simpleName}__n4denoised_metadata.py\n"
            args += " --metadata $metadata"
        }
        after_denoise += "fslmaths -dt double n4denoise.nii.gz -thr 0 ${image.simpleName}__n4denoised.nii.gz -odt double\n"

        if ( !mask.empty() )
            args += " --mask $mask"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        magic-monkey n4 $args --out n4denoise --config $config
        $after_denoise
        """
}

process normalize_inter_b0 {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("${rev_dwi.simpleName}") ? null : f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), file(rev_dwi), file(rev_bval), file(dwi_metadata), file(rev_metadata)
        val(caller_name)
        file(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__inter_b0_normalized.nii.gz"), emit: dwi
        tuple val(sid), path("${rev_dwi.simpleName}__inter_b0_normalized.nii.gz"), optional: true, emit: rev
        tuple val(sid), path("${dwi.simpleName}*_metadata.*"), optional: true, emit: dwi_metadata
        tuple val(sid), path("${rev_dwi.simpleName}*_metadata.*"), optional: true, emit: rev_metadata
    script:
        args = "--in $dwi --bvals $bval"
        after_script = ""
        if ( !rev_dwi.empty() )
            args += " --rev $rev_dwi"
        if ( !rev_bval.empty() )
            args += " --rvals $rev_bval"
        if ( !config.empty() )
            args += " --config $config"

        if (!dwi_metadata.empty())
            after_script += "cp $dwi_metadata ${dwi.simpleName}__inter_b0_normalized_metadata.py\n"
        if ( !rev_metadata.empty())
            after_script += "cp $rev_metadata ${rev_dwi.simpleName}__inter_b0_normalized_metadata.py\n"

        if (params.b0_threshold)
            args += " --ceil ${params.b0_threshold}"

        """
        magic-monkey b0 normalize $args --out ${dwi.simpleName}__inter_b0_normalized --rout ${rev_dwi.simpleName}__inter_b0_normalized --ref $params.b0_normalization_strategy
        $after_script
        """
}

process prepare_topup {
    label "res_single_cpu"

    input:
        tuple val(sid), path(b0s), path(dwi_bval), file(rev_bval), file(metadata)
        path(config)
    output:
        tuple val(sid), path("${b0s.simpleName}__topup_script.sh"), path("${b0s.simpleName}__topup_acqp.txt"), path("${b0s.simpleName}__topup_config.cnf"), val("${sid}__topup_results"), emit: config
        tuple val(sid), path("${b0s.simpleName}__topup_metadata.*"), emit: metadata
        tuple val(sid), path("{${dwi_bval.collect{ it.simpleName }.join(",")},${rev_bval.collect{ it.simpleName }.join(",")}}_topup_indexes_metadata.*"), optional: true, emit : in_metadata_w_topup
    script:
        """
        magic-monkey topup --b0s $b0s --bvals ${dwi_bval.join(',')} --rev_bvals ${rev_bval.join(',')} --out ${b0s.simpleName}__topup --config $config --verbose
        """
}

process topup {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("b0") ? null : f.contains("metadata") ? null : f.contains("topup.nii.gz") ? remove_alg_suffixes(f): null }, mode: params.publish_mode

    input:
        tuple val(sid), path(topup_script), path(topup_acqp), path(topup_cnf), path(b0), path(output_metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_b0__topup.nii.gz"), emit: image
        tuple val(sid), path("${sid}__topup_field.nii.gz"), emit: field
        tuple val(sid), path("${sid}__topup_results_movpar.txt"), path("${sid}__topup_results_fieldcoef.nii.gz"), emit: transfo
        tuple val(sid), path("${sid}_b0__topup.nii.gz"), path("${sid}__topup_field.nii.gz"), path("${sid}__topup_results_movpar.txt"), path("${sid}__topup_results_fieldcoef.nii.gz"), emit: pkg
        tuple val(sid), path(output_metadata), optional: true, emit: metadata
    script:
        """
        fslmaths $b0 -thr 0 topup_in_image.nii.gz
        ./$topup_script topup_in_image.nii.gz ${sid}__topup
        mv ${sid}__topup.nii.gz ${sid}_b0__topup.nii.gz
        """
}

process prepare_eddy {
    label "res_single_cpu"

    input:
        tuple val(sid), val(prefix), file(topup_acqp), val(rev_prefix), path(data), path(metadata)
        path(config)
    output:
        tuple val(sid), path("${prefix}__eddy_script.sh"), path("${prefix}__eddy_index.txt"), path("${prefix}__eddy_acqp.txt"), emit: config
        tuple val(sid), path("${prefix}__eddy_slspec.txt"), emit: slspec, optional: true
        tuple val(sid), path("${sid}*non_zero.bvec"), emit: bvec, optional: true
        tuple val(sid), path("${prefix}__eddy_metadata.*"), emit: metadata, optional: true
    script:
        args = "--in $prefix --debug"
        will_gen_acqp = true
        if ( !topup_acqp.empty() ) {
            args += " --acqp $topup_acqp"
            will_gen_acqp = false
        }
        if ( rev_prefix ) {
            args += " --rev $rev_prefix"
            if ( params.eddy_with_reverse )
                args += " --rev_eddy"
        }

        if ( params.use_cuda ) {
            args += " --cuda"
            if ( !params.eddy_select_gpu ) {
                args += " --dont_gpu"
            }
        }

        if ( params.eddy_force_shelled )
            args += " --shelled"

        if ( will_gen_acqp )
            """
            magic-monkey eddy $args --out ${prefix}__eddy --config $config
            """
        else
            """
            magic-monkey eddy $args --out ${prefix}__eddy --config $config && cp $topup_acqp "${prefix}__eddy_acqp.txt"
            """
}

process eddy {
    label params.use_cuda ? "res_single_cpu" : params.on_hcp ? "res_full_node_override" : "res_max_cpu"
    label params.use_cuda ? "res_gpu" : ""

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(eddy_script), path(eddy_index), path(eddy_acqp), file(eddy_slspec), path(dwi), path(bval), path(bvec), path(mask), val(topup_prefix), file(topup_package), path(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.nii.gz"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.bval"), emit: bval
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.bvec"), emit: bvec
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected_metadata.py"), optional: true, emit: metadata
    script:
        def after_script = ""
        if ( metadata )
            after_script += "cp $metadata ${dwi.simpleName}__eddy_corrected_metadata.py"

        def args = "eddy_in_image.nii.gz $bval $bvec"

        if ( mask ) {
            args += " $mask"
        }

        args += " $eddy_acqp $eddy_index"

        if ( topup_prefix ) {
            args += " --topup $topup_prefix"
        }

        if ( !eddy_slspec.empty() )
            args += " --slspec $eddy_slspec"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        fslmaths $dwi -thr 0 eddy_in_image.nii.gz
        ./$eddy_script $args eddy_corrected
        mv eddy_corrected.eddy_rotated_bvecs ${dwi.simpleName}__eddy_corrected.bvec
        cp $bval ${dwi.simpleName}__eddy_corrected.bval
        cp eddy_corrected.nii.gz ${dwi.simpleName}__eddy_corrected.nii.gz
        fslmaths ${dwi.simpleName}__eddy_corrected.nii.gz -thr 0 ${dwi.simpleName}__eddy_corrected.nii.gz
        $after_script
        """
}

process gibbs_removal {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${dwi.simpleName}__gibbs_corrected.nii.gz"), emit: image
        tuple val(sid), path("${dwi.simpleName}__gibbs_corrected_metadata.*"), optional: true, emit: metadata
    script:
    after_denoise = "fslmaths -dt double gibbs_corrected.nii.gz -thr 0 ${dwi.simpleName}__gibbs_corrected.nii.gz -odt double\n"
    if ( metadata )
        after_denoise += "cp $metadata ${dwi.simpleName}__gibbs_corrected_metadata.py"

    """
    mrdegibbs -nthreads 1 -datatype float64 $dwi gibbs_corrected.nii.gz
    $after_denoise
    """
}