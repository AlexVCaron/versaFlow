#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.eddy_with_reverse = true
params.eddy_select_gpu = true
params.use_cuda = false
params.eddy_force_shelled = true
params.b0_threshold = false
params.b0_normalization_strategy = "linear"
params.random_seed = 1234
params.verbose_outputs = false

include { remove_alg_suffixes } from '../functions.nf'

process dwi_denoise {
    label "MPCA_DENOISE"
    label params.on_hcp ? "res_full_node_override" : params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), file(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${dwi.simpleName}__dwidenoised.nii.gz"), emit: image
        tuple val(sid), path("${dwi.simpleName}__dwidenoised_metadata.*"), optional: true, emit: metadata
    script:
        def after_denoise = "fslmaths dwidenoise.nii.gz -thr 0 ${dwi.simpleName}__dwidenoised.nii.gz \n"
        if ( !metadata.empty() )
            after_denoise += "cp $metadata ${dwi.simpleName}__dwidenoised_metadata.py"

        def args = ""
        if ( !mask.empty() )
            args += " -mask $mask"

        """
        export MRTRIX_RNG_SEED=$params.random_seed
        dwidenoise -nthreads $task.cpus $args $dwi dwidenoise.nii.gz
        $after_denoise
        """
}

process nlmeans_denoise {
    label "NLMEANS_3D"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
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
        mrhardi nlmeans \
            --in $image \
            --out ${image.simpleName}__nlmeans_denoised \
            --processes $task.cpus $args
        $after_script
        """
}

process ants_gaussian_denoise {
    label "LONG"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}__ants_denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__ants_denoised_metadata.*"), optional: true, emit: metadata
    script:
        def args = ""
        if ( !mask.empty() )
            args += "--mask-image $mask"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed
        DenoiseImage --input-image $image \
            --noise-model Gaussian \
            --output [${image.simpleName}__ants_denoised.nii.gz,${image.simpleName}__ants_denoised_noise_map.nii.gz] \
            --verbose 1 $args
        """
}

process n4_denoise {
    label "N4_CORRECTION"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: {
        f -> ("$publish" == "true") ? f.contains("metadata") || f.contains("bias_field") ? null 
                                                                                         : remove_alg_suffixes(f)
                                    : null
    }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(anat), file(mask), file(metadata)
        val(caller_name)
        val(publish)
        path(config)
    output:
        tuple val(sid), path("${image.simpleName}__n4denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__n4denoised_metadata.*"), optional: true, emit: metadata
        tuple val(sid), path("${image.simpleName}_n4_bias_field.nii.gz"), emit: bias_field
    script:
        def after_denoise = ""
        def args = ""
        if ( anat.empty() ) {
            args += "--in $image"
            after_denoise += "mv n4denoise_bias_field.nii.gz ${image.simpleName}_n4_bias_field.nii.gz\n"
        }
        else {
            args += "--in $anat --apply $image"
            after_denoise += "mv tmp_n4denoised_bias_field.nii.gz ${image.simpleName}_n4_bias_field.nii.gz\n"
        }

        if ( !metadata.empty() ) {
            after_denoise += "mv n4denoise_metadata.py ${image.simpleName}__n4denoised_metadata.py\n"
            args += " --metadata $metadata"
        }
        after_denoise += "fslmaths n4denoise.nii.gz -thr 0 ${image.simpleName}__n4denoised.nii.gz\n"

        if ( !mask.empty() )
            args += " --mask $mask"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed
        mrhardi n4 $args \
            --out n4denoise \
            --config $config
        $after_denoise
        """
}

process apply_n4_bias_field {
    label "LIGHTSPEED"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), path(bias_field), file(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${image.simpleName}__n4denoised.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__n4denoised_metadata.*"), optional: true, emit: metadata
    script:
        def after_denoise = ""
        def args = ""
        if ( !metadata.empty() ) {
            after_denoise += "cp $metadata ${image.simpleName}__n4denoised_metadata.py\n"
        }

        if ( !mask.empty() )
            args += " --mask $mask"

        """
        scil_apply_bias_field_on_dwi.py $image $bias_field n4denoised.nii.gz -f $args
        fslmaths n4denoised.nii.gz -thr 0 ${image.simpleName}__n4denoised.nii.gz
        $after_denoise
        """    

}

process normalize_inter_b0 {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
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
        def args = ""
        def after_script = ""
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
        mrhardi b0 normalize \
            --in $dwi \
            --bvals $bval \
            --out ${dwi.simpleName}__inter_b0_normalized \
            --rout ${rev_dwi.simpleName}__inter_b0_normalized \
            --ref $params.b0_normalization_strategy $args 
        $after_script
        """
}

process prepare_epi_correction {
    label "LIGHTSPEED"
    label "res_single_cpu"

    input:
        tuple val(sid), path(b0s), path(dwi_bval), file(rev_bval), file(metadata)
        val(algo)
        path(config)
    output:
        tuple val(sid), path("${b0s.simpleName}__${algo}_script.sh"), emit: script
        tuple val(sid), path("${b0s.simpleName}__${algo}_acqp.txt"), emit: acqp
        tuple val(sid), path("${b0s.simpleName}__${algo}_config.cnf"), optional: true, emit: config
        tuple val(sid), val("${sid}__${algo}_results"), emit: awaited_out_name
        tuple val(sid), path("${b0s.simpleName}__${algo}_metadata.*"), emit: metadata
        tuple val(sid), path("${sid}__${algo}_metadata.*"), emit: metadata_for_corrected_dwi
        tuple val(sid), path("{${dwi_bval.collect{ it.simpleName }.join(",")},${rev_bval.collect{ it.simpleName }.join(",")}}_topup_indexes_metadata.*"), optional: true, emit : in_metadata_w_epi_correction
    script:
        """
        mrhardi epi $algo \
            --b0s $b0s \
            --bvals ${dwi_bval.join(',')} \
            --rev_bvals ${rev_bval.join(',')} \
            --out ${b0s.simpleName}__${algo} \
            --b0-thr ${params.b0_threshold ? params.b0_threshold : "0"} \
            --config $config \
            --verbose
        cp ${b0s.simpleName}__${algo}_metadata.py ${sid}__${algo}_metadata.py
        """
}

process topup {
    label "TOPUP"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
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

process bm_epi_correction {
    label "BM_EPI_CORRECTION"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    input:
        tuple val(sid), path(bm_script), path(b0), path(rev_b0), path(output_metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_b0__bm_corrected.nii.gz"), emit: image
        tuple val(sid), path("${sid}_b0__bm_field.nii.gz"), emit: displacement_field
        tuple val(sid), path("${sid}_b0__bm_fieldmap.nii.gz"), emit: fieldmap
        tuple val(sid), path(output_metadata), optional: true, emit: metadata
    script:

        """
        ./$bm_script $b0 $rev_b0 ${sid}_b0_ $task.cpus
        """
}

process prepare_eddy {
    label "LIGHTSPEED"
    label "res_single_cpu"

    input:
        tuple val(sid), val(prefix), file(acqp), val(rev_prefix), path(data), path(metadata)
        path(config)
    output:
        tuple val(sid), path("${prefix}__eddy_script.sh"), path("${prefix}__eddy_index.txt"), path("${prefix}__eddy_acqp.txt"), emit: config
        tuple val(sid), path("${prefix}__eddy_slspec.txt"), emit: slspec, optional: true
        tuple val(sid), path("${sid}*non_zero.bvec"), emit: bvec, optional: true
        tuple val(sid), path("${prefix}__eddy_metadata.*"), emit: metadata, optional: true
    script:
        def args = ""
        def after_script = ""
        def will_gen_acqp = true
        if ( !acqp.empty() ) {
            args += " --acqp $acqp"
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

        if ( params.verbose_outputs ) {
            args += " --debug"
        }

        if ( params.eddy_force_shelled )
            args += " --shelled"

        if ( !will_gen_acqp )
            after_script += "cp $acqp ${prefix}__eddy_acqp.txt\n"

        """
        mrhardi eddy \
            --in $prefix \
            --out ${prefix}__eddy \
            --config $config \
            --seed $args
        $after_script
        """
}

process eddy {
    label params.use_cuda ? "EDDY_GPU" : "EDDY_OMP"
    label params.use_cuda ? "res_single_cpu" : params.on_hcp ? "res_full_node_override" : params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"
    label params.use_cuda ? "res_gpu" : ""

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(eddy_script), path(eddy_index), path(eddy_acqp), file(eddy_slspec), file(epi_field), file(disp_field), path(dwi), path(bval), path(bvec), path(mask), val(topup_prefix), file(topup_package), path(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.nii.gz"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.bval"), emit: bval
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected.bvec"), emit: bvec
        tuple val(sid), path("${dwi.simpleName}__eddy_corrected_metadata.py"), optional: true, emit: metadata
    script:
        def after_script = ""
        def after_eddy = ""
        if ( metadata )
            after_script += "cp $metadata ${dwi.simpleName}__eddy_corrected_metadata.py"

        def args = "eddy_in_image.nii.gz $bval $bvec"
        def kwargs = ""

        if ( mask ) {
            args += " $mask"
        }

        args += " $eddy_acqp $eddy_index"

        if ( !epi_field.empty() ) {
            if ( !disp_field.empty() ) after_eddy += "animaApplyDistortionCorrection -f eddy_corrected.nii.gz -t $disp_field -o eddy_corrected.nii.gz -T $task.cpus\n"
        }
        else if ( topup_prefix ) {
            kwargs += " --topup $topup_prefix"
        }

        if ( !eddy_slspec.empty() )
            kwargs += " --slspec $eddy_slspec"

        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        fslmaths $dwi -thr 0 eddy_in_image.nii.gz
        ./$eddy_script $args eddy_corrected $kwargs
        $after_eddy
        mv eddy_corrected.eddy_rotated_bvecs ${dwi.simpleName}__eddy_corrected.bvec
        cp $bval ${dwi.simpleName}__eddy_corrected.bval
        cp eddy_corrected.nii.gz ${dwi.simpleName}__eddy_corrected.nii.gz
        fslmaths ${dwi.simpleName}__eddy_corrected.nii.gz -thr 0 ${dwi.simpleName}__eddy_corrected.nii.gz
        $after_script
        """
}

process gibbs_removal {
    label "MEDIUM"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${dwi.simpleName}__gibbs_corrected.nii.gz"), emit: image
        tuple val(sid), path("${dwi.simpleName}__gibbs_corrected_metadata.*"), optional: true, emit: metadata
    script:
    def after_denoise = "fslmaths gibbs_corrected.nii.gz -thr 0 ${dwi.simpleName}__gibbs_corrected.nii.gz\n"
    if ( metadata )
        after_denoise += "cp $metadata ${dwi.simpleName}__gibbs_corrected_metadata.py"

    """
    export MRTRIX_RNG_SEED=$params.random_seed
    mrdegibbs -nthreads $task.cpus $dwi gibbs_corrected.nii.gz
    $after_denoise
    """
}