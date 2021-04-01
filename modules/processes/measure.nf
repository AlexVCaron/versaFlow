#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.verbose_outputs = true
params.fodf_max_absolute_factor = 2.0
params.fodf_relative_thr = 0.1
params.max_fa_ventricle = 0.1
params.min_md_ventricle = 0.003


include { get_size_in_gb; swap_configurations } from '../functions.nf'

process dti_metrics {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), val(input_prefix), file(mask), path(data), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), val("${sid}_dti_metrics"), emit: prefix
        tuple val(sid), path("${sid}_dti_metrics*.nii.gz"), emit: metrics
    script:
        """
        magic-monkey dti_metrics --in $input_prefix --out ${sid}_dti_metrics --config $config
        """
}

process scil_compute_dti_fa {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(mask)
        val(processing_caller_name)
        val(measuring_caller_name)
    output:
        tuple val(sid), val("${sid}_dti"), emit: prefix
        tuple val(sid), path("${sid}_dti_dti.nii.gz"), emit: dti
        tuple val(sid), path("${sid}_dti_fa.nii.gz"), emit: fa
        tuple val(sid), path("${sid}_dti_md.nii.gz"), emit: md
    script:
        def avail_threads = Math.round(task.cpus / 3)
        def remainder_threads = task.cpus - avail_threads
        def args = "--tensor ${sid}_dti_dti.nii.gz"
        args += " --fa ${sid}_dti_fa.nii.gz --md ${sid}_dti_md.nii.gz"
        before = ""
        if ( !mask.empty() ) {
            args += " --mask $mask"
            before += "mrconvert -datatype uint8 $mask mask4scil.nii.gz\n"
        }

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${avail_threads + remainder_threads}
        export OMP_NUM_THREADS=$avail_threads
        export OPENBLAS_NUM_THREADS=1
        $before
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        scil_compute_dti_metrics.py $dwi $bval flipped_bvecs.bvec -f --not_all $args
        """
}

process scil_dti_and_metrics {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$processing_caller_name/${task.index}_${task.process.replaceAll(":", "_")}", saveAs: { f -> f.contains("dti_dti") ? f : f.contains("metadata") ? f : null }, mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/all/${sid}/$measuring_caller_name/${task.index}_${task.process.replaceAll(":", "_")}",saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f },  mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("dti_dti") ? f : null }, mode: params.publish_mode
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(mask)
        val(processing_caller_name)
        val(measuring_caller_name)
    output:
        tuple val(sid), val("${sid}_dti"), emit: prefix
        tuple val(sid), path("${sid}_dti_dti.nii.gz"), emit: dti
        tuple val(sid), path("${sid}_dti_evals.nii.gz"), path("${sid}_dti_evecs.nii.gz"), path("${sid}_dti_evals_*.nii.gz"), path("${sid}_dti_evecs_*.nii.gz"), emit: eigen
        tuple val(sid), path("${sid}_dti_fa.nii.gz"), path("${sid}_dti_ga.nii.gz"), path("${sid}_dti_rgb.nii.gz"), emit: aniso
        tuple val(sid), path("${sid}_dti_md.nii.gz"), path("${sid}_dti_ad.nii.gz"), path("${sid}_dti_rd.nii.gz"), path("${sid}_dti_mode.nii.gz"), path("${sid}_dti_norm.nii.gz"), emit: iso
        tuple val(sid), path("${sid}_dti_non_physical.nii.gz"), path("${sid}_dti_pulsation*.nii.gz"), emit: artifacts, optional: true
        tuple val(sid), path("${sid}_dti_residuals.nii.gz"), path("${sid}_dti_residuals*.nii.gz"), emit: residuals, optional: true
    script:
        def avail_threads = Math.round(task.cpus / 3)
        def remainder_threads = task.cpus - avail_threads
        def args = "--tensor ${sid}_dti_dti.nii.gz --evals ${sid}_dti_evals.nii.gz --evecs ${sid}_dti_evecs.nii.gz"
        args += " --fa ${sid}_dti_fa.nii.gz --ga ${sid}_dti_ga.nii.gz --rgb ${sid}_dti_rgb.nii.gz"
        args += " --md ${sid}_dti_md.nii.gz --ad ${sid}_dti_ad.nii.gz --rd ${sid}_dti_rd.nii.gz --mode ${sid}_dti_mode.nii.gz --norm ${sid}_dti_norm.nii.gz"
        if ( params.verbose_outputs )
            args += " --residual ${sid}_dti_residuals.nii.gz --non-physical ${sid}_dti_non_physical.nii.gz --pulsation ${sid}_dti_pulsation.nii.gz"

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${avail_threads + remainder_threads}
        export OMP_NUM_THREADS=$avail_threads
        export OPENBLAS_NUM_THREADS=1
        mrconvert -datatype uint8 $mask mask4scil.nii.gz
        magic-monkey flip2ref --in $dwi --bvecs $bvec --out flipped_bvecs
        scil_compute_dti_metrics.py $dwi $bval flipped_bvecs.bvec --mask mask4scil.nii.gz -f --not_all $args
        """
}

process diamond_metrics {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/diamond", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), val(input_prefix), file(mask), path(data), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), val("${sid}_diamond_metrics"), emit: prefix
        tuple val(sid), path("${sid}_diamond_metrics*.nii.gz"), emit: metrics
    script:
        def args = "--in $input_prefix"
        if ( !mask.empty() ) {
            args += " --mask $mask"
        }
        """
        magic-monkey diamond_metrics $args --out ${sid}_diamond_metrics --config $config
        """
}

process odf_metrics {
    label params.conservative_resources ? "res_conservative" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(odfs), path(fa), path(md), path(mask)
        val(caller_name)
        val(basis)
    output:
        tuple val(sid), val("${sid}_fodf_metrics"), emit: prefix
        tuple val(sid), path("${sid}_fodf_metrics*.nii.gz"), emit: metrics
    script:
        """
        scil_compute_fodf_max_in_ventricles.py $odfs $fa $md --max_value_output vmax.txt --sh_basis descoteaux07 --fa_t $params.max_fa_ventricle --md_t $params.min_md_ventricle -f
        abs_threshold=\$(echo $params.fodf_max_absolute_factor*\$(cat vmax.txt)|bc)
        scil_compute_fodf_metrics.py --rt $params.fodf_relative_thr --at \${abs_threshold} --sh_basis $basis --mask $mask --afd_max ${sid}_fodf_metrics_afd.nii.gz --afd_total ${sid}_fodf_metrics_afdt.nii.gz --afd_sum ${sid}_fodf_metrics_afds.nii.gz --nufo ${sid}_fodf_metrics_nufo.nii.gz --peaks ${sid}_fodf_metrics_peaks.nii.gz --rgb ${sid}_fodf_metrics_rgb.nii.gz --peak_values ${sid}_fodf_metrics_peaks_values.nii.gz --peak_indices ${sid}_fodf_metrics_peaks_indices.nii.gz $odfs
        """
}

