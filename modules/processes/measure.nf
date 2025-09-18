#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.verbose_outputs = true
params.fodf_wm_max_absolute_factor = 2.0
params.fodf_wm_relative_thr = 0.1
params.fodf_gm_max_absolute_factor = 1.0
params.fodf_gm_relative_thr = 0.6
params.ventricles_center = false
params.max_fa_ventricle = 0.1
params.min_md_ventricle = 0.003
params.max_dti_bvalue = 1300
params.random_seed = 1234
params.b0_threshold = false


process dti_metrics {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), val(input_prefix), file(mask), path(data), path(metadata)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), val("${sid}_dti_metrics"), emit: prefix
        tuple val(sid), path("${sid}_dti_metrics*.nii.gz"), emit: metrics
    script:
        """
        export MRTRIX_RNG_SEED=$params.random_seed
        mrhardi dti_metrics --in $input_prefix --out ${sid}_dti_metrics --config $config
        """
}

process scil_compute_dti_fa {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$processing_caller_name/${task.process.replaceAll(":", "/")}", saveAs: { f -> f.contains("dti_dti") ? f : f.contains("metadata") ? f : null }, mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/all/${sid}/$measuring_caller_name/${task.process.replaceAll(":", "/")}",saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f },  mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> ("$publish" == "true") ? f.contains("dti_dti") ? f : null : null }, mode: params.publish_mode, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> ("$publish" == "true") ? f.contains("dti_dti") ? null : f.contains("metadata") ? null : f : null }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(mask)
        val(processing_caller_name)
        val(measuring_caller_name)
        val(publish)
    output:
        tuple val(sid), val("${sid}_dti"), emit: prefix
        tuple val(sid), path("${sid}_dti_dti.nii.gz"), emit: dti
        tuple val(sid), path("${sid}_dti_fa.nii.gz"), emit: fa
        tuple val(sid), path("${sid}_dti_md.nii.gz"), emit: md
        tuple val(sid), path("${sid}_dti_evecs_v1.nii.gz"), emit: main_peak
        tuple val(sid), path("${sid}_dti_evecs*.nii.gz"), emit: evecs
    script:
        def args = "--tensor ${sid}_dti_dti.nii.gz"
        args += " --fa ${sid}_dti_fa.nii.gz --md ${sid}_dti_md.nii.gz"
        args += " --evecs ${sid}_dti_evecs.nii.gz"
        def before = ""
        if ( !mask.empty() ) {
            args += " --mask $mask"
            before += "scil_volume_math.py round $mask mask4scil.nii.gz --data_type uint8 -f\n"
        }

        if ( params.max_dti_bvalue ) {
            def shell_args = ""
            if (params.b0_threshold)
                shell_args += " --ceil ${params.b0_threshold}"
            before += "mrhardi shells --in $dwi --bvals $bval --bvecs $bvec --shells $params.max_dti_bvalue --keep leq --out dwi_for_dti --with_b0 $shell_args\n"
        }
        else {
            before += "cp $dwi dwi_for_dti.nii.gz\ncp $bval dwi_for_dti.bval\ncp $bvec dwi_for_dti.bvec"
        }

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        $before
        scil_dti_metrics.py dwi_for_dti.nii.gz dwi_for_dti.bval dwi_for_dti.bvec -f --not_all $args
        """
}

process scil_compute_dti_fa_np {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$processing_caller_name/${task.process.replaceAll(":", "/")}", saveAs: { f -> f.contains("dti_dti") ? f : f.contains("metadata") ? f : null }, mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/all/${sid}/$measuring_caller_name/${task.process.replaceAll(":", "/")}",saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f },  mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> ("$publish" == "true") ? f.contains("dti_dti") ? f : null : null }, mode: params.publish_mode, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> ("$publish" == "true") ? f.contains("dti_dti") ? null : f.contains("metadata") ? null : f : null }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(mask)
        val(processing_caller_name)
        val(measuring_caller_name)
        val(publish)
    output:
        tuple val(sid), val("${sid}_dti"), emit: prefix
        tuple val(sid), path("${sid}_dti_dti.nii.gz"), emit: dti
        tuple val(sid), path("${sid}_dti_fa.nii.gz"), emit: fa
        tuple val(sid), path("${sid}_dti_md.nii.gz"), emit: md
        tuple val(sid), path("${sid}_dti_evecs_v1.nii.gz"), emit: main_peak
        tuple val(sid), path("${sid}_dti_evecs*.nii.gz"), emit: evecs
        tuple val(sid), path("${sid}_dti_non_physical.nii.gz"), emit: np_outliers_mask
    script:
        def args = "--tensor ${sid}_dti_dti.nii.gz"
        args += " --fa ${sid}_dti_fa.nii.gz --md ${sid}_dti_md.nii.gz"
        args += " --evecs ${sid}_dti_evecs.nii.gz"
        args += " --non-physical ${sid}_dti_non_physical.nii.gz"
        def before = ""
        if ( !mask.empty() ) {
            args += " --mask $mask"
            before += "scil_volume_math.py round $mask mask4scil.nii.gz --data_type uint8 -f\n"
        }

        if ( params.max_dti_bvalue ) {
            def shell_args = ""
            if (params.b0_threshold)
                shell_args += " --ceil ${params.b0_threshold}"
            before += "mrhardi shells --in $dwi --bvals $bval --bvecs $bvec --shells $params.max_dti_bvalue --keep leq --out dwi_for_dti --with_b0 $shell_args\n"
        }
        else {
            before += "cp $dwi dwi_for_dti.nii.gz\ncp $bval dwi_for_dti.bval\ncp $bvec dwi_for_dti.bvec"
        }

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        $before
        scil_dti_metrics.py dwi_for_dti.nii.gz dwi_for_dti.bval dwi_for_dti.bvec -f --not_all $args
        """
}

process scil_dti_and_metrics {
    label "LONG"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$processing_caller_name/${task.process.replaceAll(":", "/")}", saveAs: { f -> f.contains("dti_dti") ? f : f.contains("metadata") ? f : null }, mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/all/${sid}/$measuring_caller_name/${task.process.replaceAll(":", "/")}",saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f },  mode: params.publish_mode, enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("dti_dti") ? f : null }, mode: params.publish_mode, overwrite: true
    publishDir "${params.output_root}/${sid}/dti", saveAs: { f -> f.contains("dti_dti") ? null : f.contains("metadata") ? null : f }, mode: params.publish_mode, overwrite: true

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
        def args = "--tensor ${sid}_dti_dti.nii.gz --evals ${sid}_dti_evals.nii.gz --evecs ${sid}_dti_evecs.nii.gz"
        args += " --fa ${sid}_dti_fa.nii.gz --ga ${sid}_dti_ga.nii.gz --rgb ${sid}_dti_rgb.nii.gz"
        args += " --md ${sid}_dti_md.nii.gz --ad ${sid}_dti_ad.nii.gz --rd ${sid}_dti_rd.nii.gz --mode ${sid}_dti_mode.nii.gz --norm ${sid}_dti_norm.nii.gz"
        args += " --residual ${sid}_dti_residuals.nii.gz"
        if ( params.verbose_outputs )
            args += " --non-physical ${sid}_dti_non_physical.nii.gz --pulsation ${sid}_dti_pulsation.nii.gz"

        before = ""
        if ( params.max_dti_bvalue ) {
            def shell_args = ""
            if (params.b0_threshold)
                shell_args += " --ceil ${params.b0_threshold}"
            before += "mrhardi shells --in $dwi --bvals $bval --bvecs $bvec --shells $params.max_dti_bvalue --keep leq --out dwi_for_dti --with_b0 $shell_args\n"
        }
        else {
            before += "cp $dwi dwi_for_dti.nii.gz\ncp $bval dwi_for_dti.bval\ncp $bvec dwi_for_dti.bvec"
        }

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_volume_math.py floor $mask mask4scil.nii.gz --data_type uint8 -f
        $before
        scil_dti_metrics.py dwi_for_dti.nii.gz dwi_for_dti.bval dwi_for_dti.bvec --mask mask4scil.nii.gz -f --not_all $args
        """
}

process diamond_metrics {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/diamond", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), val(input_prefix), file(mask), path(data), path(xml_summary), path(metadata)
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
        mrhardi diamond_metrics $args --out ${sid}_diamond_metrics --xml-config $xml_summary --config $config
        """
}

process odf_metrics {
    label "MEDIUM"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(wm_odfs), file(gm_odfs), file(csf_odfs), path(fa), path(md), path(wm_mask), path(gm_mask)
        val(caller_name)
        val(basis)
    output:
        tuple val(sid), val("${sid}_fodf_metrics"), emit: prefix
        tuple val(sid), path("${sid}_fodf_metrics*.nii.gz"), emit: metrics
        tuple val(sid), path("${sid}_ventricles_mask.nii.gz"), path("${sid}_ventricles_fodf_max.txt"), emit: ventricles
    script:
        def args = ""
        def csf_f = csf_odfs.empty() ? "$wm_odfs" : "$csf_odfs"
        def gm_f = gm_odfs.empty() ? "$wm_odfs" : "$gm_odfs"
        if ( params.ventricles_center )
            args += " --center ${ params.ventricles_center.join(' ') }"
        """
        scil_fodf_max_in_ventricles.py $csf_f $fa $md \
            --max_value_output ${sid}_ventricles_fodf_max.txt \
            --sh_basis descoteaux07 \
            --fa_t $params.max_fa_ventricle \
            --md_t $params.min_md_ventricle \
            --mask_output ${sid}_ventricles_mask.nii.gz \
            $args \
            -f

        wm_abs_threshold=\$(echo $params.fodf_wm_max_absolute_factor*\$(cat ${sid}_ventricles_fodf_max.txt)|bc)
        gm_abs_threshold=\$(echo $params.fodf_gm_max_absolute_factor*\$(cat ${sid}_ventricles_fodf_max.txt)|bc)

        scil_fodf_metrics.py $wm_odfs \
            --rt $params.fodf_wm_relative_thr \
            --at \${wm_abs_threshold} \
            --sh_basis $basis \
            --mask $wm_mask \
            --afd_max ${sid}_fodf_metrics_wm_afd.nii.gz \
            --afd_total ${sid}_fodf_metrics_wm_afdt.nii.gz \
            --afd_sum ${sid}_fodf_metrics_wm_afds.nii.gz \
            --nufo ${sid}_fodf_metrics_wm_nufo.nii.gz \
            --peaks ${sid}_fodf_metrics_wm_peaks.nii.gz \
            --rgb ${sid}_fodf_metrics_wm_rgb.nii.gz \
            --peak_values ${sid}_fodf_metrics_wm_peaks_values.nii.gz \
            --peak_indices ${sid}_fodf_metrics_wm_peaks_indices.nii.gz \
            --processes $task.cpus

        scil_fodf_metrics.py $gm_f \
            --rt $params.fodf_gm_relative_thr \
            --at \${gm_abs_threshold} \
            --sh_basis $basis \
            --mask $gm_mask \
            --afd_max ${sid}_fodf_metrics_gm_afd.nii.gz \
            --afd_total ${sid}_fodf_metrics_gm_afdt.nii.gz \
            --afd_sum ${sid}_fodf_metrics_gm_afds.nii.gz \
            --nufo ${sid}_fodf_metrics_gm_nufo.nii.gz \
            --peaks ${sid}_fodf_metrics_gm_peaks.nii.gz \
            --rgb ${sid}_fodf_metrics_gm_rgb.nii.gz \
            --peak_values ${sid}_fodf_metrics_gm_peaks_values.nii.gz \
            --peak_indices ${sid}_fodf_metrics_gm_peaks_indices.nii.gz \
            --processes $task.cpus
        """
}

