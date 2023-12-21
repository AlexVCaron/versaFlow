#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.segmentation_classes = ["csf", "wm", "gm", "dgm", "pdgm", "blood"]
params.atropos_n4_tissues = ["wm", "dgm", "pdgm"]
params.atropos_prior_weight = 0.2
params.atropos_mrf_weight = 0.3
params.atropos_mrf_neighborhood = 1
params.atropos_n4_iterations = [200, 100, 50, 50]
params.atropos_n4_convergence_eps = 1E-10
params.atropos_n4_shrink_factor = 2
params.atropos_n4_bspline_fitting = 100
params.csf_distance_lambda = 0.1
params.csf_distance_probability = 0.5
params.default_distance_lambda = 0.1
params.default_distance_probability = 0.8
params.default_blur_factor = 0.4
params.csf_blur_factor = 0.25
params.random_seed = 1234

include { remove_alg_suffixes } from '../functions.nf'


/*process prepare_atropos {
    label "LIGHTSPEED"
    label "res_single_cpu"

    input:
        tuple val(sid), path(anat_images), path(segmentation)
    output:
        tuple val(sid), path("${sid}_run_atropos.sh"), emit: atropos_script
        tuple val(sid), path("${sid}_*_priors.nii.gz"), emit: segmentation_priors
    script:
        def args = ""
        if (!tissue_mappings.empty()) args += " --mappings $tissue_mappings"
        """
        mrhardi seg2pvf $args \
            --in $segmentation \
            --wm-label $params.segmentation_classes.indexOf("wm") \
            --gm-label $params.segmentation_classes.indexOf("gm") \
            
            --prefix ${segmentation.simpleName}_
        """
}*/


process atropos {
    label "ATROPOS"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/segmentation", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(t1_image), path(mask), path(segmentation)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_segmentation.nii.gz"), emit: segmentation
        tuple val(sid), path("${sid}_{${params.segmentation_classes.join(',')}}_pvf.nii.gz"), emit: vol_fractions
    script:
        def fractions_rename = []
        def tissue_label = params.segmentation_classes.indexOf("csf") + 1
        def prior_filename = "${segmentation.simpleName}_0${tissue_label}.nii.gz"
        def blurring_lines = ["scil_image_math.py blur $prior_filename $params.csf_blur_factor $prior_filename --data_type float32 -f"]
        def dist_priors = ["-l ${tissue_label}[$params.csf_distance_lambda,$params.csf_distance_probability]"]
        def i = 1
        for (cl in params.segmentation_classes) {
            if (cl != "csf") {
                tissue_label = params.segmentation_classes.indexOf(cl) + 1
                prior_filename = "${segmentation.simpleName}_0${tissue_label}.nii.gz"
                dist_priors += ["-l ${tissue_label}[$params.default_distance_lambda,$params.default_distance_probability]"]
                blurring_lines += ["scil_image_math.py blur $prior_filename $params.default_blur_factor $prior_filename --data_type float32 -f"]
            }

            fractions_rename += ["mv ${sid}_SegmentationPosteriors0${i}.nii.gz ${sid}_${cl}_pvf.nii.gz"]
            i += 1
        }
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed

        mrhardi seg2mask \
            --in $segmentation \
            --values ${(1..params.segmentation_classes.size()).join(',')} \
            --labels ${(1..params.segmentation_classes.size()).collect{ "0$it" }.join(',')} \
            --out ${segmentation.simpleName}

        ${blurring_lines.join('\n')}

        spacing=\$(mrinfo -spacing $t1_image | awk '{print \$1}')
        antsAtroposN4.sh -u 0 -d 3 \
            -a $t1_image \
            -x $mask \
            -b Aristotle[1] \
            ${dist_priors.join(' ')} \
            -r [$params.atropos_mrf_weight,$params.atropos_mrf_neighborhood] \
            ${params.atropos_n4_tissues.collect{ "-y ${params.segmentation_classes.indexOf(it)}" }.join(' ')} \
            -e [${params.atropos_n4_iterations.join('x')},$params.atropos_n4_convergence_eps] \
            -f $params.atropos_n4_shrink_factor \
            -c ${params.segmentation_classes.size()} \
            -p ${segmentation.simpleName}_%02d.nii.gz \
            -q [\$(echo "\${spacing[0]}*$params.atropos_n4_bspline_fitting" | bc)] \
            -o ${sid}_ \
            -w $params.atropos_prior_weight

        mv ${sid}_Segmentation.nii.gz tmp.nii.gz
        mv tmp.nii.gz ${sid}_segmentation.nii.gz

        ${fractions_rename.join('\n')}
        """
}
