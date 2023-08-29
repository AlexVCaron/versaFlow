#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.segmentation_classes = ["csf", "gm", "wm", "dgm"]
params.atropos_prior_weight = 0.0
params.atropos_n4_bspline_spacing = 100.0
params.random_seed = 1234

include { remove_alg_suffixes } from '../functions.nf'

process atropos {
    label "ATROPOS"
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/segmentation", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(t1_image), path(mask), path(segmentation)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_segmentation.nii.gz"), emit: segmentation
        tuple val(sid), path("${sid}_{${params.segmentation_classes.join(',')}}_pvf.nii.gz"), emit: vol_fractions
    script:
        def after_script = ""
        def i = 1
        for (cl in params.segmentation_classes) {
            after_script += "mv ${sid}_SegmentationPosteriors0${i}.nii.gz ${sid}_${cl}_pvf.nii.gz\n"
            i += 1
        }
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=$params.random_seed
        mrhardi seg2mask --in $segmentation --values 1,2,3,4,5 --labels 01,02,04,03,05 --out ${segmentation.simpleName}
        scil_image_math.py addition ${segmentation.simpleName}_01.nii.gz ${segmentation.simpleName}_05.nii.gz ${segmentation.simpleName}_01.nii.gz --data_type float32 -f
        rm ${segmentation.simpleName}_05.nii.gz
        scil_image_math.py blur ${segmentation.simpleName}_01.nii.gz 1 ${segmentation.simpleName}_01.nii.gz --data_type float32 -f
        scil_image_math.py blur ${segmentation.simpleName}_02.nii.gz 1 ${segmentation.simpleName}_02.nii.gz --data_type float32 -f
        scil_image_math.py blur ${segmentation.simpleName}_03.nii.gz 1 ${segmentation.simpleName}_03.nii.gz --data_type float32 -f
        scil_image_math.py blur ${segmentation.simpleName}_04.nii.gz 1 ${segmentation.simpleName}_04.nii.gz --data_type float32 -f
        spacing=\$(mrinfo -spacing $t1_image | awk '{print \$1}')
        antsAtroposN4.sh -u 0 -d 3 \
            -a $t1_image \
            -x $mask \
            -c ${params.segmentation_classes.size()} \
            -p ${segmentation.simpleName}_%02d.nii.gz \
            -q \$(echo "\${spacing[0]}*$params.atropos_n4_bspline_spacing" | bc) \
            -o ${sid}_ \
            -w $params.atropos_prior_weight
        mv ${sid}_Segmentation.nii.gz tmp.nii.gz
        mv tmp.nii.gz ${sid}_segmentation.nii.gz
        $after_script
        """
}
