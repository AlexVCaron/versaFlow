#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.segmentation_classes = ["csf", "gm", "wm"]
params.segmentation_mask_threshold = 1e-3

include { remove_alg_suffixes } from '../functions.nf'

process atropos {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/segmentation", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(t1_image), path(mask), path(segmentation)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_segmentation.nii.gz"), emit: segmentation
        tuple val(sid), path("${sid}_{${params.segmentation_classes.join(',')}}_mask.nii.gz"), emit: masks
        tuple val(sid), path("${sid}_{${params.segmentation_classes.join(',')}}_vf.nii.gz"), emit: vol_fractions
    script:
        after_script = ""
        i = 1
        for (cl in params.segmentation_classes) {
            after_script += "mv ${sid}_SegmentationPosteriors0${i}.nii.gz ${sid}_${cl}_vf.nii.gz\n"
            after_script += "scil_image_math.py lower_threshold_eq ${sid}_${cl}_vf.nii.gz $params.segmentation_mask_threshold ${sid}_${cl}_mask.nii.gz --data_type uint8 -f\n"
            i += 1
        }
        """
        magic-monkey seg2mask --in $segmentation --values 1,2,3,4 --labels 01,02,04,03 --out ${segmentation.simpleName}
        scil_image_math.py addition ${segmentation.simpleName}_02.nii.gz ${segmentation.simpleName}_04.nii.gz ${segmentation.simpleName}_02.nii.gz --data_type uint8 -f
        antsAtroposN4.sh -d 3 -a $t1_image -x $mask -c ${params.segmentation_classes.size()} -p ${segmentation.simpleName}_%02d.nii.gz -o ${sid}_
        mv ${sid}_Segmentation.nii.gz ${sid}_segmentation.nii.gz
        $after_script
        """
}