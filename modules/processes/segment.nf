#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.classes = ["csf", "gm", "wm"]

include { remove_alg_suffixes } from '../functions.nf'

process atropos {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(t1_image), path(mask), path(priors)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_segmentation.nii.gz"), emit: segmentation
        tuple val(sid), path("${sid}_{${params.classes.join(',')}}_mask.nii.gz"), emit: masks
    script:
        mask_renaming = ""
        i = 1
        for (cl in params.classes) {
            mask_renaming += "mv ${sid}_SegmentationPosteriors0${i}.nii.gz ${sid}_${cl}_mask.nii.gz\n"
        }
        """
        antsAtroposN4.sh -d 3 -a $t1_image -x $mask -c ${params.classes.size()} -p $priors -o ${sid}_
        mv ${sid}_Segmentation.nii.gz ${sid}_segmentation.nii.gz
        $mask_renaming
        """
}
