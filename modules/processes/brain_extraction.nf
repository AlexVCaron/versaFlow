#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.deep_bet_model = "Site-All-T-epoch_36_update_with_Site_6_plus_7-epoch_09.model"

process deepbet_t1 {
    label "BET"
    label params.use_cuda ? "res_single_cpu" : params.on_hcp ? "res_full_node_override" : "res_max_cpu"
    label params.use_cuda ? "res_gpu" : ""

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(t1_image)
        val(caller_name)
    output:
        tuple val(sid), path("${t1_image.simpleName}_mask.nii.gz"), emit: mask
    script:
        """
        muSkullStrip.py -in $t1_image -model ${params.deep_bet_model} -suffix mask
        """

}

process bet_mask {
    label "BET"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), suffix ? "$suffix" : "_bet_mask") }, mode: params.publish_mode

    input:
        tuple val(sid), path(img)
        val(caller_name)
        val(suffix)
    output:
        tuple val(sid), path("${img.simpleName}_bet_mask.nii.gz"), emit: mask
    script:
        """
        fslmaths $img -Tmean mean_image.nii.gz
        bet mean_image.nii.gz "${img.simpleName}_bet.nii.gz" -m -R -f $params.bet.f
        """
}