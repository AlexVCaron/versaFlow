#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.streamline_compression_factor = 0.2
params.pve_threshold = 0.05
params.pft_seeding_strategy = "npv"
params.pft_number_of_seeds = 10
params.pft_step_size = 0.5
params.pft_theta_max = 20
params.pft_sfthres = 0.1
params.pft_sfthres_init = 0.5
params.pft_min_tract_length = 20
params.pft_max_tract_length = 200
params.pft_number_of_particles = 15
params.pft_back_tracking_length = 2
params.pft_forward_tracking_length = 1

include { remove_alg_suffixes } from "../functions.nf"

process PFT_maps {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/tracking", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(wm_vf), path(gm_vf), path(csf_vf)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_map_include.nii.gz"), path("${sid}_map_exclude.nii.gz"), emit: maps
        tuple val(sid), path("${sid}_wm_gm_interface.nii.gz"), emit: wm_gm_interface
    script:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_maps_for_particle_filter_tracking.py $wm_vf $gm_vf $csf_vf --include ${sid}_map_include.nii.gz --exclude ${sid}_map_exclude.nii.gz --interface ${sid}_wm_gm_interface.nii.gz -t $params.pve_threshold -f
        """
}

process PFT_tracking {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/tracking", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(fodf), path(map_include), path(map_exclude), path(seeding_mask)
        val(caller_name)
        each seed
        each algo
    output:
        tuple val(sid), path("${sid}_pft_${algo}_seed_${seed}_tracking.trk"), emit: tractogram
    script:
        def compress = params.streamline_compression_factor ? '--compress ' + params.streamline_compression_factor : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_pft.py $fodf $seeding_mask $map_include $map_exclude ${sid}_pft_${algo}_seed_${seed}_tracking.trk --algo $algo --${params.pft_seeding_strategy} $params.pft_number_of_seeds --seed $seed --step $params.pft_step_size --theta $params.pft_theta_max --sfthres $params.pft_sfthres --sfthres_init $params.pft_sfthres_init --min_length $params.pft_min_tract_length --max_length $params.pft_max_tract_length --particles $params.pft_number_of_particles --back $params.pft_back_tracking_length --forward $params.pft_forward_tracking_length $compress --sh_basis descoteaux07
        """
}