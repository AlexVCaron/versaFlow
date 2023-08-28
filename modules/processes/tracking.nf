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
params.pft_min_tract_length = 5
params.pft_max_tract_length = 150
params.pft_number_of_particles = 15
params.pft_back_tracking_length = 2
params.pft_forward_tracking_length = 1

params.local_sfthres = 0.1
params.local_min_len = 5
params.local_max_len = 150

params.commit_frf_n_directions = 500
params.commit_n_iterations = 1000
params.commit_run_commit2 = false
params.commit_lambda_commit2 = 1E-3
params.commit_use_ball_stick = false
params.commit_parallel_diffusivity = 1.7E-3
params.commit_perpendicular_diffusivity = 0.51E-3
params.commit_isotropic_diffusivity = [1.7E-3, 3.0E-3]

include { remove_alg_suffixes; asArray; remove_extension } from "../functions.nf"

process PFT_maps {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${params.output_root}/${sid}/tracking", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode, overwrite: true

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
        scil_compute_maps_for_particle_filter_tracking.py \
            $wm_vf $gm_vf $csf_vf \
            --include ${sid}_map_include.nii.gz \
            --exclude ${sid}_map_exclude.nii.gz \
            --interface ${sid}_wm_gm_interface.nii.gz \
            -t $params.pve_threshold -f
        """
}

process PFT_tracking {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(fodf), path(map_include), path(map_exclude), path(seeding_mask), val(mask_type)
        val(caller_name)
        each seed
        each algo
        each seeding_strategy
        each n_seeds
        each step_length
        each theta
        each n_particles
        each back_length
    output:
        tuple val(sid), path("${sid}_pft_${algo}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${mask_type}_step_${step_length}_theta_${theta}_np_${n_particles}_back_${back_length}_tracking.trk"), emit: tractogram
    script:
        def compress = params.streamline_compression_factor ? '--compress ' + params.streamline_compression_factor : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_pft.py $fodf $seeding_mask $map_include $map_exclude \
            tmp.trk \
            --algo $algo \
            --${seeding_strategy} $n_seeds \
            --seed $seed \
            --step $step_length \
            --theta $theta \
            --sfthres $params.pft_sfthres \
            --sfthres_init $params.pft_sfthres_init \
            --min_length $params.pft_min_tract_length \
            --max_length $params.pft_max_tract_length \
            --particles $n_particles \
            --back $back_length \
            --forward $params.pft_forward_tracking_length $compress \
            --sh_basis descoteaux07 \
            --sphere symmetric724 \
            --subdivide_sphere 2
        scil_remove_invalid_streamlines.py tmp.trk\
            ${sid}_pft_${algo}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${mask_type}_step_${step_length}_theta_${theta}_np_${n_particles}_back_${back_length}_tracking.trk \
            --remove_single_point
        """
}

process Local_prob_tracking_opencl {
    label params.use_cuda ? "res_single_cpu" : params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"
    label params.use_cuda ? "res_gpu" : ""

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(fodf), path(seeding_mask, stageAs: "seeding_mask*.nii.gz"), val(seeding_mask_type), path(tracking_mask, stageAs: "tracking_mask*.nii.gz"), val(tracking_mask_type)
        val(caller_name)
        each seed
        each seeding_strategy
        each n_seeds
        each step_length
        each theta
    output:
        tuple val(sid), path("${sid}_local_gpu_prob_in_${tracking_mask_type}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${seeding_mask_type}_step_${step_length}_theta_${theta}_tracking.trk"), emit: tractogram
    script:
        def compress = params.streamline_compression_factor ? '--compress ' + params.streamline_compression_factor : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_local_tracking_gpu.py $fodf $seeding_mask $tracking_mask \
            tmp.trk \
            --step $step_length \
            --theta $theta \
            --min_length $params.local_min_len \
            --max_length $params.local_max_len \
            --${seeding_strategy} $n_seeds $compress \
            --rng_seed $seed
        scil_remove_invalid_streamlines.py tmp.trk \
            ${sid}_local_gpu_prob_in_${tracking_mask_type}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${seeding_mask_type}_step_${step_length}_theta_${theta}_tracking.trk \
            --remove_single_point
        """
}

process Local_tracking {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(fodf), path(seeding_mask, stageAs: "seeding_mask*.nii.gz"), val(seeding_mask_type), path(tracking_mask, stageAs: "tracking_mask*.nii.gz"), val(tracking_mask_type)
        val(caller_name)
        each seed
        each algo
        each seeding_strategy
        each n_seeds
        each step_length
        each theta
    output:
        tuple val(sid), path("${sid}_local_${algo}_in_${tracking_mask_type}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${seeding_mask_type}_step_${step_length}_theta_${theta}_tracking.trk"), emit: tractogram
    script:
        def compress = params.streamline_compression_factor ? '--compress ' + params.streamline_compression_factor : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_local_tracking.py $fodf $seeding_mask $tracking_mask \
            tmp.trk \
            --algo $algo \
            --${seeding_strategy} $n_seeds \
            --seed $seed \
            --step $step_length \
            --theta $theta \
            --sfthres $params.local_sfthres \
            --min_length $params.local_min_len \
            --max_length $params.local_max_len $compress \
            --sh_basis descoteaux07 \
            --sphere symmetric724 \
            --subdivide_sphere 2
        scil_remove_invalid_streamlines.py tmp.trk \
            ${sid}_local_${algo}_in_${tracking_mask_type}_seed_${seed}_${seeding_strategy}${n_seeds}_in_${seeding_mask_type}_step_${step_length}_theta_${theta}_tracking.trk \
            --remove_single_point
        """
}

process Ensemble_Tractograms {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(tractograms), path(reference)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}_ensemble_tractogram.trk"), emit: tractogram
    script:
        """
        scil_tractogram_math.py union \
            ${tractograms.join(" ")} \
            ${sid}_ensemble_tractogram.trk \
            --reference $reference
        """
}

process Commit {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "link", enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: {
        f -> f.contains(".trk")
            ? "tractograms/${remove_alg_suffixes(f)}"
            : f.contains(".txt")
                ? "weights/${remove_alg_suffixes(f)}"
                : "fractions/${remove_alg_suffixes(f)}"
    }, mode: params.publish_mode

    input:
        tuple val(sid), path(tractogram), path(dwi), path(bval), path(bvec), file(tracking_mask), file(peaks)
        val(caller_name)
    output:
        tuple val(sid), path("${remove_extension(tractogram, 1)}_filtered.trk"), emit: filtered_tractogram
        tuple val(sid), path("${remove_extension(tractogram, 1)}_outliers.trk"), emit: outlier_streamlines
        tuple val(sid), path("${remove_extension(tractogram, 1)}_EC_fraction.nii.gz"), emit: extracellular_fraction
        tuple val(sid), path("${remove_extension(tractogram, 1)}_IC_fraction.nii.gz"), emit: intracellular_fraction
        tuple val(sid), path("${remove_extension(tractogram, 1)}_FW_fraction.nii.gz"), emit: freewater_fraction
        tuple val(sid), path("${remove_extension(tractogram, 1)}_commit_weights.txt"), emit: commit_weights
        tuple val(sid), path("${remove_extension(tractogram, 1)}_total_commit_weights.txt"), emit: total_commit_weights
        tuple val(sid), path("${remove_extension(tractogram, 1)}_streamlines_length.txt"), emit: streamlines_length
    script:
        def opt_params = ""
        if ( params.commit_run_commit2 ) opt_params += " --commit2"
        if ( params.commit_use_ball_stick ) opt_params += " --ball_stick"
        if ( !peaks.empty() ) opt_params += " --in_peaks $peaks"
        if ( !tracking_mask.empty() ) opt_params += " --in_tracking_mask $tracking_mask"
        if ( params.commit_perpendicular_diffusivity ) opt_params += " --perp_diff ${asArray(params.commit_perpendicular_diffusivity).join(" ")}"
        """
        scil_run_commit.py $tractogram $dwi $bval $bvec results \
            --b_thr $params.b0_threshold \
            --nbr_dir $params.commit_frf_n_directions \
            --nbr_iter $params.commit_n_iterations \
            --lambda_commit_2 $params.commit_lambda_commit2 \
            --para_diff $params.commit_parallel_diffusivity \
            --iso_diff ${asArray(params.commit_isotropic_diffusivity).join(" ")} \
            --processes $task.cpus $opt_params

        cp results/commit_1/essential_tractogram.trk ${remove_extension(tractogram, 1)}_filtered.trk
        cp results/commit_1/nonessential_tractogram.trk ${remove_extension(tractogram, 1)}_outliers.trk
        cp results/commit_1/compartment_EC.nii.gz ${remove_extension(tractogram, 1)}_EC_fraction.nii.gz
        cp results/commit_1/compartment_IC.nii.gz ${remove_extension(tractogram, 1)}_IC_fraction.nii.gz
        cp results/commit_1/compartment_ISO.nii.gz ${remove_extension(tractogram, 1)}_FW_fraction.nii.gz
        cp results/commit_1/streamline_weights.txt ${remove_extension(tractogram, 1)}_commit_weights.txt
        cp results/commit_1/streamline_weights_by_length.txt ${remove_extension(tractogram, 1)}_total_commit_weights.txt
        cp results/commit_1/streamlines_length.txt ${remove_extension(tractogram, 1)}_streamlines_length.txt
        """
}
