#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false

//include { print_channel } from "./modules/debug.nf" // For debugging purpose only
include { load_dataset } from "./workflows/io.nf"
include { preprocess_wkf } from "./workflows/preprocess.nf"
include { reconstruct_wkf } from "./workflows/reconstruct.nf"
include { measure_wkf } from "./workflows/measure.nf"
include { tracking_wkf } from "./workflows/tracking.nf"

workflow {
    if (params.help) display_usage()
    else {
        validate_required_parameters()
        dataloader = load_dataset()
        preprocess_wkf(dataloader.dwi, dataloader.rev, dataloader.anat, dataloader.seg, dataloader.metadata, dataloader.rev_metadata, dataloader.dwi_mask, dataloader.anat_mask)
        reconstruct_wkf(preprocess_wkf.out.dwi, preprocess_wkf.out.mask, preprocess_wkf.out.segmentation, preprocess_wkf.out.metadata)
        measure_wkf(preprocess_wkf.out.dwi, reconstruct_wkf.out.all, preprocess_wkf.out.mask, preprocess_wkf.out.metadata)
        tracking_wkf(reconstruct_wkf.out.csd, preprocess_wkf.out.segmentation)
    }
}

def validate_required_parameters () {
    if ( !params.data_root ) error "Error ~ Input data root not specified, use --data_root"
    if ( params.resample_data && !params.resampling_resolution ) error "Error ~ Resampling is enabled, but resampling resolution is not defined, use --resampling_resolution"
    if ( params.pft_tracking && !params.recons_csd ) error "Error ~ CSD reconstruction is required for tracking step, use --recons_csd to enable fodf reconstruction, or disable tracking with --pft_tracking = false"
}

def display_usage () {
    usage = file("$projectDir/USAGE")

    bindings = [
            "use_cuda" : "$params.use_cuda",
            "conservative_resources" : "$params.conservative_resources",
            "free_processes" : "$params.free_processes",
            "max_cpu_per_process" : "$params.max_cpu_per_process",
            "max_attempts" : "$params.max_attempts",
            "output_root" : "$params.output_root",
            "publish_all" : "$params.publish_all",
            "publish_mode" : "$params.publish_mode",
            "verbose_outputs" : "$params.verbose_outputs",
            "resample_data" : "$params.resample_data",
            "force_resampling_sequential" : "$params.force_resampling_sequential",
            "b0_threshold" : "$params.b0_threshold",
            "b0_normalization_strategy" : "$params.b0_normalization_strategy",
            "bet_f" : "$params.bet_f",
            "segmentation_mask_threshold" : "$params.segmentation_mask_threshold ",
            "t1_intensity_normalization" : "$params.t1_intensity_normalization",
            "t1mask2dwi_registration" : "$params.t1mask2dwi_registration",
            "register_t12b0_denoised" : "$params.register_t12b0_denoised",
            "register_syn_t12b0" : "$params.register_syn_t12b0",
            "register_syn_t12b0_with_mask" : "$params.register_syn_t12b0_with_mask",
            "denoise_t1" : "$params.denoise_t1",
            "nlmeans_t1" : "$params.nlmeans_t1",
            "generate_tissue_segmentation" : "$params.generate_tissue_segmentation",
            "generate_wm_segmentation" : "$params.generate_wm_segmentation",
            "gaussian_noise_correction" : "$params.gaussian_noise_correction",
            "gibbs_ringing_correction" : "$params.gibbs_ringing_correction",
            "normalize_inter_b0" : "$params.normalize_inter_b0",
            "topup_correction" : "$params.topup_correction",
            "eddy_correction" : "$params.eddy_correction",
            "eddy_force_shelled" : "$params.eddy_force_shelled",
            "eddy_with_reverse" : "$params.eddy_with_reverse",
            "eddy_select_gpu" : "$params.eddy_select_gpu",
            "dwi_intensity_normalization" : "$params.dwi_intensity_normalization",
            "reconstruct_use_mrtrix" : "$params.reconstruct_use_mrtrix",
            "recons_dti" : "$params.recons_dti",
            "recons_csd" : "$params.recons_csd",
            "max_dti_bvalue": "$params.max_dti_bvalue",
            "msmt_odf" : "$params.msmt_odf",
            "convert_tournier2descoteaux" : "$params.convert_tournier2descoteaux",
            "frf_fa" : "$params.frf_fa",
            "frf_min_fa" : "$params.frf_min_fa",
            "frf_min_nvox" : "$params.frf_min_nvox",
            "frf_radii" : "$params.frf_radii",
            "frf_center" : "$params.frf_center",
            "max_fa_ventricle" : "$params.max_fa_ventricle",
            "min_md_ventricle" : "$params.min_md_ventricle",
            "ventricles_center" : "$params.ventricles_center",
            "fodf_max_absolute_factor" : "$params.fodf_max_absolute_factor",
            "fodf_relative_thr" : "$params.fodf_relative_thr",
            "sh_order" : "$params.sh_order",
            "recons_diamond" : "$params.recons_diamond",
            "n_fascicles" : "$params.n_fascicles",
            "fascicle_model" : "$params.fascicle_model",
            "model_selection_with_tensor" : "$params.model_selection_with_tensor",
            "estimate_restriction" : "$params.estimate_restriction",
            "normalized_fractions" : "$params.normalized_fractions",
            "free_water_tensor" : "$params.free_water_tensor",
            "pft_random_seed" : "$params.pft_random_seed",
            "tracking_algorithm" : "$params.tracking_algorithm",
            "streamline_compression_factor" : "$params.streamline_compression_factor",
            "pft_tracking": "$params.pft_tracking",
            "pve_threshold": "$params.pve_threshold",
            "pft_seeding_strategy" : "$params.pft_seeding_strategy",
            "pft_number_of_seeds" : "$params.pft_number_of_seeds",
            "pft_step_size" : "$params.pft_step_size",
            "pft_theta_max" : "$params.pft_theta_max",
            "pft_sfthres" : "$params.pft_sfthres",
            "pft_sfthres_init" : "$params.pft_sfthres_init",
            "pft_min_tract_length" : "$params.pft_min_tract_length",
            "pft_max_tract_length" : "$params.pft_max_tract_length",
            "pft_number_of_particles" : "$params.pft_number_of_particles",
            "pft_back_tracking_length" : "$params.pft_back_tracking_length",
            "pft_forward_tracking_length" : "$params.pft_forward_tracking_length",
            "raw_to_processed_space": "$params.raw_to_processed_space",
            "cuda_max_parallel": "$params.cuda_max_parallel"
    ]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
}