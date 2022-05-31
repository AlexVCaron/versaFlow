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
        display_run_info()
        dataloader = load_dataset()
        preprocess_wkf(dataloader.dwi, dataloader.rev, dataloader.anat, dataloader.pvf, dataloader.metadata, dataloader.rev_metadata, dataloader.dwi_mask, dataloader.anat_mask)
        reconstruct_wkf(preprocess_wkf.out.dwi, preprocess_wkf.out.mask, preprocess_wkf.out.pvf, preprocess_wkf.out.safe_wm_mask, preprocess_wkf.out.metadata)
        measure_wkf(preprocess_wkf.out.dwi, reconstruct_wkf.out.all, preprocess_wkf.out.mask, preprocess_wkf.out.diamond_summary, preprocess_wkf.out.metadata)
        if ( params.pft_tracking ) {
            tracking_wkf(reconstruct_wkf.out.csd, preprocess_wkf.out.pvf)
        }
    }
}

def validate_required_parameters () {
    if ( !params.data_root ) error "Error ~ Input data root not specified, use --data_root"
    if ( params.resample_data && !params.resampling_resolution ) error "Error ~ Resampling is enabled, but resampling resolution is not defined, use --resampling_resolution"
    if ( params.pft_tracking && !params.recons_csd ) error "Error ~ CSD reconstruction is required for tracking step, use --recons_csd to enable fodf reconstruction, or disable tracking with --pft_tracking = false"
}

def display_run_info () {
    log.info "mrHARDIflow pipeline"
    log.info "==================="
    log.info ""
    log.info "Start time : $workflow.start"
    log.info ""
    log.info "Run parameters"
    log.info "=============="
    log.info ""
    log.info "I/O :"
    log.info " - Input root      : $params.data_root"
    log.info "    - b0 threshold : ${params.b0_threshold ? params.b0_threshold : 0}"
    log.info " - Output root     : $params.output_root"
    log.info " - Publish mode    : $params.publish_mode"
    log.info " - Publish all outputs ${params.publish_all ? "(enabled)" : "(disabled)"}"
    log.info " - Verbose ${params.verbose_outputs ? "(enabled)" : "(disabled)"}"
    log.info "Resources allocation :"
    log.info " - Use GPU ${params.use_cuda ? "(enabled)" : "(disabled)"}"
    if (params.use_cuda) {
        log.info "    - Auto-select GPU ${params.eddy_select_gpu ? "(enabled)" : "(disabled)"}"
        log.info "    - Max parallel GPU : $params.cuda_max_parallel"
    }
    log.info " - Max CPU per process : ${params.max_cpu_per_process ? params.max_cpu_per_process : params.processes}"
    log.info ""
    log.info "DWI preprocessing :"
    log.info " - Background denoising ${params.gaussian_noise_correction ? "(enabled)" : "(disabled)"}"
    log.info " - Gibbs ringing ${params.gibbs_ringing_correction ? "(enabled)" : "(disabled)"}"
    log.info " - B0 normalization ${params.normalize_inter_b0 ? "(enabled)" : "(disabled)"}"
    if (params.normalize_inter_b0) {
        log.info "    - Normalization strategy : $params.b0_normalization_strategy"
    }
    log.info " - Topup ${params.topup_correction ? "(enabled)" : "(disabled)"}"
    log.info " - Eddy ${params.eddy_correction ? "(enabled)" : "(disabled)"}"
    if (params.eddy_correction) {
        log.info "    - DWI shells check ${params.eddy_force_shelled ? "(disabled)" : "(enabled)"}"
        log.info "    - Use reverse phase ${params.eddy_with_reverse ? "(enabled)" : "(disabled)"}"
    }
    log.info " - N4 normalization ${params.dwi_intensity_normalization ? "(enabled)" : "(disabled)"}"
    log.info "T1 preprocessing :"
    log.info " - Background denoising ${params.denoise_t1 ? "(enabled)" : "(disabled)"}"
    if (params.denoise_t1) {
        log.info "    - Use nlmeans ${params.nlmeans_t1 ? "(enabled)" : "(disabled)"}"
    }
    log.info " - N4 normalization ${params.t1_intensity_normalization ? "(enabled)" : "(disabled)"}"
    log.info "T1 to DWI registration :"
    log.info " - Register T1 mask ${params.t1mask2dwi_registration ? "(enabled)" : "(disabled)"}"
    log.info " - Register after denoising ${params.register_t12b0_denoised ? "(enabled)" : "(disabled)"}"
    if (params.t1mask2dwi_registration || params.register_t12b0_denoised) {
        log.info " - Use SyN deformation ${params.register_syn_t12b0 ? "(enabled)" : "(disabled)"}"
        if (params.register_syn_t12b0) {
            log.info "    - Use masked images ${params.register_syn_t12b0_with_mask ? "(enabled)" : "(disabled)"}"
        }
    }
    log.info "Upscaling :"
    log.info " - Resample T1 and DWI ${params.resample_data ? "(enabled)" : "(disabled)"}"
    if (params.resample_data) {
        log.info "    - Sequential processing ${params.force_resampling_sequential ? "(enabled)" : "(disabled)"}"
        log.info "    - Resampling resolution : $params.resampling_resolution"
    }
    log.info "Segmentation :"
    log.info " - Segment WM/GM/CSF from T1 ${params.generate_tissue_segmentation ? "(enabled)" : "(disabled)"}"
    if (params.generate_tissue_segmentation) {
        log.info "     - Template directory : ${params.tissue_segmentation_root}"
    }
    log.info " - Segment WM parcellation ${params.generate_wm_segmentation ? "(enabled)" : "(disabled)"}"
    if (params.generate_wm_segmentation) {
        log.info "     - Atlas directory : ${params.wm_segmentation_root}"
    }
    log.info "Diffusion modeling :"
    log.info " - DTI ${params.recons_dti ? "(enabled)" : "(disabled)"}"
    if (params.recons_dti) {
        log.info "    - Maximal b-value : $params.max_dti_bvalue"
    }
    if (params.msmt_odf) log.info " - MSMT CSD ${params.recons_csd ? "(enabled)" : "(disabled)"}"
    else log.info " - SSST CSD ${params.recons_csd ? "(enabled)" : "(disabled)"}"
    if (params.recons_csd) {
        log.info "    - Compute FRF on DTI shells ${params.frf_on_dti_shell ? "(enabled)" : "(disabled)"}"
        log.info "    - Spherical harmonics order           : $params.sh_order"
        log.info "    - FRF - FA upper threshold            : $params.frf_fa"
        log.info "    - FRF - FA lower threshold            : $params.frf_min_fa"
        log.info "    - FRF - Minimal sample size           : $params.frf_min_nvox"
        log.info "    - FRF - Search radius                 : ${params.frf_radii ? params.frf_radii: "none"}"
        log.info "    - FRF - Search center                 : ${params.frf_center ? params.frf_center: "none"}"
        log.info "    - Ventricles - FA upper threshold     : $params.max_fa_ventricle"
        log.info "    - Ventricles - MD lower threshold     : $params.min_md_ventricle"
        log.info "    - Ventricles - Search center          : ${params.ventricles_center ? params.frf_center: "none"}"
        log.info "    - Peak filtering - Absolute factor    : $params.fodf_max_absolute_factor"
        log.info "    - Peak filtering - Relative threshold : $params.fodf_relative_thr"
    }
    log.info " - DIAMOND ${params.recons_diamond ? "(enabled)" : "(disabled)"}"
    if (params.recons_diamond) {
        log.info "    - Model selection on tensor ${params.model_selection_with_tensor ? "(enabled)" : "(disabled)"}"
        log.info "    - Estimate ISOR - restricted isotropic diffusion fraction ${params.estimate_restriction ? "(enabled)" : "(disabled)"}"
        log.info "    - Estimate ICVF - hindered diffusion fraction per fascicle ${params.estimate_hindered ? "(enabled)" : "(disabled)"}"
        log.info "    - Use FW tensor model ${params.free_water_tensor ? "(enabled)" : "(disabled)"}"
        log.info "    - Use ISOR tensor model ${params.restriction_tensor ? "(enabled)" : "(disabled)"}"
        log.info "    - Normalize fractions ${params.normalized_fractions ? "(enabled)" : "(disabled)"}"
        log.info "    - Optimize # parameters ${params.strict_n_parameters ? "(disabled)" : "(enabled)"}"
        log.info "    - Number of fascicles : $params.n_fascicles"
        log.info "    - Fascicle model      : $params.fascicle_model"
    }

    workflow.onComplete {
        log.info "Pipeline completed at : $workflow.complete"
        log.info "Execution status : ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Execution duration : $workflow.duration"
    }
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
            "frf_on_dti_shell": "$params.frf_on_dti_shell",
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
            "estimate_hindered": "$params.estimate_hindered",
            "normalized_fractions" : "$params.normalized_fractions",
            "free_water_tensor" : "$params.free_water_tensor",
            "strict_n_parameters": "$params.strict_n_parameters",
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