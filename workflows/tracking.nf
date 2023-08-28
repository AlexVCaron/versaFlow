#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def isCollectionOrArray ( object ) {    
    return [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def asArray ( object ) {
    return isCollectionOrArray(object)
        ? object
        : object instanceof String
            ? object.tokenize(',')
            : [ object ]
}

params.pft_tracking = true
params.local_tracking = true
params.run_commit = true

params.use_opencl_tracking = false

params.pft_random_seed = 0
params.local_random_seed = 0

params.pft_tracking_algorithm = "prob"
params.local_tracking_algorithm = "prob"

params.pft_seeding_mask = "interface"
params.local_seeding_mask = "interface"
params.local_tracking_mask = "wm"

params.pft_seeding_strategy = "npv"
params.pft_number_of_seeds = 20
params.pft_step_size = 0.25
params.pft_theta_max = 20
params.pft_number_of_particles = 15
params.pft_back_tracking_length = 1

params.local_seeding_strategy = "npv"
params.local_number_of_seeds = 20
params.local_step_size = 0.25
params.local_theta_max = 20

params.pft_seeding_mask_fa_threshold = 0.1
params.local_tracking_mask_fa_threshold = 0.1
params.local_seeding_mask_fa_threshold = 0.1

pft_random_seed = asArray(params.pft_random_seed)
pft_tracking_algorithm = asArray(params.pft_tracking_algorithm)
pft_seeding_strategy = asArray(params.pft_seeding_strategy)
pft_number_of_seeds = asArray(params.pft_number_of_seeds)
pft_step_size = asArray(params.pft_step_size)
pft_theta_max = asArray(params.pft_theta_max)
pft_number_of_particles = asArray(params.pft_number_of_particles)
pft_back_tracking_length = asArray(params.pft_back_tracking_length) 

local_random_seed = asArray(params.local_random_seed)
local_tracking_algorithm = asArray(params.local_tracking_algorithm)
local_seeding_strategy = asArray(params.local_seeding_strategy)
local_number_of_seeds = asArray(params.local_number_of_seeds)
local_step_size = asArray(params.local_step_size)
local_theta_max = asArray(params.local_theta_max)

local_tracking_mask = asArray(params.local_tracking_mask)
pft_seeding_mask = asArray(params.pft_seeding_mask)
local_seeding_mask = asArray(params.local_seeding_mask)


include {
    PFT_maps;
    PFT_tracking;
    Local_tracking;
    Ensemble_Tractograms;
    Commit
} from "../modules/processes/tracking.nf"

include {
    upper_threshold_image as pft_seeding_mask_threshold_fa;
    upper_threshold_image as local_tracking_mask_threshold_fa;
    upper_threshold_image as local_seeding_mask_threshold_fa
} from "../modules/processes/utils.nf"

workflow tracking_wkf {
    take:
        fodfs
        volume_fractions
    main:
        out_tractogram = Channel.empty()
        out_maps = Channel.empty()
        out_interface = Channel.empty()

        wm_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_wm") }] }
        gm_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_gm") }] }
        csf_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_csf") }] }

        PFT_maps(wm_vf.join(gm_vf).join(csf_vf), "tracking")
        PFT_tracking(
            fodfs.join(PFT_maps.out.maps).join(PFT_maps.out.wm_gm_interface),
            "tracking",
            pft_random_seed,
            pft_tracking_algorithm,
            pft_seeding_strategy,
            pft_number_of_seeds,
            pft_step_size,
            pft_theta_max,
            pft_number_of_particles,
            pft_back_tracking_length
        )

        out_tractogram = PFT_tracking.out.tractogram
        out_maps = PFT_maps.out.maps
        out_interface = PFT_maps.out.wm_gm_interface
    emit:
        tractogram = out_tractogram
        maps = out_maps
        wm_gm_interface = out_interface
}


workflow ensemble_tracking_wkf {
    take:
        dwi_channel
        fodf_channel
        tensor_channel
        fa_channel
        wm_mask_channel
        pvf_channel
    main:
        out_tractograms = Channel.empty()
        out_ensemble_tractogram = Channel.empty()

        wm_vf_channel = pvf_channel
            .map { [it[0], it[1].find{ i -> i.simpleName.contains("_wm") }] }
        gm_vf_channel = pvf_channel
            .map { [it[0], it[1].find{ i -> i.simpleName.contains("_gm") }] }
        csf_vf_channel = pvf_channel
            .map { [it[0], it[1].find{ i -> i.simpleName.contains("_csf") }] }

        PFT_maps(
            wm_vf_channel.join(gm_vf_channel).join(csf_vf_channel),
            "tracking"
        )

        if ( params.pft_tracking ) {
            pft_seeding_masks_channel = Channel.empty()

            if ( pft_seeding_mask.contains("interface") ) {
                pft_seeding_masks_channel = pft_seeding_masks_channel.mix(
                    PFT_maps.out.wm_gm_interface.map{ it + ["wm_gm_interface"] }
                )
            }
            if ( pft_seeding_mask.contains("wm") ) {
                pft_seeding_masks_channel = pft_seeding_masks_channel.mix(
                    wm_mask_channel.map{ it + ["wm_mask"] }
                )
            }
            if ( pft_seeding_mask.contains("fa") ) {
                pft_seeding_mask_threshold_fa(
                    fa_channel,
                    params.pft_seeding_mask_fa_threshold,
                    "tracking/pft_tracking", "false"
                )
                pft_seeding_masks_channel = pft_seeding_masks_channel.mix(
                    pft_seeding_mask_threshold_fa.out.mask.map{ it + ["fa_mask"] }
                )
            }

            PFT_tracking(
                fodf_channel.join(PFT_maps.out.maps)
                    .combine(pft_seeding_masks_channel, by: 0),
                "tracking/pft_tracking",
                pft_random_seed,
                pft_tracking_algorithm,
                pft_seeding_strategy,
                pft_number_of_seeds,
                pft_step_size,
                pft_theta_max,
                pft_number_of_particles,
                pft_back_tracking_length
            )

            out_tractograms = out_tractograms.mix(
                PFT_tracking.out.tractogram
            )
        }

        if ( params.local_tracking ) {
            local_tracking_masks_channel = Channel.empty()
            local_seeding_masks_channel = Channel.empty()

            if ( local_tracking_mask.contains("wm") ) {
                local_tracking_masks_channel = local_tracking_masks_channel.mix(
                    wm_mask_channel.map{ it + ["wm_mask"] }
                )
            }
            if ( local_tracking_mask.contains("fa") ) {
                local_tracking_mask_threshold_fa(
                    fa_channel,
                    params.local_tracking_mask_fa_threshold,
                    "tracking/local_tracking", "false"
                )
                local_tracking_masks_channel = local_tracking_masks_channel.mix(
                    local_tracking_mask_threshold_fa.out.mask.map{ it + ["fa_mask"] }
                )
            }

            if ( local_seeding_mask.contains("interface") ) {
                local_seeding_masks_channel = local_seeding_masks_channel.mix(
                    PFT_maps.out.wm_gm_interface.map{ it + ["wm_gm_interface"] }
                )
            }
            if ( local_seeding_mask.contains("wm") ) {
                local_seeding_masks_channel = local_seeding_masks_channel.mix(
                    wm_mask_channel.map{ it + ["wm_mask"] }
                )
            }
            if ( local_seeding_mask.contains("fa") ) {
                local_seeding_mask_threshold_fa(
                    fa_channel,
                    params.local_seeding_mask_fa_threshold,
                    "tracking/local_tracking", "false"
                )
                local_seeding_masks_channel = local_seeding_masks_channel.mix(
                    local_seeding_mask_threshold_fa.out.mask.map{ it + ["fa_mask"] }
                )
            } 

            if ( params.use_opencl_tracking ) {
                Local_prob_tracking_opencl(
                    fodf_channel
                        .combine(local_seeding_masks_channel, by: 0)
                        .combine(local_tracking_masks_channel, by: 0),
                    "tracking/local_tracking",
                    local_random_seed,
                    local_seeding_strategy,
                    local_number_of_seeds,
                    local_step_size,
                    local_theta_max
                )

                out_tractograms = out_tractograms.mix(
                    Local_prob_tracking_opencl.out.tractogram
                )
            }
            else {
                Local_tracking(
                    fodf_channel
                        .combine(local_seeding_masks_channel, by: 0)
                        .combine(local_tracking_masks_channel, by: 0),
                    "tracking/local_tracking",
                    local_random_seed,
                    local_tracking_algorithm,
                    local_seeding_strategy,
                    local_number_of_seeds,
                    local_step_size,
                    local_theta_max
                )

                out_tractograms = out_tractograms.mix(
                    Local_tracking.out.tractogram
                )
            }
        }

        if (params.run_commit) {
            Commit(
                out_tractograms
                    .combine(dwi_channel, by: 0)
                    .map{ it + ["", ""] },
                "tracking/commit_filtered"
            )

            out_tractograms = Commit.out.filtered_tractogram
        }
        
        Ensemble_Tractograms(out_tractograms.groupTuple().join(fa_channel), "tracking")
        out_ensemble_tractogram = Ensemble_Tractograms.out.tractogram

    emit:
        tractograms = out_tractograms
        ensemble_tractogram = out_ensemble_tractogram
}
