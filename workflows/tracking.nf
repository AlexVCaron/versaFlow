#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.pft_tracking = true
params.pft_random_seed = 0
params.tracking_algorithm = "prob"

pft_random_seed = params.pft_random_seed instanceof String ? params.pft_random_seed.tokenize(',') : params.pft_random_seed
tracking_algorithm = params.tracking_algorithm instanceof String ? params.tracking_algorithm.tokenize(',') : params.tracking_algorithm

include { PFT_maps; PFT_tracking } from "../modules/processes/tracking.nf"

workflow tracking_wkf {
    take:
        fodfs
        volume_fractions
    main:
        if ( params.pft_tracking ) {
            wm_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_wm") }] }
            gm_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_gm") }] }
            csf_vf = volume_fractions.map { [it[0], it[1].find{ i -> i.simpleName.contains("_csf") }] }

            PFT_maps(wm_vf.join(gm_vf).join(csf_vf), "tracking")
            PFT_tracking(
                fodfs.join(PFT_maps.out.maps).join(PFT_maps.out.wm_gm_interface),
                "tracking",
                pft_random_seed,
                tracking_algorithm
            )
        }
    emit:
        tractogram = PFT_tracking.out.tractogram
        maps = PFT_maps.out.maps
        wm_gm_interface = PFT_maps.out.wm_gm_interface
}
