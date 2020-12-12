#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reconstruct_use_mrtrix = false
params.recons_dti = true
params.recons_csd = true
params.recons_diamond = true


params.config.workflow.dti_for_odf_metrics = file("$projectDir/.config/.workflow/dti_for_odf_metrics.py")
params.config.measure.diamond = file("$projectDir/.config/diamond_metrics.py")
params.config.measure.dti = file("$projectDir/.config/dti_metrics.py")

include { dti_metrics; dti_metrics as dti_for_odfs_metrics; diamond_metrics; odf_metrics } from '../modules/processes/measure.nf'
include { uniformize_naming; replace_naming_to_underscore } from '../modules/functions.nf'
include { dti_wkf } from './reconstruct.nf'

workflow measure_wkf {
    take:
        dwi_channel
        data_channel
        mask_channel
        metadata_channel
    main:
        dti_channel = Channel.empty()
        diamond_channel = Channel.empty()
        odfs_channel = Channel.empty()
        data_dti = Channel.empty()

        if ( params.recons_dti && params.reconstruct_use_mrtrix ) {
            data_dti = data_channel.map{ [it[0], it[1]] }
            metadata_dti = uniformize_naming(metadata_channel, "dti_metadata", false)
            mask_dti = uniformize_naming(mask_channel, "dti_mask", false)
            prefix_dti = data_dti.map{ [it[0], "${it[0]}__dti"] }
            dti_metrics(prefix_dti.join(mask_dti).join(data_dti).join(metadata_dti), "measure", params.config.measure.dti)
            dti_channel = dti_metrics.out.metrics
        }

        if ( params.recons_diamond ) {
            data = data_channel.map{ [it[0], it[3]] }
            metadata = uniformize_naming(metadata_channel, "diamond_metadata", false)
            mask_diamond = uniformize_naming(mask_channel, "diamond_mask", false)
            prefix_channel = data.map{ [it[0], "${it[0]}__diamond"] }
            diamond_metrics(prefix_channel.join(mask_diamond).join(data).join(metadata), "measure", params.config.measure.diamond)
            diamond_channel = diamond_metrics.out.metrics
        }

        if ( params.recons_csd ) {
            dti_wkf(dwi_channel, mask_channel)
            data_dti = dti_wkf.out.dti

            metadata_dti = uniformize_naming(metadata_channel, "dti_metadata", false)
            mask_dti = uniformize_naming(mask_channel, "dti_mask", false)
            prefix_dti = data_dti.map{ [it[0], "${it[0]}__dti"] }
            dti_for_odfs_metrics(
                prefix_dti.join(mask_dti).join(data_dti).join(metadata_dti),
                "measure", params.config.workflow.dti_for_odf_metrics
            )
            data = data_channel.map{ [it[0], it[2]] }.join(dti_for_odfs_metrics.out.metrics).map{ it.flatten() }

            mask_odfs = uniformize_naming(mask_channel, "desc07_odf_mask", false)
            if ( !params.reconstruct_use_mrtrix )
                mask_odfs = uniformize_naming(mask_channel, "fodf_mask", false)

            basis = "tournier07"
            if ( !params.reconstruct_use_mrtrix )
                basis = "descoteaux07"

            odf_metrics(data.join(mask_odfs), "measure", basis)
            odfs_channel = odf_metrics.out.metrics
        }
    emit:
        metrics = dti_channel.join(diamond_channel)
        dti = dti_channel
        diamond = diamond_channel
        odfs = odfs_channel
}
