#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reconstruct_use_mrtrix = false
params.recons_dti = true
params.recons_csd = true
params.recons_diamond = true


params.config.workflow.dti_for_odf_metrics = file("$projectDir/.config/.workflow/dti_for_odf_metrics.py")
params.config.measure.diamond = file("$projectDir/.config/diamond_metrics.py")
params.config.measure.dti = file("$projectDir/.config/dti_metrics.py")

include { dti_metrics; dti_metrics as dti_for_odfs_metrics; diamond_metrics; odf_metrics; scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { uniformize_naming; replace_naming_to_underscore; rename_according_to; rename } from '../modules/functions.nf'
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
            mask_dti = rename_according_to(mask_channel, data_dti.map{ it.subList(0, 2) }, "dti_mask", false)
            prefix_dti = data_dti.map{ [it[0], "${it[0]}__dti"] }
            dti_metrics(prefix_dti.join(mask_dti).join(data_dti).join(metadata_dti), "measure", params.config.measure.dti)
            dti_channel = dti_metrics.out.metrics
        }

        if ( params.recons_diamond ) {
            data_diamond = data_channel.map{ [it[0], it[3]] }
            metadata = rename(metadata_channel, "diamond_metadata")
            mask_diamond = rename(mask_channel, "diamond_mask")
            prefix_channel = data_diamond.map{ [it[0], "${it[0]}_diamond"] }
            diamond_metrics(prefix_channel.join(mask_diamond).join(data_diamond).join(metadata), "measure", params.config.measure.diamond)
            diamond_channel = diamond_metrics.out.metrics
        }

        if ( params.recons_csd ) {
            scil_compute_dti_fa(
                dwi_channel.join(mask_channel),
                "preprocess", "measure"
            )
            data_odfs = data_channel.map{ [it[0], it[2]] }.join(scil_compute_dti_fa.out.fa).join(scil_compute_dti_fa.out.md)

            mask_odfs = uniformize_naming(mask_channel, "desc07_odf_mask", false)
            if ( !params.reconstruct_use_mrtrix )
                mask_odfs = uniformize_naming(mask_channel, "fodf_mask", false)

            basis = "tournier07"
            if ( !params.reconstruct_use_mrtrix )
                basis = "descoteaux07"

            odf_metrics(data_odfs.join(mask_odfs).filter{ !it.contains(null) }, "measure", basis)
            odfs_channel = odf_metrics.out.metrics
        }
    emit:
        metrics = dti_channel.join(diamond_channel)
        dti = dti_channel
        diamond = diamond_channel
        odfs = odfs_channel
}
