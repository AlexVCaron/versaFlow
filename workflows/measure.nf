#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { dti_metrics; dti_metrics as dti_for_odfs_metrics; diamond_metrics; odf_metrics; scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { uniformize_naming; replace_naming_to_underscore; rename_according_to; rename; get_config_path } from '../modules/functions.nf'
include { dti_wkf } from './reconstruct.nf'

params.reconstruct_use_mrtrix = false
params.recons_dti = true
params.recons_csd = true
params.recons_diamond = true
params.msmt_odf = false

params.measures_on_diamond_config = file("${get_config_path()}/measures_on_diamond_config.py")
params.measures_on_dti_config = file("${get_config_path()}/measures_on_dti_config.py")


workflow measure_wkf {
    take:
        dwi_channel
        data_channel
        mask_channel
        diamond_summary_channel
        metadata_channel
    main:
        dti_channel = Channel.empty()
        diamond_channel = Channel.empty()
        odfs_channel = Channel.empty()
        data_dti = Channel.empty()

        if ( params.recons_dti && params.reconstruct_use_mrtrix ) {
            data_dti = data_channel.map{ [it[0], it[1]] }
            metadata_dti = uniformize_naming(metadata_channel, "dti_metadata", "false", "false")
            mask_dti = rename_according_to(mask_channel, data_dti.map{ it.subList(0, 2) }, "dti_mask", false)
            prefix_dti = data_dti.map{ [it[0], "${it[0]}__dti"] }
            dti_metrics(prefix_dti.join(mask_dti).join(data_dti).join(metadata_dti), "measure", params.measures_on_dti_config)
            dti_channel = dti_metrics.out.metrics
        }

        if ( params.recons_diamond ) {
            data_diamond = data_channel.map{ [it[0], it[3]] }
            metadata = rename(metadata_channel, "diamond_metadata")
            mask_diamond = rename(mask_channel, "diamond_mask")
            prefix_channel = data_diamond.map{ [it[0], "${it[0]}_diamond"] }
            diamond_metrics(
                prefix_channel.join(mask_diamond).join(data_diamond).join(diamond_summary_channel).join(metadata),
                "measure",
                params.measures_on_diamond_config
            )
            diamond_channel = diamond_metrics.out.metrics
        }

        if ( params.recons_csd ) {
            scil_compute_dti_fa(
                dwi_channel.join(mask_channel),
                "preprocess", "measure", false
            )

            if ( params.msmt_odf )
                data_odfs = data_channel.map{ [it[0], it[2][0], it[2][2]] }
            else
                data_odfs = data_channel.map{ [it[0], it[2][0], ""] }

            data_odfs = data_odfs.join(scil_compute_dti_fa.out.fa).join(scil_compute_dti_fa.out.md)

            mask_odfs = uniformize_naming(mask_channel, "desc07_odf_mask", "false", "false")
            if ( !params.reconstruct_use_mrtrix )
                mask_odfs = uniformize_naming(mask_channel, "fodf_mask", "false", "false")

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
