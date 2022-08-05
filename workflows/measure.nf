#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { dti_metrics; dti_metrics as dti_for_odfs_metrics; diamond_metrics; odf_metrics; scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { uniformize_naming; replace_naming_to_underscore; rename_according_to; rename; get_config_path } from '../modules/functions.nf'
include { dti_wkf } from './reconstruct.nf'

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
        tissue_masks_channel
        diamond_summary_channel
        metadata_channel
    main:
        diamond_channel = Channel.empty()
        odfs_channel = Channel.empty()

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
            scil_compute_dti_fa(dwi_channel.join(mask_channel), "preprocess", "measure", false)

            if ( params.msmt_odf )
                data_odfs = data_channel.map{ [it[0]] + it[2] }
            else
                data_odfs = data_channel.map{ [it[0], it[2][0], "", ""] }

            data_odfs = data_odfs.join(scil_compute_dti_fa.out.fa).join(scil_compute_dti_fa.out.md)
            odf_metrics(data_odfs.join(tissue_masks_channel.map{ [it[0]] + it[1][0..-2] }).filter{ !it.contains(null) }, "measure", "descoteaux07")
            odfs_channel = odf_metrics.out.metrics
        }
    emit:
        diamond = diamond_channel
        odfs = odfs_channel
}
