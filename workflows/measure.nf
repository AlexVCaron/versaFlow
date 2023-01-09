#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    diamond_metrics;
    odf_metrics;
    scil_compute_dti_fa
} from '../modules/processes/measure.nf'
include {
    change_name as rename_metadata_for_diamond;
    change_name as rename_mask_for_diamond
} from '../modules/processes/io.nf'
include {
    get_config_path;
    collect_paths
} from '../modules/functions.nf'

params.measures_on_diamond_config = file("${get_config_path()}/measures_on_diamond_config.py")


workflow measure_wkf {
    take:
        dwi_channel
        mask_channel
        tissue_masks_channel
        odfs_channel
        diamond_channel
        diamond_summary_channel
        metadata_channel
    main:
        diamond_ids_channel = diamond_channel.map{ [it[0]] }

        diamond_metadata_channel = rename_metadata_for_diamond(
            collect_paths(diamond_ids_channel.join(metadata_channel))
                .filter{ it[1] },
            "diamond_metadata"
        ).map{ it.flatten() }
        diamond_mask_channel = rename_mask_for_diamond(
            collect_paths(diamond_ids_channel.join(mask_channel)),
            "diamond_mask"
        )

        diamond_metrics(
            diamond_ids_channel
                .map{ [it[0], "${it[0]}_diamond"] }
                .join(diamond_mask_channel)
                .join(diamond_channel)
                .join(diamond_summary_channel)
                .join(diamond_metadata_channel),
            "measure",
            params.measures_on_diamond_config
        )

        odfs_ids_channel = odfs_channel.map{ [it[0]] }

        scil_compute_dti_fa(
            odfs_ids_channel
                .join(dwi_channel)
                .join(mask_channel),
            "preprocess", "measure",
            false
        )

        odf_metrics(
            odfs_channel
                .map{ it.flatten() }
                .map{ it.size() > 2 ? it : it + ["", ""] }
                .join(scil_compute_dti_fa.out.fa)
                .join(scil_compute_dti_fa.out.md)
                .join(tissue_masks_channel.map{ [it[0]] + it[1][0..-2] })
                .filter{ !it.contains(null) },
            "measure",
            "descoteaux07"
        )

    emit:
        diamond = diamond_metrics.out.metrics
        odfs = odf_metrics.out.metrics
}
