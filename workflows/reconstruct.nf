#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.recons_dti = true
params.recons_csd = true
params.recons_diamond = true

include {
    diamond_wkf;
    dti_wkf;
    csd_wkf
} from '../modules/workflows/reconstruct.nf'

workflow reconstruct_wkf {
    take:
        dwi_channel
        mask_channel
        tissue_masks_channel
        safe_wm_mask_channel
        metadata_channel
    main:
        dti_channel = Channel.empty()
        csd_channel = Channel.empty()
        diamond_channel = Channel.empty()
        diamond_summary_channel = Channel.empty()

        if ( params.recons_dti  ) {
            dti_wkf(dwi_channel, mask_channel)

            dti_channel = dti_channel.mix(dti_wkf.out.dti)
        }

        if ( params.recons_csd  ) {
            csd_wkf(
                dwi_channel,
                mask_channel,
                tissue_masks_channel,
                safe_wm_mask_channel
            )

            csd_channel = csd_wkf.out.odfs
                .map{ [it[0], it[1..-1]] }
                .mix(csd_channel)
        }

        if ( params.recons_diamond  ) {
            diamond_wkf(dwi_channel, mask_channel)

            diamond_channel = diamond_channel.mix(diamond_wkf.out.data)
            diamond_summary_channel = diamond_summary_channel
                .mix(diamond_wkf.out.xml_summary)
        }

    emit:
        dti = dti_channel
        csd = csd_channel
        diamond = diamond_channel
        diamond_summary = diamond_summary_channel
}
