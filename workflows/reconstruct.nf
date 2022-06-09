#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.recons_dti = true
params.recons_csd = true
params.recons_diamond = true

include { diamond_wkf; dti_wkf; csd_wkf } from '../modules/workflows/reconstruct.nf'

workflow reconstruct_wkf {
    take:
        dwi_channel
        mask_channel
        pvf_channel
        safe_wm_mask_channel
        metadata_channel
    main:
        out_channels = []

        if ( params.recons_dti  ) {
            dti_wkf(dwi_channel, mask_channel)
            out_channels += [dti_wkf.out.dti]
        }
        else
            out_channels += [null]

        if ( params.recons_csd  ) {
            csd_wkf(dwi_channel, mask_channel, pvf_channel, safe_wm_mask_channel)
            out_channels += [csd_wkf.out.odfs.map{ [it[0], it.subList(1, it.size())] }]
        }
        else
            out_channels += [null]

        if ( params.recons_diamond  ) {
            diamond_wkf(dwi_channel, mask_channel)
            out_channels += [diamond_wkf.out.data, diamond_wkf.out.xml_summary]
        }
        else
            out_channels += [null, null]

    emit:
        dti = out_channels[0]
        csd = out_channels[1]
        diamond = out_channels[2]
        diamond_summary = out_channels[3]
        all = out_channels[1..-1].inject(out_channels[0]){ c, n -> n ? c.join(n, remainder: true) : c }
}
