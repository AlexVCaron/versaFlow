#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reconstruct_use_mrtrix = false
params.convert_tournier2descoteaux = true

include { diamond; mrtrix_dti; csd; response; scilpy_response; scilpy_csd } from '../processes/reconstruct.nf'
include { scil_dti_and_metrics } from '../processes/measure.nf'
include { tournier2descoteaux_odf } from '../processes/utils.nf'

workflow csd_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        response_channel = Channel.empty()
        odfs_channel = Channel.empty()
        if ( params.reconstruct_use_mrtrix ) {
            response(dwi_channel.join(mask_channel), "reconstruct")
            csd(response.out.responses.join(dwi_channel.join(mask_channel)), "reconstruct")
            response_channel = response.out.responses
            odfs_channel = csd.out.odfs
            if ( params.convert_tournier2descoteaux ) {
                tournier2descoteaux_odf(csd.out.odfs, "reconstruct")
                csd_channel = tournier2descoteaux_odf.out.odfs
            }
        }
        else {
            scilpy_response(dwi_channel.join(mask_channel), "reconstruct")
            scilpy_csd(dwi_channel.join(scilpy_response.out.response).join(mask_channel), "reconstruct")
            response_channel = scilpy_response.out.response
            odfs_channel = scilpy_csd.out.odfs
        }
    emit:
        odfs = odfs_channel
        responses = response_channel
}

workflow dti_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        dti_output = Channel.empty()
        if ( params.reconstruct_use_mrtrix ) {
            mrtrix_dti(dwi_channel.join(mask_channel), "reconstruct")
            dti_output = mrtrix_dti.out.dti
        }
        else {
            scil_dti_and_metrics(dwi_channel.join(mask_channel), "reconstruct", "measure")
            dti_output = scil_dti_and_metrics.out.dti
        }
    emit:
        dti = dti_output
}

workflow diamond_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        dwi_channel = dwi_channel.groupTuple()
        dwi_image = dwi_channel.map{ [it[0], it[1]] }
        other_files = dwi_channel.map{ [it[0], it.subList(2, it.size()).inject([]){ c, t -> c + t }] }
        diamond(dwi_image.join(mask_channel).join(other_files), "reconstruct")
    emit:
        diamond = diamond.out.diamond
}