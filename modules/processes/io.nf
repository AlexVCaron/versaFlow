#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.default_readout = null
params.default_multiband_factor = null
params.default_is_interleaved = null
params.default_slicing_direction = null
params.default_phase_direction = null
params.default_acquisition_tensor_type = null

include { extract_extension } from '../functions.nf'

def metadata_from_params ( reverse ) {
    if ([
        params.default_readout,
        params.default_multiband_factor,
        params.default_is_interleaved,
        params.default_slicing_direction,
        params.default_phase_direction,
        params.default_acquisition_tensor_type
    ].any{ it == null })
        error "Some default acquisition parameters are not set, but are required. Set their values in the nextflow.config file."

    def direction = "${params.default_phase_direction}"
    if ( "$reverse" == "true" ) direction = direction.reverse()

    def margs = "--acq ${params.default_acquisition_tensor_type} --dir $direction --readout ${params.default_readout} --sd ${params.default_slicing_direction}"
    if ( params.default_is_interleaved )
        margs += " --interleaved"

    if ( params.default_multiband_factor && params.default_multiband_factor > 1 )
        margs += " --mb ${params.default_multiband_factor}"

    return margs
}

process prepare_metadata {
    label "LIGHTSPEED"
    label "res_single_cpu"
    input:
        tuple val(sid), path(image), file(metadata), val(reverse)
    output:
        tuple val(sid), path("${image.simpleName}_metadata.py")
    script:
        def args = ""
        if ( !metadata.empty() )
            args += "--json $metadata"
        else
            args = metadata_from_params(reverse)

        """
        mrhardi metadata --in $image $args
        """
}

process enforce_sid_convention {
    label "LIGHTSPEED"
    label "res_single_cpu"
    input:
        tuple val(sid), path(images), val(suffix)
    output:
        tuple val(sid), path("${task.workDir}/${sid}_*"), emit: image
    script:
        if ( (images instanceof Path ? images.getNameCount() : images.size()) == 1 ) {
            """
            ln -s $images ${task.workDir}/${sid}_${suffix}.${extract_extension(images)}
            """
        }
        else {
            def cmd = ""
            images.eachWithIndex{ img, i -> cmd += "ln -s $img ${task.workDir}/${sid}_${suffix[i]}.${extract_extension(img)}\n" }
            """
            $cmd
            """
        }
}