#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def metadata_from_params ( reverse ) {
    def direction = "${params.metadata.direction}"
    if ( "$reverse" == "false" )
        direction = direction.reverse()

    def margs = "--acq ${params.metadata.acquisition} --dir $direction --dwell ${params.metadata.dwell}"
    if ( params.metadata.multiband && params.metadata.multiband > 1 ) {
        margs += " --mb ${params.metadata.multiband} --sd ${params.metadata.slice_direction}"
    }
    args += margs

    return margs
}

process prepare_metadata {
    label "res_single_cpu"
    input:
        tuple val(sid), path(image), file(metadata), val(reverse)
    output:
        tuple val(sid), path("${image.simpleName}_metadata.py")
    script:
        args = ""
        if ( !metadata.empty() )
            args += "--json $metadata"
        else {
            args = metadata_from_params(reverse)
        }
        """
        magic-monkey metadata --in $image $args
        """
}
