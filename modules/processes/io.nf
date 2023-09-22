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
    cache 'lenient'

    input:
        tuple val(sid), path(images), val(suffix)
    output:
        tuple val(sid), path("${sid}_*", includeInputs: true), emit: image
    script:
        def name = ""
        if ( (images instanceof Path ? images.getNameCount() : images.size()) == 1 ) {
            name = "${sid}_${suffix}.${extract_extension(images)}"
            if ( name != "${images.simpleName}.${extract_extension(images)}" ) {
                """
                ln -s $images $name
                """
            }
            else {
                """
                """
            }
        }
        else {
            def cmd = ""
            images.eachWithIndex{ img, i -> cmd += ( "${img.simpleName}.${extract_extension(img)}" == "${sid}_${suffix[i]}.${extract_extension(img)}" )
                ? ""
                : "ln -s $img ${sid}_${suffix[i]}.${extract_extension(img)}\n" 
            }
            """
            $cmd
            """
        }
}

process change_name {
    label "LIGHTSPEED"
    label "res_single_cpu"
    cache 'lenient'

    input:
        tuple val(sid), file(files)
        val(suffix)
    output:
        tuple val(sid), path("*__${suffix}*", includeInputs: true)
    script:
        def extension = ""
        def name = ""
        if ( (files instanceof Path ? files.getNameCount() : files.size()) == 1 ) {
            name = "${files.simpleName.split("__")[0]}__${suffix}.${extract_extension(files)}"
            if ( "${files.simpleName}.${extract_extension(files)}" != name ) {
                """
                ln -sf $files $name
                """
            }
            else {
                """
                """
            }
        }
        else {
            def cmd = ""
            for (f in files) {
                if ( !f.empty() ) {
                    name = "${f.simpleName.split("__")[0]}__${suffix}.${extract_extension(f)}"
                    if ( "${f.simpleName}.${extract_extension(f)}" != name )
                        cmd += "ln -sf $f $name\n"
                }
            }
            """
            $cmd
            """
        }
}

process rename_sequentially {
    label "LIGHTSPEED"
    label "res_single_cpu"
    cache 'lenient'

    input:
        tuple val(sid), path(files)
        val(suffix)
        val(start_character)
    output:
        tuple val(sid), path("*_${suffix}*", includeInputs: true)
    script:
        def commands = ""
        def name = ""
        if ( (files instanceof Path ? files.getNameCount() : files.size()) == 1 ) {
            name = "${sid}_${start_character}_${suffix}.${extract_extension(files)}"
            if ( name != "${files.simpleName}.${extract_extension(files)}" ) {
                commands += "ln -sf $files $name\n"
            }
        }
        else {
            for (f in files) {
                if ( !f.empty() ) {
                    name = "${sid}_${start_character}_${suffix}.${extract_extension(f)}"
                    if ( name != "${f.simpleName}.${extract_extension(f)}" ) {
                        commands += "ln -sf $f $name\n"
                    }
                    start_character = start_character.next()
                }
            }
        }
        """
        $commands
        """
}