#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_root = false
params.include_sid = false
params.exclude_sid = false

include {
    prepare_metadata as pmeta_dwi;
    prepare_metadata as pmeta_rev;
    enforce_sid_convention as enforce_sid_convention_dwi;
    enforce_sid_convention as enforce_sid_convention_dwi_mask;
    enforce_sid_convention as enforce_sid_convention_anat;
    enforce_sid_convention as enforce_sid_convention_anat_mask;
    enforce_sid_convention as enforce_sid_convention_rev;
    enforce_sid_convention as enforce_sid_convention_rev_bval_bvec;
    enforce_sid_convention as enforce_sid_convention_pvf;
    enforce_sid_convention as enforce_sid_convention_metadata;
    enforce_sid_convention as enforce_sid_convention_rev_metadata
} from "../modules/processes/io.nf"

include {
    fill_missing_datapoints;
    exclude_missing_datapoints
} from "../modules/functions.nf"

def get_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}

workflow load_dataset {
    main:
        if ( !params.data_root )
            error "You must supply an input data root using --data_root"
        root = file(params.data_root)

        // Load all images
        dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { get_id(it.parent, root) }
        rev_channel = Channel.fromFilePairs("$root/**/*rev.nii.gz", size: 1, flat: true)
            { get_id(it.parent, root) }
        rev_bval_bvec_channel = Channel.fromFilePairs("$root/**/*rev.{bval,bvec}", size: 2, flat: true)
            { get_id(it.parent, root) }
        pvf_channel = Channel.fromFilePairs("$root/**/*{wm,gm,csf}_pvf.nii.gz", size: 3, flat: true)
            { get_id(it.parent, root) }
        anat_channel = Channel.fromFilePairs("$root/**/*t1.nii.gz", size: 1, flat: true)
            { get_id(it.parent, root) }
        dwi_mask_channel = Channel.fromFilePairs("$root/**/*dwi_mask.nii.gz", size: 1, flat: true)
            { get_id(it.parent, root) }
        anat_mask_channel = Channel.fromFilePairs("$root/**/*t1_mask.nii.gz", size: 1, flat: true)
            { get_id(it.parent, root) }
        dwi_json_channel = Channel.fromFilePairs("$root/**/*dwi.json", size: 1, flat: true)
            { get_id(it.parent, root) }
        rev_json_channel = Channel.fromFilePairs("$root/**/*rev.json", size: 1, flat: true)
            { get_id(it.parent, root) }

        // Enforce sid as suffix for all files
        enforce_sid_convention_anat(anat_channel.map{ it + ["t1"] })
        enforce_sid_convention_dwi(dwi_channel.map{ [it[0], it[1..-1], ["dwi"] * (it.size() - 1)] })
        enforce_sid_convention_rev(rev_channel.map{ it + ["rev"] })
        enforce_sid_convention_rev_bval_bvec(rev_bval_bvec_channel.map{ [it[0], it[1..-1], ["rev"] * (it.size() - 1)] })
        enforce_sid_convention_pvf(pvf_channel.map{ [it[0], it[1..-1], it[1..-1].collect{ i -> i.simpleName.tokenize("_")[-2] + "_pvf"} ] })
        enforce_sid_convention_dwi_mask(dwi_mask_channel.map{ it + ["dwi_mask"] })
        enforce_sid_convention_anat_mask(anat_mask_channel.map{ it + ["t1_mask"] })
        enforce_sid_convention_metadata(dwi_json_channel.map{ it + ["dwi"] })
        enforce_sid_convention_rev_metadata(rev_json_channel.map{ it + ["rev"] })

        // Filter subjects to get only those specified by the user
        ref_id_channel = enforce_sid_convention_anat.out.image
            .map{ [it[0]] }

        if (params.include_sid) {
            ref_id_channel = ref_id_channel.filter{ params.include_sid.contains(it[0]) }
        }
        if (params.exclude_sid) {
            ref_id_channel = ref_id_channel.filter{ !params.exclude_sid.contains(it[0]) }
        }

        // Unpack base images (T1 + DWI)
        anat_channel = ref_id_channel
            .join(enforce_sid_convention_anat.out.image)

        dwi_channel = ref_id_channel
            .join(enforce_sid_convention_dwi.out.image)
            .map{ [it[0], it[1][2], it[1][0], it[1][1]] }

        // Unpack reverse phase images
        rev_bval_bvec_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_rev_bval_bvec.out.image)
                .map{ it.flatten() },
            ref_id_channel,
            1, ["", ""]
        )

        rev_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_rev.out.image),
            ref_id_channel,
            1, [""]
        ).join(rev_bval_bvec_channel)

        // Unpack image metadata
        dwi_json_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_metadata.out.image),
            ref_id_channel,
            1, [""]
        )
        rev_json_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_rev_metadata.out.image),
            ref_id_channel,
            1, [""]
        )

        dwi_meta_channel = pmeta_dwi(
            dwi_channel
                .map{ it[0..1] }
                .join(dwi_json_channel)
                .map{ it + ["false"] }
        )
        rev_meta_channel = pmeta_rev(
            exclude_missing_datapoints(rev_channel, 1, "")
                .map{ it[0..1] }
                .join(rev_json_channel)
                .map{ it + ["true"] }
        )

        rev_meta_channel = fill_missing_datapoints(
            rev_meta_channel,
            ref_id_channel,
            1, [""]
        )

        // Unpack tissues PVF images
        pvf_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_pvf.out.image)
                .map{ [it[0], it[1].reverse()] },
            ref_id_channel,
            1, [[]]
        )

        // Load available masks (T1 and/or DWI)
        dwi_mask_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_dwi_mask.out.image),
            ref_id_channel,
            1, [""]
        )
        anat_mask_channel = fill_missing_datapoints(
            ref_id_channel
                .join(enforce_sid_convention_anat_mask.out.image),
            ref_id_channel,
            1, [""]
        )

    emit:
        dwi = dwi_channel
        dwi_mask = dwi_mask_channel
        anat = anat_channel
        anat_mask = anat_mask_channel
        rev = rev_channel 
        pvf = pvf_channel
        metadata = dwi_meta_channel
        rev_metadata = rev_meta_channel
}
