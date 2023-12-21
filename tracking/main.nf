#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    rng_sampler
} from "../modules/functions.nf"

params.data_root = false

params.pft_tracking = true
params.local_tracking = true

params.pft_tracking_seeds = false
params.pft_n_tracking_seeds = 1
params.local_tracking_seeds = false
params.local_n_tracking_seeds = 1

params.pft_tracking_algorithm = ["prob", "det"]
params.local_tracking_algorithm = ["prob", "det", "eudx"]

params.local_tracking_mask = ["wm", "fa"]
params.ptf_seeding_mask = ["wm", "fa", "interface"]
params.local_seeding_mask = ["wm", "fa", "interface"]

params.pft_tracking_mask_fa_threshold = 0.1
params.local_tracking_mask_fa_threshold = 0.1

params.pft_random_seed = rng_sampler(params.pft_tracking_seeds, params.pft_n_tracking_seeds)
params.local_random_seed = rng_sampler(params.local_tracking_seeds, params.local_n_tracking_seeds)

include {
    enforce_sid_convention as enforce_sid_convention_fodf;
    enforce_sid_convention as enforce_sid_convention_dt_tensor;
    enforce_sid_convention as enforce_sid_convention_dt_fa;
    enforce_sid_convention as enforce_sid_convention_wm_mask;
    enforce_sid_convention as enforce_sid_convention_pvf;
    enforce_sid_convention as enforce_sid_convention_dwi
} from "../modules/processes/io.nf"

include {
    fill_missing_datapoints;
    exclude_missing_datapoints
} from "../modules/functions.nf"

include {
    ensemble_tracking_wkf
} from "../workflows/tracking.nf"


workflow {
    dataloader = load_tracking_data()
    ensemble_tracking_wkf(
        dataloader.dwi,
        dataloader.fodf,
        dataloader.tensor,
        dataloader.fa,
        dataloader.wm_mask,
        dataloader.pvf
    )
}


def get_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}


workflow load_tracking_data {
    main:
        if ( !params.data_root )
            error "You must supply an input data root using --data_root"
        root = file(params.data_root)

        dwi_channel = Channel.fromFilePairs("$root/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { get_id(it.parent, root) }

        fodf_channel = Channel.fromFilePairs("$root/**/fodf/*fodf.nii.gz", size: 1, flat: true)
            { get_id(it.parent.parent, root) }

        dt_tensor_channel = Channel.fromFilePairs("$root/**/dti/*dti_dti.nii.gz", size: 1, flat: true)
            { get_id(it.parent.parent, root) }

        dt_fa_channel = Channel.fromFilePairs("$root/**/dti/*dti_fa.nii.gz", size: 1, flat: true)
            { get_id(it.parent.parent, root) }

        wm_mask_channel = Channel.fromFilePairs("$root/**/segmentation/*_wm_mask.nii.gz", size: 1, flat: true)
            { get_id(it.parent.parent, root) }
        tissue_pvf_channel = Channel.fromFilePairs("$root/**/segmentation/*_3t_{csf,gm,wm}_pvf.nii.gz", size: 3, flat: true)
            { get_id(it.parent.parent, root) }

        enforce_sid_convention_dwi(dwi_channel.map{ [it[0], it[1..-1], ["dwi"] * (it.size() - 1)] })
        enforce_sid_convention_fodf(fodf_channel.map{ it + ["fodf"] })
        enforce_sid_convention_dt_tensor(dt_tensor_channel.map{ it + ["dti_dti"] })
        enforce_sid_convention_dt_fa(dt_fa_channel.map{ it + ["dti_fa"] })
        enforce_sid_convention_wm_mask(
            wm_mask_channel
                .filter{ !it[1].simpleName.contains("safe_wm_mask") }
                .map{ it + ["wm_mask"] }
        )
        enforce_sid_convention_pvf(
            tissue_pvf_channel.map{[
                it[0],
                it[1..-1],
                it[1..-1].collect{ i -> i.simpleName.tokenize("_")[-2] + "_pvf"}
            ]}
        )

        dwi_channel = enforce_sid_convention_dwi.out.image
            .map{ [it[0], it[1][2], it[1][0], it[1][1]] }

        ref_id_channel = enforce_sid_convention_fodf.out.image
            .map{ [it[0]] }

        dt_tensor_channel = fill_missing_datapoints(
            enforce_sid_convention_dt_tensor.out.image,
            ref_id_channel,
            1, [""]
        )
        dt_fa_channel = fill_missing_datapoints(
            enforce_sid_convention_dt_fa.out.image,
            ref_id_channel,
            1, [""]
        )
        wm_mask_channel = fill_missing_datapoints(
            enforce_sid_convention_wm_mask.out.image,
            ref_id_channel,
            1, [""]
        )
        tissue_pvf_channel = fill_missing_datapoints(
            enforce_sid_convention_pvf.out.image
                .map{ [it[0], it[1].reverse()] },
            ref_id_channel,
            1, [[]]
        )

    emit:
        dwi = dwi_channel
        fodf = fodf_channel
        tensor = dt_tensor_channel
        fa = dt_fa_channel
        wm_mask = wm_mask_channel
        pvf = tissue_pvf_channel
}