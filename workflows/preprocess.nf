#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    merge_channels_non_blocking;
    replace_dwi_file;
    exclude_missing_datapoints;
    fill_missing_datapoints;
    filter_datapoints;
    get_config_path;
    collect_paths;
    is_path_list
} from '../modules/functions.nf'
include {
    extract_b0 as dwi_b0;
    extract_b0 as extract_topup_b0;
    extract_b0 as extract_b0_preprocessed
} from '../modules/processes/preprocess.nf'
include {
    scil_compute_dti_fa
} from '../modules/processes/measure.nf'
include {
    ants_transform as ants_transform_base_t1;
    ants_transform as ants_transform_base_dwi;
    ants_transform as ants_transform_syn_t1;
    ants_transform as ants_transform_syn_dwi;
    ants_transform as ants_transform_base_wm;
    ants_transform as ants_transform_base_gm;
    ants_transform as ants_transform_base_csf;
    ants_transform as ants_transform_syn_wm;
    ants_transform as ants_transform_syn_gm;
    ants_transform as ants_transform_syn_csf;
    ants_transform as ants_transform_base_raw_t1;
    ants_transform as ants_transform_syn_raw_t1;
    ants_transform as ants_transform_wm_mask;
    ants_transform as ants_transform_gm_mask;
    ants_transform as ants_transform_csf_mask;
    ants_transform as ants_transform_safe_wm_mask;
    ants_transform as ants_transform_raw_t1_mask
} from '../modules/processes/register.nf'
include {
    convert_float_to_integer as convert_wm_segmentation;
    convert_float_to_integer as convert_gm_segmentation;
    convert_float_to_integer as convert_csf_segmentation;
    convert_float_to_integer as dwi_mask_convert_datatype;
    convert_float_to_integer as t1_mask_convert_datatype;
    crop_image as crop_dwi;
    crop_image as crop_t1;
    crop_image as crop_wm;
    crop_image as crop_gm;
    crop_image as crop_csf;
    crop_image as crop_wm_mask;
    crop_image as crop_gm_mask;
    crop_image as crop_csf_mask;
    crop_image as crop_safe_wm_mask;
    crop_image as crop_raw_dwi;
    crop_image as crop_raw_t1;
    bet_mask;
    fit_bounding_box;
    merge_masks;
    check_odd_dimensions;
    pvf_to_mask;
    validate_gradients
} from '../modules/processes/utils.nf'
include {
    gibbs_removal as dwi_gibbs_removal;
    gibbs_removal as rev_gibbs_removal;
    nlmeans_denoise;
    ants_gaussian_denoise;
    normalize_inter_b0
} from '../modules/processes/denoise.nf'
include {
    scilpy_resample_to_reference as resample_wm;
    scilpy_resample_to_reference as resample_gm;
    scilpy_resample_to_reference as resample_csf;
    scilpy_resample_to_reference as resample_t1;
    scilpy_resample_to_reference as resample_dwi;
    scilpy_resample_to_reference as resample_raw_dwi;
    scilpy_resample_to_reference as resample_raw_t1;
    resampling_reference
} from '../modules/processes/upsample.nf'
include {
    registration_wkf as dwi_mask_registration_wkf;
    dwi_denoise_wkf;
    dwi_denoise_wkf as rev_denoise_wkf;
    n4_denoise_wkf;
    n4_denoise_wkf as n4_denoise_t1_to_b0_wkf;
    squash_wkf;
    squash_wkf as squash_raw_wkf;
    squash_wkf as squash_for_topup_wkf;
    topup_wkf;
    apply_topup_wkf;
    apply_topup_wkf as raw_apply_topup_wkf;
    eddy_wkf
} from "../modules/workflows/preprocess.nf"
include {
    t1_mask_to_b0;
    t12b0_registration as t1_registration_wkf
} from '../modules/workflows/t1_registration.nf'
include {
    segment_nmt_wkf;
    segment_wm_wkf
} from '../modules/workflows/segment.nf'
include {
    change_name as rename_dwi_for_topup;
    change_name as rename_rev_for_topup;
    change_name as rename_dwi_meta_for_topup;
    change_name as rename_rev_meta_for_topup;
    change_name as rename_topup_corrected_dwi;
    change_name as rename_topup_corrected_metadata;
    change_name as rename_transformed_raw_dwi;
    change_name as rename_transformed_raw_rev;
    change_name as rename_transformed_raw_metadata;
    change_name as rename_transformed_raw_rev_metadata;
    change_name as rename_processed_dwi;
    change_name as rename_processed_dwi_metadata;
    change_name as rename_processed_dwi_mask;
    change_name as rename_dwi_mask;
    change_name as rename_dwi_for_eddy;
    change_name as rename_rev_for_eddy;
    change_name as rename_dwi_metadata_for_eddy;
    change_name as rename_rev_metadata_for_eddy
} from '../modules/processes/io.nf'

// Preprocess workflow parameters
params.gaussian_noise_correction = true
params.gibbs_ringing_correction = true
params.dwi_mask_from_t1_mask = true
params.topup_correction = true
params.eddy_correction = true
params.dwi_intensity_normalization = true
params.resample_data = true
params.register_t1_to_dwi = true
params.generate_tissue_segmentation = false
params.generate_wm_segmentation = true
params.raw_to_processed_space = false

// T1 preprocess workflow parameters
params.denoise_t1 = true
params.nlmeans_t1 = true
params.t1_intensity_normalization = true

params.quick_t1_mask_registration = true
params.quick_denoised_t1_registration = false
params.t1_registration_in_subject_space = false

params.b02t1_mask_registration_config = file("${get_config_path()}/b02t1_mask_registration_config.py")
params.t1_mask_to_topup_b0_registration_config = file("${get_config_path()}/t1_mask_to_topup_b0_registration_config.py")
params.ants_transform_base_config = file("${get_config_path()}/ants_transform_base_config.py")
params.ants_transform_mask_config = file("${get_config_path()}/ants_transform_mask_config.py")
params.extract_mean_b0_base_config = file("${get_config_path()}/extract_mean_b0_base_config.py")
params.dwi_n4_normalization_config = file("${get_config_path()}/dwi_n4_normalization_config.py")
params.dwi_n4_normalization_quick_config = file("${get_config_path()}/dwi_n4_normalization_quick_config.py")
params.t1_n4_normalization_config = file("${get_config_path()}/t1_n4_normalization_config.py")
params.b0_to_b0_normalization_config = file("${get_config_path()}/b0_to_b0_normalization_config.py")


workflow preprocess_wkf {
    take:
        dwi_channel
        rev_channel
        t1_channel
        pvf_channel
        meta_channel
        rev_meta_channel
        dwi_mask_channel
        t1_mask_channel
    main:
        // Keep sid references for channel management
        ref_id_channel = dwi_channel.map{ [it[0]] }
        absent_dwi_mask_id_channel = filter_datapoints(
            dwi_mask_channel,
            { it[1] == "" }
        ).map{ [it[0]] }

        // T1 preprocessing
        t1_preprocess_wkf(t1_channel, t1_mask_channel)
        t1_channel = t1_preprocess_wkf.out.t1
        t1_mask_channel = t1_preprocess_wkf.out.mask

        // Fix odd number of slices in phase direction for Topup
        check_odd_dimensions(
            dwi_channel
                .join(rev_channel)
                .join(dwi_mask_channel)
                .join(collect_paths(meta_channel.join(rev_meta_channel))),
            "preprocess"
        )

        rev_bval_bvec_channel = fill_missing_datapoints(
            check_odd_dimensions.out.rev_bval_bvec,
            ref_id_channel,
            1, ["", ""]
        )

        dwi_channel = check_odd_dimensions.out.dwi
        rev_channel = fill_missing_datapoints(
            check_odd_dimensions.out.rev,
            ref_id_channel,
            1, [""]
        ).join(rev_bval_bvec_channel)
        dwi_mask_channel = fill_missing_datapoints(
            check_odd_dimensions.out.mask,
            ref_id_channel,
            1, [""]
        )

        meta_channel = check_odd_dimensions.out.metadata
            .map{ [it[0], it[1] instanceof Path ? it[1] : it[1].findAll{ i -> !i.simpleName.contains("_rev") }].flatten() }
        rev_meta_channel = check_odd_dimensions.out.metadata
            .map{ [it[0], it[1].findAll{ i -> i.simpleName.contains("_rev") }].flatten() }

        // Copy input channels for later
        dwi_channel.tap{ raw_dwi_channel }
        rev_channel.tap{ raw_rev_channel }
        t1_channel.tap{ raw_t1_channel }
        t1_mask_channel.tap{ raw_t1_mask_channel }
        meta_channel.tap{ raw_meta_channel }
        rev_meta_channel.tap{ raw_rev_meta_channel }

        // Perform DWI and b0 denoising
        if ( params.gaussian_noise_correction ) {
            dwi_denoise_wkf(dwi_channel, dwi_mask_channel, meta_channel, "true")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_denoise_wkf.out.image)
            meta_channel = dwi_denoise_wkf.out.metadata

            rev_denoise_wkf(rev_channel, dwi_mask_channel, rev_meta_channel, "false")
            rev_channel = replace_dwi_file(rev_channel, rev_denoise_wkf.out.image)
            rev_meta_channel = rev_denoise_wkf.out.metadata
        }

        // Perform DWI and b0 gibbs correction
        if ( params.gibbs_ringing_correction ) {
            dwi_gibbs_removal(dwi_channel.map{ it[0..1] }.join(meta_channel), "preprocess", "true")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_gibbs_removal.out.image)
            meta_channel = dwi_gibbs_removal.out.metadata

            rev_gibbs_removal(
                exclude_missing_datapoints(rev_channel.map{ it[0..1] }.join(rev_meta_channel), 1, ""),
                "preprocess", "false"
            )
            rev_channel = replace_dwi_file(
                rev_channel,
                fill_missing_datapoints(rev_gibbs_removal.out.image, ref_id_channel, 1, [""])
            )
            rev_meta_channel = fill_missing_datapoints(rev_gibbs_removal.out.metadata, ref_id_channel, 1, [""])
        }

        // Perform DWI signal normalization between b0 volumes
        if ( params.normalize_inter_b0 ) {
            normalize_inter_b0(
                dwi_channel
                    .map{ it[0..2] }
                    .join(rev_channel.map{ it[0..2] })
                    .join(meta_channel)
                    .join(rev_meta_channel),
                "preprocess",
                params.b0_to_b0_normalization_config
            )
            dwi_channel = replace_dwi_file(dwi_channel, normalize_inter_b0.out.dwi)
            meta_channel = normalize_inter_b0.out.dwi_metadata
            rev_channel = replace_dwi_file(
                rev_channel,
                fill_missing_datapoints(normalize_inter_b0.out.rev, ref_id_channel, 1, [""])
            )
            rev_meta_channel = fill_missing_datapoints(
                normalize_inter_b0.out.rev_metadata,
                ref_id_channel,
                1, [""]
            )
        }

        // Average consecutive b0 volumes just like for Topup
        squash_wkf(
            dwi_channel,
            rev_channel,
            meta_channel.join(rev_meta_channel),
            ""
        )

        dwi_channel = squash_wkf.out.dwi
        rev_channel = squash_wkf.out.rev

        meta_channel = squash_wkf.out.metadata
            .map{ it.flatten() }
            .map{ [it[0], it[1..-1]] }
            .map{ [it[0], it[1].findAll{ i -> !i.simpleName.contains("_rev") }].flatten() }
        rev_meta_channel = squash_wkf.out.metadata
            .map{ it.flatten() }
            .map{ [it[0], it[1..-1]] }
            .map{ [it[0], it[1].findAll{ i -> i.simpleName.contains("_rev") }].flatten() }

        // Extract mean b0
        dwi_b0(
            dwi_channel
                .map{ it[0..2] }
                .join(meta_channel.map{ [it[0], it[1..-1]] }),
            "preprocess",
            "false",
            params.extract_mean_b0_base_config
        )

        b0_channel = dwi_b0.out.b0
        b0_metadata_channel = dwi_b0.out.metadata

        // Topup correction
        topup2eddy_channel = ref_id_channel.map{ it + ["", "", []] }
        dwi_after_topup_channel = dwi_channel
        meta_after_topup_channel = meta_channel
        if ( params.topup_correction ) {
            ref_rev_id_channel = exclude_missing_datapoints(
                raw_rev_channel, 1, ""
            ).map{ [it[0]] }
            excluded_id_channel = filter_datapoints(
                raw_rev_channel, { it[1] == "" }
            ).map{ [it[0]] }

            // Average consecutive b0 volumes to speed up Topup
            squash_for_topup_wkf(
                ref_rev_id_channel.join(raw_dwi_channel),
                ref_rev_id_channel.join(raw_rev_channel),
                ref_rev_id_channel.join(raw_meta_channel).join(raw_rev_meta_channel),
                ""
            )

            // Run Topup sub-workflow
            topup_dwi_channel = squash_for_topup_wkf.out.dwi
            topup_rev_channel = squash_for_topup_wkf.out.rev
            topup_wkf(
                topup_dwi_channel,
                topup_rev_channel,
                squash_for_topup_wkf.out.metadata.map{ it.flatten() }
            )

            topup2eddy_channel = topup_wkf.out.param
                .join(topup_wkf.out.prefix)
                .join(topup_wkf.out.topup.map{ [it[0], it.subList(1, it.size())] })

            dwi2topup_channel = rename_dwi_for_topup(
                collect_paths(topup_wkf.out.topupable_indexes.join(dwi_channel)),
                "dwi_to_topup"
            ).map{ [it[0], it[1][2], it[1][0], it[1][1]] }

            rev2topup_channel = rename_rev_for_topup(
                collect_paths(topup_wkf.out.topupable_indexes.join(rev_channel)),
                "dwi_to_topup_rev"
            ).map{ it.flatten() }.map{ it.size() == 4 ? [it[0], it[3], it[1], it[2]] : it + ["", ""] }

            meta2topup_channel = rename_dwi_meta_for_topup(
                topup_wkf.out.in_metadata_w_topup
                    .map{ [it[0], it[1][(0..<it[1].size()).step(2)]] },
                "dwi_to_topup_metadata"
            ).map{ it.flatten() }
            rev_meta2topup_channel = rename_rev_meta_for_topup(
                topup_wkf.out.in_metadata_w_topup
                    .map{ [it[0], it[1][(1..<it[1].size()).step(2)]] },
                "dwi_to_topup_rev_metadata"
            ).map{ it.flatten() }

            // Applied estimated susceptibility correction to DWI
            apply_topup_wkf(
                dwi2topup_channel,
                rev2topup_channel,
                topup2eddy_channel,
                meta2topup_channel
                    .join(rev_meta2topup_channel)
                    .map{ [it[0], it[1..-1]] },
                ""
            )

            topup_corrected_dwi_channel = rename_topup_corrected_dwi(
                collect_paths(apply_topup_wkf.out.dwi),
                "topup_corrected"
            ).map{ [it[0], it[1][2], it[1][0], it[1][1]] }
            topup_corrected_dwi_meta_channel = rename_topup_corrected_metadata(
                collect_paths(apply_topup_wkf.out.metadata),
                "topup_corrected_metadata"
            ).map{ it.flatten() }

            dwi_after_topup_channel = excluded_id_channel
                .join(dwi_channel)
                .mix(topup_corrected_dwi_channel)
            meta_after_topup_channel = excluded_id_channel
                .join(meta_channel)
                .mix(topup_corrected_dwi_meta_channel)

            // Get average susceptibility corrected b0
            extract_topup_b0(
                topup_corrected_dwi_channel
                    .map{ it[0..2] }
                    .join(collect_paths(topup_corrected_dwi_meta_channel)),
                "preprocess",
                "true",
                params.extract_mean_b0_base_config
            )

            b0_channel = excluded_id_channel
                .join(b0_channel)
                .mix(extract_topup_b0.out.b0)
            b0_metadata_channel = excluded_id_channel
                .join(b0_metadata_channel)
                .mix(extract_topup_b0.out.metadata)

            if ( !params.eddy_correction ) {
                dwi_channel = dwi_after_topup_channel
                meta_channel = meta_after_topup_channel
            }
            else {
                topup2eddy_channel = fill_missing_datapoints(
                    topup2eddy_channel,
                    ref_id_channel,
                    1, ["", "", []]
                )
            }

            // Apply susceptibility corrections to raw images (for comparison)
            if ( params.raw_to_processed_space ) {
                raw_dwi_channel = rename_transformed_raw_dwi(
                    collect_paths(raw_dwi_channel),
                    "raw"
                ).map{ [it[0], it[1][2], it[1][0], it[1][1]] }
                raw_rev_channel = rename_transformed_raw_rev(
                    collect_paths(raw_rev_channel).filter{ it[1] },
                    "raw"
                )
                raw_rev_channel = raw_rev_channel
                    .map{ it.flatten() }
                    .map{ it.size() == 4 ? [it[0], it[3], it[1], it[2]] : it + ["", ""] }

                raw_meta_channel = excluded_id_channel.join(raw_meta_channel)
                    .mix(
                        rename_transformed_raw_metadata(
                            collect_paths(meta2topup_channel).filter{ it[1] },
                            "raw_metadata"
                        ).map{ it.flatten() }
                    )

                raw_rev_meta_channel = excluded_id_channel.join(raw_meta_channel)
                    .mix(
                        rename_transformed_raw_rev_metadata(
                            collect_paths(rev_meta2topup_channel),
                            "raw_metadata"
                        ).map{ it.flatten() }
                    )

                raw_apply_topup_wkf(
                    topup_wkf.out.topupable_indexes.join(raw_dwi_channel),
                    topup_wkf.out.topupable_indexes
                        .join(raw_rev_channel)
                        .map{ it[0..1] },
                    topup2eddy_channel,
                    collect_paths(raw_meta_channel.join(raw_rev_meta_channel)),
                    "raw"
                )

                raw_dwi_channel = excluded_id_channel
                    .join(raw_dwi_channel)
                    .mix(raw_apply_topup_wkf.out.dwi)

                raw_meta_channel = excluded_id_channel
                    .join(raw_meta_channel)
                    .mix(raw_apply_topup_wkf.out.metadata)
            }
        }

        empty_dwi_mask_id_channel = filter_datapoints(
            dwi_mask_channel,
            { it[1] == "" }
        ).map{ [it[0]] }
        dwi_mask_channel = exclude_missing_datapoints(dwi_mask_channel, 1, "")

        // Compute brain mask for the DWI (when missing)
        bet_mask(
            empty_dwi_mask_id_channel.join(b0_channel),
            "preprocess",
            "${!params.dwi_mask_from_t1_mask}",
            "dwi_mask"
        )
        dwi_mask_channel = dwi_mask_channel.mix(bet_mask.out.mask)

        // Get better mask for the DWI from the T1 (when missing and if present)
        if ( params.dwi_mask_from_t1_mask ) {
            existing_t1_mask_id_channel = exclude_missing_datapoints(
                t1_mask_channel,
                1, ""
            ).map{ [it[0]] }
            absent_t1_mask_id_channel = filter_datapoints(
                t1_mask_channel,
                { it[1] == "" }
            ).map{ [it[0]] }

            n4_denoise_t1_to_b0_wkf(
                existing_t1_mask_id_channel.join(dwi_after_topup_channel.map{ it[0..1] }),
                existing_t1_mask_id_channel.join(b0_channel),
                existing_t1_mask_id_channel.join(dwi_mask_channel),
                existing_t1_mask_id_channel.join(meta_after_topup_channel),
                params.dwi_n4_normalization_quick_config,
                false
            )

            t1_mask_to_b0(
                replace_dwi_file(dwi_after_topup_channel, n4_denoise_t1_to_b0_wkf.out.image),
                existing_t1_mask_id_channel.join(t1_channel),
                existing_t1_mask_id_channel.join(t1_mask_channel),
                "false"
            )

            t1_mask_convert_datatype(
                t1_mask_to_b0.out.mask,
                "uint8", "preprocess",
                !params.register_t1_to_dwi,
                "dwi_mask", ""
            )

            dwi_mask_channel = t1_mask_convert_datatype.out.image
                .mix(absent_t1_mask_id_channel.join(dwi_mask_channel))
        }

        // Perform Eddy currents and motion correction on the DWI
        if ( params.eddy_correction ) {
            dwi_channel = rename_dwi_for_eddy(
                collect_paths(dwi_channel),
                "to_eddy"
            ).map{ [it[0], it[1][2], it[1][0], it[1][1]] }

            rev_channel = rename_rev_for_eddy(
                collect_paths(rev_channel).filter{ it[1] },
                "to_eddy"
            )
            rev_channel = rev_channel
                .map{ it.flatten() }
                .map{ it.size() == 4 ? [it[0], it[3], it[1], it[2]] : it + ["", ""] }

            rev_channel = fill_missing_datapoints(
                rev_channel,
                ref_id_channel,
                1, ["", "", ""]
            )

            meta_channel = rename_dwi_metadata_for_eddy(
                collect_paths(meta_channel).filter{ it[1] },
                "to_eddy_metadata"
            ).map{ it.flatten() }

            rev_meta_channel = rename_rev_metadata_for_eddy(
                collect_paths(rev_meta_channel).filter{ it[1] },
                "to_eddy_metadata"
            ).map{ it.flatten() }

            rev_meta_channel = fill_missing_datapoints(
                rev_meta_channel,
                ref_id_channel,
                1, [""]
            )

            // Run Eddy sub-workflow
            eddy_wkf(
                dwi_channel,
                dwi_mask_channel,
                topup2eddy_channel,
                b0_channel,
                rev_channel,
                meta_channel.join(rev_meta_channel)
            )

            dwi_channel = eddy_wkf.out.dwi
                .join(eddy_wkf.out.bval)
                .join(eddy_wkf.out.bvec)
            meta_channel = eddy_wkf.out.metadata
        }

        // Perform intensity normalization on the DWI
        if ( params.dwi_intensity_normalization ) {
            // Run N4 sub-workflow
            n4_denoise_wkf(
                dwi_channel
                    .map{ it.subList(0, 2) },
                b0_channel,
                dwi_mask_channel,
                meta_channel,
                params.dwi_n4_normalization_config,
                true
            )

            dwi_channel = replace_dwi_file(dwi_channel, n4_denoise_wkf.out.image)
            meta_channel = n4_denoise_wkf.out.metadata
        }

        absent_t1_mask_id_channel = filter_datapoints(
            t1_mask_channel,
            { it[1] == "" }
        ).map{ [it[0]] }
        existing_t1_mask_id_channel = exclude_missing_datapoints(
            t1_mask_channel,
            1, ""
        ).map{ [it[0]] }

        // Get DWI mask in T1 space (for missing T1 masks)
        dwi_mask_registration_wkf(
            absent_t1_mask_id_channel.join(t1_channel).map{ [it[0], [it[1]]] },
            absent_t1_mask_id_channel.join(b0_channel).map{ [it[0], [it[1]]] },
            absent_t1_mask_id_channel.join(dwi_mask_channel),
            null,
            null,
            absent_t1_mask_id_channel
                .join(b0_metadata_channel)
                .map{ it[0..1] + [""] },
            "",
            false, "", "",
            params.b02t1_mask_registration_config,
            null
        )

        dwi_mask_convert_datatype(
            dwi_mask_registration_wkf.out.image,
            "uint8", "preprocess",
            true,
            "t1_mask", ""
        )

        t1_mask_channel = existing_t1_mask_id_channel
            .join(t1_mask_channel)
            .mix(dwi_mask_convert_datatype.out.image)
        raw_t1_mask_channel = t1_mask_channel
        raw_dwi_mask_channel = dwi_mask_channel

        // Compute best resampling reference
        resampling_reference(
            collect_paths(dwi_channel.map{ it[0..1] }.join(t1_channel)),
            "preprocess"
        )

        reference_channel = resampling_reference.out.reference
        pvf_to_resample_channel = pvf_channel.filter{ !it[1].isEmpty() }

        // Resample all volumes
        resample_dwi(
            dwi_channel
                .map{ it[0..1] }
                .join(reference_channel)
                .join(dwi_mask_channel)
                .join(meta_channel),
            "preprocess", "lin",
            true, true,
            "dwi_mask", ""
        )

        resample_t1(
            t1_channel
                .join(reference_channel)
                .join(t1_mask_channel)
                .map{ it + [""] },
            "preprocess", "lin",
            true, true,
            "t1_mask", ""
        )

        resample_wm(
            pvf_to_resample_channel
                .map{ [it[0], it[1][0]] }
                .join(reference_channel)
                .join(t1_mask_channel)
                .map{ it + [""] },
            "preprocess", "nn",
            true, false,
            "", "segmentation"
        )
        resample_gm(
            pvf_to_resample_channel
                .map{ [it[0], it[1][1]] }
                .join(reference_channel)
                .join(t1_mask_channel)
                .map{ it + [""] },
            "preprocess", "nn",
            true, false,
            "", "segmentation"
        )
        resample_csf(
            pvf_to_resample_channel
                .map{ [it[0], it[1][2]] }
                .join(reference_channel)
                .join(t1_mask_channel)
                .map{ it + [""] },
            "preprocess", "nn",
            true, false,
            "", "segmentation"
        )

        dwi_channel = replace_dwi_file(dwi_channel, resample_dwi.out.image)
        dwi_mask_channel = resample_dwi.out.mask
        meta_channel = resample_dwi.out.metadata
        t1_channel = resample_t1.out.image
        t1_mask_channel = resample_t1.out.mask

        pvf_channel = resample_wm.out.image
            .join(resample_gm.out.image)
            .join(resample_csf.out.image)
            .map{ [it[0], it[1..-1]] }
            .mix(pvf_channel.filter{ it[1].isEmpty() })

        if ( params.raw_to_processed_space ) {
            resample_raw_dwi(
                raw_dwi_channel
                    .map{ it[0..1] }
                    .join(reference_channel)
                    .join(raw_dwi_mask_channel)
                    .join(raw_meta_channel),
                "preprocess", "lin",
                true, true,
                "dwi_mask", "raw"
            )

            resample_raw_t1(
                raw_t1_channel
                    .join(reference_channel)
                    .join(raw_t1_mask_channel)
                    .map{ it + [""] },
                "preprocess", "lin",
                true, true,
                "t1_mask", "raw"
            )

            raw_dwi_channel = replace_dwi_file(raw_dwi_channel, resample_raw_dwi.out.image)
            raw_meta_channel = resample_raw_dwi.out.metadata
            raw_t1_channel = resample_raw_t1.out.image
            raw_t1_mask_channel = resample_raw_t1.out.mask
            raw_dwi_mask_channel = resample_raw_dwi.out.mask
        }

        // Get tissue masks from PVF
        pvf_to_mask(
             pvf_channel
                .filter{ !it[1].isEmpty() }
                .join(dwi_mask_channel),
            "preprocess",
            "segmentation"
        )

        tissue_mask_channel = collect_paths(
            pvf_to_mask.out.wm_mask
                .join(pvf_to_mask.out.gm_mask)
                .join(pvf_to_mask.out.csf_mask)
        ).mix(pvf_channel.filter{ it[1].isEmpty() })

        safe_wm_mask_channel = pvf_channel
            .filter{ it[1].isEmpty() }
            .map{ [it[0], ""] }
            .mix(pvf_to_mask.out.safe_wm_mask)

        // Generate tissue segmentation from T1
        if ( params.generate_tissue_segmentation ) {
            absent_pvf_id_channel = pvf_channel
                .filter{ it[1].isEmpty() }
                .map{ [it[0]] }

            segment_nmt_wkf(
                absent_pvf_id_channel.join(t1_channel),
                absent_pvf_id_channel.join(t1_mask_channel)
            )

            pvf_channel = segment_nmt_wkf.out.volume_fractions
                .map{ [it[0], it[1].reverse()] }
                .mix(pvf_channel.filter{ !it[1].isEmpty() })

            tissue_mask_channel = segment_nmt_wkf.out.tissue_masks
                .mix(tissue_mask_channel.filter{ it[1] })
            safe_wm_mask_channel = segment_nmt_wkf.out.safe_wm_mask
                .mix(safe_wm_mask_channel.filter{ it[1] })
        }

        // Register T1 to diffusion space (DWI) with masks and segmentations
        if ( params.register_t1_to_dwi ) {
            t1_registration_wkf(
                dwi_channel,
                t1_channel,
                t1_mask_channel,
                dwi_mask_channel,
                meta_channel,
                true,
                true,
                params.quick_denoised_t1_registration,
                params.t1_registration_in_subject_space
            )

            t1_channel = t1_registration_wkf.out.t1
            t1_mask_channel = t1_registration_wkf.out.mask
            dwi_mask_channel = rename_dwi_mask(
                collect_paths(t1_registration_wkf.out.mask),
                "dwi_mask"
            ).map{ it.flatten() }

            pvf_to_register_channel = pvf_channel
                .filter{ !it[1].isEmpty() }
            tissue_mask_to_register_channel = tissue_mask_channel
                .filter{ !it[1].isEmpty() }

            ants_transform_base_wm(
                pvf_to_register_channel.map{ [it[0], it[1][0]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_base_config
            )
            ants_transform_base_gm(
                pvf_to_register_channel
                    .map{ [it[0], it[1][1]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_base_config
            )
            ants_transform_base_csf(
                pvf_to_register_channel
                    .map{ [it[0], it[1][2]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_base_config
            )
            ants_transform_wm_mask(
                tissue_mask_to_register_channel
                    .map{ [it[0], it[1][1]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_mask_config
            )
            ants_transform_gm_mask(
                tissue_mask_to_register_channel
                    .map{ [it[0], it[1][2]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_mask_config
            )
            ants_transform_csf_mask(
                tissue_mask_to_register_channel
                    .map{ [it[0], it[1][2]] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_mask_config
            )
            ants_transform_safe_wm_mask(
                safe_wm_mask_channel
                    .filter{ it[1] }
                    .join(t1_channel)
                    .join(t1_registration_wkf.out.transform)
                    .map{ it + ["", ""] },
                "preprocess", "segmentation", "true", "",
                params.ants_transform_mask_config
            )

            pvf_channel = collect_paths(
                ants_transform_base_wm.out.image
                    .join(ants_transform_base_gm.out.image)
                    .join(ants_transform_base_csf.out.image)
            ).mix(pvf_channel.filter{ it[1].isEmpty() })

            tissue_mask_channel = collect_paths(
                ants_transform_wm_mask.out.image
                    .join(ants_transform_gm_mask.out.image)
                    .join(ants_transform_csf_mask.out.image)
            ).mix(pvf_channel.filter{ it[1].isEmpty() })

            safe_wm_mask_channel = pvf_channel
                .filter{ it[1].isEmpty() }
                .map{ [it[0], ""] }
                .mix(ants_transform_safe_wm_mask.out.image)

            if ( params.raw_to_processed_space ) {
                ants_transform_base_raw_t1(
                    raw_t1_channel
                        .join(t1_registration_wkf.out.reference)
                        .join(t1_registration_wkf.out.transform)
                        .map{ it + ["", ""] },
                    "preprocess",
                    "raw",
                    "true",
                    "t1",
                    params.ants_transform_base_config
                )

                ants_transform_raw_t1_mask(
                    raw_t1_mask_channel
                        .join(t1_registration_wkf.out.reference)
                        .join(t1_registration_wkf.out.transform)
                        .map{ it + ["", ""] },
                    "preprocess",
                    "raw",
                    "true",
                    "t1_mask",
                    params.ants_transform_mask_config
                )

                raw_t1_channel = ants_transform_base_raw_t1.out.image
                raw_t1_mask_channel = ants_transform_raw_t1_mask.out.image
            }
        }

        crop_dwi(
            dwi_channel
                .map{ it[0..1] }
                .join(dwi_mask_channel)
                .map{ it + [""] }
                .join(collect_paths(meta_channel)),
            "preprocess", true, "dwi_mask", ""
        )

        fit_bounding_box(
            t1_channel
                .join(dwi_channel.map{ it[0..1] })
                .join(crop_dwi.out.bbox),
            "preprocess"
        )

        dwi_bbox_channel = crop_dwi.out.bbox
        t1_bbox_channel = fit_bounding_box.out.bbox
        pvf_to_crop_channel = pvf_channel
            .filter{ !it[1].isEmpty() }
        tissue_mask_to_crop_channel = tissue_mask_channel
            .filter{ !it[1].isEmpty() }

        crop_t1(
            t1_channel
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", true, "t1_mask", ""
        )

        crop_wm(
            pvf_to_crop_channel
                .map{ [it[0], it[1][0]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )
        crop_gm(
            pvf_to_crop_channel
                .map{ [it[0], it[1][1]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )
        crop_csf(
            pvf_to_crop_channel
                .map{ [it[0], it[1][2]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )

        crop_wm_mask(
            tissue_mask_to_crop_channel
                .map{ [it[0], it[1][0]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )
        crop_gm_mask(
            tissue_mask_to_crop_channel
                .map{ [it[0], it[1][1]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )
        crop_csf_mask(
            tissue_mask_to_crop_channel
                .map{ [it[0], it[1][2]] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )
        crop_safe_wm_mask(
            safe_wm_mask_channel
                .filter{ it[1] }
                .join(t1_mask_channel)
                .join(t1_bbox_channel)
                .map{ it + [""] },
            "preprocess", false, "", "segmentation"
        )

        dwi_channel = replace_dwi_file(dwi_channel, crop_dwi.out.image)
        dwi_mask_channel = crop_dwi.out.mask
        t1_channel = crop_t1.out.image
        t1_mask_channel = crop_t1.out.mask

        pvf_channel = collect_paths(
            crop_wm.out.image
                .join(crop_gm.out.image)
                .join(crop_csf.out.image)
        ).mix(pvf_channel.filter{ it[1].isEmpty() })

        tissue_mask_channel = collect_paths(
            crop_wm_mask.out.image
                .join(crop_gm_mask.out.image)
                .join(crop_csf_mask.out.image)
        ).mix(pvf_channel.filter{ it[1].isEmpty() })

        safe_wm_mask_channel = pvf_channel
            .filter{ it[1].isEmpty() }
            .map{ [it[0], ""] }
            .mix(crop_safe_wm_mask.out.image)

        if ( params.raw_to_processed_space ) {
            crop_raw_dwi(
                raw_dwi_channel
                    .map{ it[0..1] }
                    .join(raw_dwi_mask_channel)
                    .join(dwi_bbox_channel)
                    .join(collect_paths(raw_meta_channel)),
                "preprocess", true, "dwi_mask", "raw"
            )
            crop_raw_t1(
                raw_t1_channel
                    .join(raw_t1_mask_channel)
                    .join(t1_bbox_channel)
                    .map{ it + [""] },
                "preprocess", true, "t1_mask", "raw"
            )

            raw_dwi_channel = replace_dwi_file(raw_dwi_channel, crop_raw_dwi.out.image)
            raw_dwi_mask_channel = crop_raw_dwi.out.mask
            raw_meta_channel = crop_raw_dwi.out.metadata
            raw_t1_channel = crop_raw_t1.out.image
            raw_t1_mask_channel = crop_raw_dwi.out.mask
        }

        extract_b0_preprocessed(
            dwi_channel
                .map{ it[0..2] }
                .join(meta_channel),
            "preprocess", "true",
            params.extract_mean_b0_base_config
        )

        validate_gradients_wkf(dwi_channel, dwi_mask_channel)

        dwi_channel = rename_processed_dwi(
            collect_paths(validate_gradients_wkf.out.dwi),
            "dwi_preprocessed"
        ).map{ [it[0], it[1][2], it[1][0], it[1][1]] }
        meta_channel = rename_processed_dwi_metadata(
            collect_paths(meta_channel),
            "dwi_preprocessed_metadata"
        ).map{ it.flatten() }
        dwi_mask_channel = rename_processed_dwi_mask(
            collect_paths(dwi_mask_channel),
            "mask_preprocessed"
        ).map{ it.flatten() }

        wm_segmentation_channel = Channel.empty()
        if ( params.generate_wm_segmentation ) {
            segment_wm_wkf(dwi_channel, dwi_mask_channel)
            wm_segmentation_channel = segment_wm_wkf.out.segmentation
        }

    emit:
        t1 = t1_channel
        dwi = dwi_channel
        mask = dwi_mask_channel
        pvf = pvf_channel
        tissue_masks = tissue_mask_channel
        safe_wm_mask = safe_wm_mask_channel
        wm_segmentation = wm_segmentation_channel
        metadata = meta_channel
}

workflow validate_gradients_wkf {
    take:
        dwi_channel
        mask_channel
    main:
        scil_compute_dti_fa(dwi_channel.join(mask_channel), "preprocess", "preprocess", "false")
        validate_gradients(
            dwi_channel.map{ [it[0], it[3]]}
                .join(scil_compute_dti_fa.out.main_peak)
                .join(scil_compute_dti_fa.out.fa)
                .join(mask_channel)
                .map{ it + [""] },
            "preprocess"
        )
        dwi_channel = dwi_channel.map{ it[0..-2] }.join(validate_gradients.out.bvecs)
    emit:
        dwi = dwi_channel
}

workflow t1_preprocess_wkf {
    take:
        t1_channel
        mask_channel
    main:
        def ref_id_channel = t1_channel.map{ [it[0]] }
        mask_channel = fill_missing_datapoints(mask_channel, ref_id_channel, 1, [""])

        if ( params.denoise_t1 ) {
            if ( params.nlmeans_t1 ) {
                nlmeans_denoise(t1_channel.join(mask_channel).map{ it + [""] }, "preprocess", "true")
                t1_channel = nlmeans_denoise.out.image
            }
            else {
                ants_gaussian_denoise(t1_channel, "preprocess")
                t1_channel = ants_gaussian_denoise.out.image
            }
        }

        if ( params.t1_intensity_normalization ) {
            n4_denoise_wkf(
                t1_channel,
                t1_channel.map{ [it[0], ""] },
                mask_channel,
                t1_channel.map{ [it[0], ""] },
                params.t1_n4_normalization_config,
                true
            )
            t1_channel = n4_denoise_wkf.out.image
        }
    emit:
        t1 = t1_channel
        mask = mask_channel
}
