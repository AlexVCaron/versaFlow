#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Preprocess workflow parameters
params.gaussian_noise_correction = true
params.gibbs_ringing_correction = true
params.t1mask2dwi_registration = true
params.topup_correction = true
params.eddy_correction = true
params.dwi_intensity_normalization = true
params.resample_data = true
params.register_t12b0_denoised = true
params.register_syn_t12b0 = false
params.generate_tissue_segmentation = false
params.generate_wm_segmentation = true
params.raw_to_processed_space = false

// T1 preprocess workflow parameters
params.denoise_t1 = true
params.nlmeans_t1 = true
params.t1_intensity_normalization = true

params.b02t1_mask_registration_config = file("$projectDir/.config/b02t1_mask_registration_config.py")
params.t1_mask_to_topup_b0_registration_config = file("$projectDir/.config/t1_mask_to_topup_b0_registration_config.py")
params.ants_transform_base_config = file("$projectDir/.config/ants_transform_base_config.py")
params.extract_mean_b0_base_config = file("$projectDir/.config/extract_mean_b0_base_config.py")
params.dwi_n4_normalization_config = file("$projectDir/.config/dwi_n4_normalization_config.py")
params.t1_n4_normalization_config = file("$projectDir/.config/t1_n4_normalization_config.py")
params.b0_to_b0_normalization_config = file("$projectDir/.config/b0_to_b0_normalization_config.py")

include {
    merge_channels_non_blocking; replace_dwi_file; uniformize_naming; exclude_missing_datapoints; fill_missing_datapoints; filter_datapoints
} from '../modules/functions.nf'
include { extract_b0 as dwi_b0; extract_b0 as extract_topup_b0; extract_b0 as extract_b0_preprocessed } from '../modules/processes/preprocess.nf'
include { scil_compute_dti_fa } from '../modules/processes/measure.nf'
include {
    ants_transform as ants_transform_base_t1; ants_transform as ants_transform_base_dwi; ants_transform as ants_transform_syn_t1;
    ants_transform as ants_transform_syn_dwi; ants_transform as ants_transform_base_wm; ants_transform as ants_transform_base_gm;
    ants_transform as ants_transform_base_csf; ants_transform as ants_transform_syn_wm; ants_transform as ants_transform_syn_gm;
    ants_transform as ants_transform_syn_csf; ants_transform as ants_transform_base_raw_t1; ants_transform as ants_transform_syn_raw_t1
} from '../modules/processes/register.nf'
include {
    convert_float_to_integer as convert_wm_segmentation; convert_float_to_integer as convert_gm_segmentation;
    convert_float_to_integer as convert_csf_segmentation; convert_float_to_integer as dwi_mask_convert_datatype;
    convert_float_to_integer as t1_mask_convert_datatype;
    crop_image as crop_dwi; crop_image as crop_t1; crop_image as crop_wm; crop_image as crop_gm; crop_image as crop_csf;
    crop_image as crop_raw_dwi; crop_image as crop_raw_t1;
    apply_mask as apply_mask_to_t1_for_reg; apply_mask as apply_mask_to_b0_for_reg;
    dilate_mask as dilate_t1_mask; dilate_mask as dilate_b0_mask;
    bet_mask;fit_bounding_box; merge_masks; check_odd_dimensions; pvf_to_mask
} from '../modules/processes/utils.nf'
include { gibbs_removal as dwi_gibbs_removal; gibbs_removal as rev_gibbs_removal; nlmeans_denoise; ants_gaussian_denoise; normalize_inter_b0 } from '../modules/processes/denoise.nf'
include {
    scilpy_resample as scilpy_resample_wm; scilpy_resample as scilpy_resample_gm; scilpy_resample as scilpy_resample_csf;
    scilpy_resample as scilpy_resample_t1; scilpy_resample as scilpy_resample_dwi;
    scilpy_resample as scilpy_resample_raw_dwi; scilpy_resample as scilpy_resample_raw_t1
} from '../modules/processes/upsample.nf'
include {
    registration_wkf as dwi_mask_registration_wkf; registration_wkf as t1_mask_registration_wkf;
    dwi_denoise_wkf; dwi_denoise_wkf as rev_denoise_wkf; n4_denoise_wkf;
    squash_wkf; squash_wkf as squash_raw_wkf; topup_wkf; apply_topup_wkf; apply_topup_wkf as raw_apply_topup_wkf; eddy_wkf
} from "../modules/workflows/preprocess.nf"
include { t12b0_registration as mask_registration_wkf; t12b0_registration as t1_registration_wkf } from '../modules/workflows/t1_registration.nf'
include { segment_nmt_wkf; segment_wm_wkf } from '../modules/workflows/segment.nf'

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
        def ref_id_channel = dwi_channel.map{ [it[0]] }

        raw_dwi_channel = dwi_channel
        raw_rev_channel = rev_channel
        raw_t1_channel = t1_channel
        raw_t1_mask_channel = t1_mask_channel
        raw_meta_channel = meta_channel
        raw_rev_meta_channel = rev_meta_channel

        if ( params.gaussian_noise_correction ) {
            dwi_denoise_wkf(dwi_channel, dwi_mask_channel, meta_channel, "true")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_denoise_wkf.out.image)
            meta_channel = dwi_denoise_wkf.out.metadata

            rev_denoise_wkf(rev_channel, dwi_mask_channel, rev_meta_channel, "false")
            rev_channel = replace_dwi_file(rev_channel, rev_denoise_wkf.out.image)
            rev_meta_channel = rev_denoise_wkf.out.metadata
        }

        if ( params.gibbs_ringing_correction ) {
            dwi_gibbs_removal(dwi_channel.map{ it.subList(0, 2) }.join(meta_channel), "preprocess", "true")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_gibbs_removal.out.image)
            meta_channel = dwi_gibbs_removal.out.metadata

            rev_gibbs_removal(
                exclude_missing_datapoints(rev_channel.map{ it.subList(0, 2) }.join(rev_meta_channel), 1, ""),
                "preprocess", "false"
            )
            rev_channel = replace_dwi_file(
                rev_channel,
                fill_missing_datapoints(rev_gibbs_removal.out.image, ref_id_channel, 1, [""])
            )
            rev_meta_channel = fill_missing_datapoints(rev_gibbs_removal.out.metadata, ref_id_channel, 1, [""])
        }

        if ( params.normalize_inter_b0 ) {
            normalize_inter_b0(
                dwi_channel.map{ it.subList(0, 3) }.join(rev_channel.map{ it.subList(0, 3) }).join(meta_channel).join(rev_meta_channel),
                "preprocess",
                params.b0_to_b0_normalization_config
            )
            dwi_channel = replace_dwi_file(dwi_channel, normalize_inter_b0.out.dwi)
            meta_channel = normalize_inter_b0.out.dwi_metadata
            rev_channel = replace_dwi_file(rev_channel, fill_missing_datapoints(normalize_inter_b0.out.rev, ref_id_channel, 1, [""]))
            rev_meta_channel = fill_missing_datapoints(normalize_inter_b0.out.rev_metadata, ref_id_channel, 1, [""])
        }

        squash_wkf(dwi_channel, rev_channel, meta_channel.join(rev_meta_channel), "")
        dwi_channel = squash_wkf.out.dwi
        rev_channel = squash_wkf.out.rev
        meta_channel = squash_wkf.out.metadata

        topup_corrected_dwi = Channel.empty()
        excluded_dwi_channel = dwi_channel
        if ( params.topup_correction ) {
            check_odd_dimensions(dwi_channel.join(rev_channel).join(dwi_mask_channel).join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess")
            dwi_channel = check_odd_dimensions.out.dwi
            rev_channel = fill_missing_datapoints(
                check_odd_dimensions.out.rev.join(check_odd_dimensions.out.rev_bval_bvec, remainder: true),
                ref_id_channel,
                2, ["", ""]
            )
            dwi_mask_channel = fill_missing_datapoints(
                check_odd_dimensions.out.mask,
                ref_id_channel,
                1, [""]
            )

            meta_channel = check_odd_dimensions.out.metadata.map{ it.flatten() }

            topup_wkf(dwi_channel, rev_channel, meta_channel)

            topup2eddy_channel = topup_wkf.out.param.join(topup_wkf.out.prefix).join(topup_wkf.out.topup.map{ [it[0], it.subList(1, it.size())] })

            dwi2topup_channel = uniformize_naming(topup_wkf.out.topupable_indexes.join(dwi_channel), "dwi_to_topup", "false", "false")
            rev2topup_channel = uniformize_naming(topup_wkf.out.topupable_indexes.join(rev_channel), "dwi_to_topup_rev", "false", "false")
            meta2topup_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(0..<it[1].size()).step(2)] }, "dwi_to_topup_metadata", "false", "false")
            rev_meta2topup_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(1..<it[1].size()).step(2)] }, "dwi_to_topup_rev_metadata", "false", "false")

            apply_topup_wkf(dwi2topup_channel, rev2topup_channel, topup2eddy_channel, meta2topup_channel.join(rev_meta2topup_channel).map{ [it[0], it.subList(1, it.size())] }, "")
            topup_corrected_dwi = uniformize_naming(apply_topup_wkf.out.dwi, "topup_corrected", "false", "false")
            topup_corrected_dwi_meta = uniformize_naming(apply_topup_wkf.out.metadata, "topup_corrected_metadata", "false", "false")

            excluded_dwi_channel = topup_wkf.out.excluded_dwi
            excluded_meta_channel = topup_wkf.out.excluded_dwi_metadata

            extract_topup_b0(
                topup_corrected_dwi.mix(excluded_dwi_channel).map{ it.subList(0, 3) }.join(
                    topup_corrected_dwi_meta.mix(excluded_meta_channel).map{ [it[0], it.subList(1, it.size())] }
                ),
                "preprocess",
                "false",
                params.extract_mean_b0_base_config
            )
            b0_channel = extract_topup_b0.out.b0
            b0_metadata = extract_topup_b0.out.metadata

            if ( !params.eddy_correction ) {
                dwi_channel = topup_corrected_dwi.mix(excluded_dwi_channel)
                meta_channel = topup_corrected_dwi_meta.mix(excluded_meta_channel)
            }
            else {
                topup2eddy_channel = fill_missing_datapoints(topup2eddy_channel, ref_id_channel, 1, ["", "", []])
            }

            if ( params.raw_to_processed_space ) {
                squash_raw_wkf(raw_dwi_channel, raw_rev_channel, raw_meta_channel.join(raw_rev_meta_channel), "raw")
                raw_dwi_channel = uniformize_naming(raw_dwi_channel, "raw", "false", "false")
                raw_rev_channel = uniformize_naming(raw_rev_channel, "raw", "false", "false")
                raw_meta_channel = uniformize_naming(meta2topup_channel, "raw_metadata", "false", "false")
                raw_rev_meta_channel = uniformize_naming(rev_meta2topup_channel, "raw_metadata", "false", "false")
                raw_apply_topup_wkf(
                    topup_wkf.out.topupable_indexes.join(raw_dwi_channel),
                    topup_wkf.out.topupable_indexes.join(raw_rev_channel),
                    topup2eddy_channel, raw_meta_channel.join(raw_rev_meta_channel).map{ [it[0], it.subList(1, it.size())] },
                    "raw"
                )
                raw_dwi_channel = excluded_dwi_channel.map{ [it[0]] }.join(raw_dwi_channel).mix(raw_apply_topup_wkf.out.dwi)
                raw_meta_channel = excluded_dwi_channel.map{ [it[0]] }.join(raw_meta_channel).mix(raw_apply_topup_wkf.out.metadata)
            }
        }
        else {
            dwi_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess", "false", params.extract_mean_b0_base_config)
            b0_channel = dwi_b0.out.b0
            b0_metadata = dwi_b0.out.metadata
            topup2eddy_channel = ref_id_channel.map{ it + ["", "", []] }
        }

        empty_dwi_mask_ids = filter_datapoints(dwi_mask_channel, { it[1] == "" }).map{ [it[0]] }
        dwi_mask_channel = exclude_missing_datapoints(dwi_mask_channel, 1, "")
        dwi_mask_channel = dwi_mask_channel.mix(bet_mask(empty_dwi_mask_ids.join(b0_channel), "preprocess", "${!params.t1mask2dwi_registration}"))

        if ( params.t1mask2dwi_registration ) {
            existing_t1_mask_ids = exclude_missing_datapoints(t1_mask_channel, 1, "").map{ [it[0]] }
            absent_t1_mask_ids = filter_datapoints(t1_mask_channel, { it[1] == "" }).map{ [it[0]] }

            apply_mask_to_b0_for_reg(existing_t1_mask_ids.join(b0_channel).join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
            reg_b0_channel = apply_mask_to_b0_for_reg.out.image

            scil_compute_dti_fa(existing_t1_mask_ids.join(topup_corrected_dwi.mix(excluded_dwi_channel)).join(dwi_mask_channel), "preprocess", "preprocess")

            apply_mask_to_t1_for_reg(existing_t1_mask_ids.join(t1_channel).join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")
            reg_t1_channel = apply_mask_to_t1_for_reg.out.image

            reg_t1_mask_channel = dilate_t1_mask(existing_t1_mask_ids.join(t1_mask_channel), 3, "preprocess")
            reg_b0_mask_channel = dilate_b0_mask(existing_t1_mask_ids.join(dwi_mask_channel), 3, "preprocess")

            t1_mask_registration_wkf(
                merge_channels_non_blocking(reg_b0_channel, scil_compute_dti_fa.out.fa),
                reg_t1_channel.map{ [it[0], [it[1]]] },
                existing_t1_mask_ids.join(t1_mask_channel),
                reg_b0_mask_channel.join(reg_t1_mask_channel).map{ [it[0], [it[1], it[2]]] },
                null,
                meta_channel.map{ [it[0], it[1], ""] },
                "",
                !params.register_t12b0_denoised,
                "", "mask",
                params.t1_mask_to_topup_b0_registration_config,
                null
            )

            t1_mask_convert_datatype(t1_mask_registration_wkf.out.image, "uint8", "preprocess", !params.register_t12b0_denoised, "mask", "")
            dwi_mask_channel = t1_mask_convert_datatype.out.image.mix(absent_t1_mask_ids.join(dwi_mask_channel))
        }

        if ( params.eddy_correction ) {
            eddy_wkf(dwi_channel, dwi_mask_channel, topup2eddy_channel, b0_channel, rev_channel, meta_channel)

            dwi_channel = eddy_wkf.out.dwi.join(eddy_wkf.out.bval).join(eddy_wkf.out.bvec)
            meta_channel = eddy_wkf.out.metadata
        }

        if ( params.dwi_intensity_normalization ) {
            n4_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, b0_channel, dwi_mask_channel, meta_channel, params.dwi_n4_normalization_config)
            dwi_channel = replace_dwi_file(dwi_channel, n4_denoise_wkf.out.image)
            meta_channel = n4_denoise_wkf.out.metadata
        }

        absent_t1_mask_ids = filter_datapoints(t1_mask_channel, { it[1] == "" }).map{ [it[0]] }

        dwi_mask_registration_wkf(
            absent_t1_mask_ids.join(t1_channel).map{ [it[0], [it[1]]] },
            absent_t1_mask_ids.join(b0_channel).map{ [it[0], [it[1]]] },
            absent_t1_mask_ids.join(dwi_mask_channel),
            null,
            null,
            absent_t1_mask_ids.join(b0_metadata).map{ it.subList(0, 2) + [""] },
            "",
            false, "", "",
            params.b02t1_mask_registration_config,
            null
        )

        dwi_mask_convert_datatype(dwi_mask_registration_wkf.out.image, "uint8", "preprocess", false, "", "")
        t1_mask_channel = exclude_missing_datapoints(t1_mask_channel, 1, "").mix(dwi_mask_convert_datatype.out.image)
        raw_t1_mask_channel = t1_mask_channel

        t1_preprocess_wkf(t1_channel.map{ it.subList(0, 2) }, t1_mask_channel)
        t1_channel = t1_preprocess_wkf.out.t1
        t1_mask_channel = t1_preprocess_wkf.out.mask

        if ( params.resample_data ) {
            scilpy_resample_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).join(meta_channel), "preprocess", "lin", true, "mask", "")

            pvf_to_resample = pvf_channel.filter{ !it[1].isEmpty() }

            scilpy_resample_wm(pvf_to_resample.map{ [it[0], it[1][0]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn", false, "", "")
            scilpy_resample_gm(pvf_to_resample.map{ [it[0], it[1][1]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn", false, "", "")
            scilpy_resample_csf(pvf_to_resample.map{ [it[0], it[1][2]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn", false, "", "")
            pvf_channel = scilpy_resample_wm.out.image.join(scilpy_resample_gm.out.image).join(scilpy_resample_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(pvf_channel.filter{ it[1].isEmpty() })

            if ( params.raw_to_processed_space ) {
                scilpy_resample_raw_dwi(raw_dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).join(raw_meta_channel), "preprocess", "lin", true, "mask", "raw")
                raw_dwi_channel = replace_dwi_file(raw_dwi_channel, scilpy_resample_raw_dwi.out.image)
                raw_meta_channel = scilpy_resample_raw_dwi.out.metadata

                scilpy_resample_raw_t1(raw_t1_channel.join(raw_t1_mask_channel).map{ it + [""] }, "preprocess", "lin", false, "", "raw")
                raw_t1_mask_channel = scilpy_resample_raw_t1.out.image
                raw_t1_mask_channel = t1_mask_channel
            }

            dwi_channel = replace_dwi_file(dwi_channel, scilpy_resample_dwi.out.image)
            dwi_mask_channel = scilpy_resample_dwi.out.mask
            meta_channel = scilpy_resample_dwi.out.metadata
        }

        if ( params.register_t12b0_denoised ) {
            t1_registration_wkf(
                dwi_channel,
                t1_channel,
                t1_mask_channel,
                dwi_mask_channel,
                meta_channel,
                false,
                true
            )
            t1_channel = t1_registration_wkf.out.t1
            ants_transform_base_t1(
                t1_mask_channel.join(t1_channel).join(t1_registration_wkf.out.transform_base.map{ [it[0], it[2]] }).map{ it + ["", ""] },
                "preprocess",
                "",
                "${!params.register_syn_t12b0}",
                "",
                params.ants_transform_base_config
            )
            ants_transform_base_dwi(
                t1_mask_channel.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] },
                "preprocess",
                "",
                "${!params.register_syn_t12b0}",
                "mask",
                params.ants_transform_base_config
            )
            t1_mask_channel = ants_transform_base_t1.out.image
            dwi_mask_channel = uniformize_naming(ants_transform_base_dwi.out.image, "dwi_mask", "false", "true")

            pvf_to_register = pvf_channel.filter{ !it[1].isEmpty() }
            ants_transform_base_wm(
                pvf_to_register.map{ [it[0], it[1][0]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] },
                "preprocess", "", "true", "",
                params.ants_transform_base_config
            )
            ants_transform_base_gm(
                pvf_to_register.map{ [it[0], it[1][1]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] },
                "preprocess", "", "true", "",
                params.ants_transform_base_config
            )
            ants_transform_base_csf(
                pvf_to_register.map{ [it[0], it[1][2]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] },
                "preprocess", "", "true", "",
                params.ants_transform_base_config
            )
            pvf_channel = ants_transform_base_wm.out.image.join(ants_transform_base_gm.out.image).join(ants_transform_base_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(pvf_channel.filter{ it[1].isEmpty() })

            if ( params.raw_to_processed_space ) {
                ants_transform_base_raw_t1(
                    raw_t1_channel.join(t1_channel).join(t1_registration_wkf.out.transform_base.map{ [it[0], it[2]] }).map{ it + ["", ""] },
                    "preprocess",
                    "raw",
                    "${!params.register_syn_t12b0}",
                    "",
                    params.ants_transform_base_config
                )
                raw_t1_channel = ants_transform_base_raw_t1.out.image
                raw_t1_mask_channel = t1_mask_channel
            }

            if ( params.register_syn_t12b0 ) {
                ants_transform_syn_t1(
                    t1_mask_channel.join(t1_channel).join(t1_registration_wkf.out.transform_syn.map{ [it[0], it[2]] }).map{ it + ["", ""] },
                    "preprocess", "", "false", "",
                    params.ants_transform_base_config
                )
                ants_transform_syn_dwi(
                    t1_mask_channel.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] },
                    "preprocess", "", "true", "mask",
                    params.ants_transform_base_config
                )
                t1_mask_channel = ants_transform_syn_t1.out.image
                dwi_mask_channel = uniformize_naming(ants_transform_syn_dwi.out.image, "dwi_mask", "false", "true")

                pvf_to_register = pvf_channel.filter{ !it[1].isEmpty() }
                ants_transform_syn_wm(
                    pvf_to_register.map{ [it[0], it[1][0]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] },
                    "preprocess", "", "true", "",
                    params.ants_transform_base_config
                )
                ants_transform_syn_gm(
                    pvf_to_register.map{ [it[0], it[1][1]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] },
                    "preprocess", "", "true", "",
                    params.ants_transform_base_config
                )
                ants_transform_syn_csf(
                    pvf_to_register.map{ [it[0], it[1][2]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] },
                    "preprocess", "", "true", "",
                    params.ants_transform_base_config
                )
                pvf_channel = ants_transform_syn_wm.out.image.join(ants_transform_syn_gm.out.image).join(ants_transform_syn_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(pvf_channel.filter{ it[1].isEmpty() })

                if ( params.raw_to_processed_space ) {
                    ants_transform_syn_raw_t1(
                        raw_t1_channel.join(t1_channel).join(t1_registration_wkf.out.transform_base.map{ [it[0], it[2]] }).map{ it + ["", ""] },
                        "preprocess",
                        "raw",
                        "true",
                        "",
                        params.ants_transform_base_config
                    )
                    raw_t1_channel = ants_transform_syn_raw_t1.out.image
                    raw_t1_mask_channel = t1_mask_channel
                }
            }
        }

        crop_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).map{ it + [""] }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess", true, "mask", "")
        dwi_bbox_channel = crop_dwi.out.bbox
        fit_bounding_box(t1_channel.join(dwi_channel.map{ it.subList(0, 2) }).join(dwi_bbox_channel), "preprocess")
        dwi_bbox_channel = fit_bounding_box.out.bbox
        crop_t1(t1_channel.join(t1_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess", false, "", "")

        if ( params.raw_to_processed_space ) {
            crop_raw_dwi(raw_dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).map{ it + [""] }.join(raw_meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess", true, "mask", "raw")
            raw_dwi_channel = replace_dwi_file(raw_dwi_channel, crop_raw_dwi.out.image)
            raw_meta_channel = crop_raw_dwi.out.metadata
            crop_raw_t1(raw_t1_channel.join(raw_t1_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess", false, "", "")
            raw_t1_channel = crop_raw_t1.out.image
            raw_t1_mask_channel = crop_raw_dwi.out.mask
        }

        pvf_to_crop = pvf_channel.filter{ !it[1].isEmpty() }
        crop_wm(pvf_to_crop.map{ [it[0], it[1][0]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess", true, "", "segmentation")
        crop_gm(pvf_to_crop.map{ [it[0], it[1][1]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess", true, "", "segmentation")
        crop_csf(pvf_to_crop.map{ [it[0], it[1][2]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess", true, "", "segmentation")

        pvf_channel = crop_wm.out.image.join(crop_gm.out.image).join(crop_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(pvf_channel.filter{ it[1].isEmpty() })

        dwi_channel = replace_dwi_file(dwi_channel, crop_dwi.out.image)
        meta_channel = crop_dwi.out.metadata
        dwi_mask_channel = crop_dwi.out.mask
        t1_channel = crop_t1.out.image
        t1_mask_channel = crop_dwi.out.mask

        seg_to_tissue_masks = pvf_channel.filter{ !it[1].isEmpty() }
        pvf_to_mask(seg_to_tissue_masks.join(dwi_mask_channel), "preprocess", "segmentation")
        tissue_masks = pvf_to_mask.out.wm_mask.join(
            pvf_to_mask.out.gm_mask
        ).join(
            pvf_to_mask.out.csf_mask
        ).mix(
            pvf_channel.filter{ it[1].isEmpty() }.map{ [it[0], "", "", ""] }
        )
        safe_wm_mask = pvf_to_mask.out.safe_wm_mask.mix(pvf_channel.filter{ it[1].isEmpty() }.map{ [it[0], ""] })

        if ( params.generate_tissue_segmentation ) {
            empty_segmentations = pvf_channel.filter{ it[1].isEmpty() }.map{ [it[0]] }
            segment_nmt_wkf(empty_segmentations.join(t1_channel), empty_segmentations.join(t1_mask_channel))
            pvf_channel = pvf_channel.filter{ !it[1].isEmpty() }.mix(segment_nmt_wkf.out.volume_fractions.map{ [it[0], it[1].reverse()] })
            tissue_masks = segment_nmt_wkf.out.tissue_masks.mix(tissue_masks.filter{ it[1] })
            safe_wm_mask = segment_nmt_wkf.out.safe_wm_mask.mix(safe_wm_mask.filter{ it[1] })
        }

        wm_segmentation = Channel.empty()
        if ( params.generate_wm_segmentation ) {
            segment_wm_wkf(dwi_channel, dwi_mask_channel)
            wm_segmentation = segment_wm_wkf.out.segmentation
        }

        extract_b0_preprocessed(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "preprocess", "true", params.extract_mean_b0_base_config)

        dwi_channel = uniformize_naming(dwi_channel.map{ it.subList(0, 4) }, "dwi_preprocessed", "false", "false")
        meta_channel = uniformize_naming(meta_channel, "dwi_preprocessed_metadata", "false", "false")
        dwi_mask_channel = uniformize_naming(dwi_mask_channel, "mask_preprocessed", "false", "false")
    emit:
        t1 = t1_channel
        dwi = dwi_channel
        mask = dwi_mask_channel
        pvf = pvf_channel
        tissue_masks = tissue_masks
        safe_wm_mask = safe_wm_mask
        wm_segmentation = wm_segmentation
        metadata = meta_channel
}

workflow t1_preprocess_wkf {
    take:
        t1_channel
        mask_channel
    main:
        def ref_id_channel = t1_channel.map{ [it[0]] }
        if ( params.denoise_t1 ) {
            if ( params.nlmeans_t1 ) {
                nlmeans_denoise(t1_channel.join(fill_missing_datapoints(mask_channel, ref_id_channel, 1, [""])).map{ it + [""] }, "preprocess", "true")
                t1_channel = nlmeans_denoise.out.image
            }
            else {
                ants_gaussian_denoise(t1_channel, "preprocess")
                t1_channel = ants_gaussian_denoise.out.image
            }
        }

        if ( params.t1_intensity_normalization ) {
            n4_denoise_wkf(t1_channel, Channel.empty(), Channel.empty(), Channel.empty(), params.t1_n4_normalization_config)
            t1_channel = n4_denoise_wkf.out.image
        }

        if ( params.resample_data ) {
            scilpy_resample_t1(t1_channel.join(mask_channel).map{ it + [""] }, "preprocess", "lin", false, "", "")
            t1_channel = scilpy_resample_t1.out.image
            mask_channel = scilpy_resample_t1.out.mask
        }
    emit:
        t1 = t1_channel
        mask = mask_channel
}
