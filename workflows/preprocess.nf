#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Preprocess workflow parameters
params.gaussian_noise_correction = true
params.rev_is_b0 = true
params.has_reverse = true
params.gibbs_ringing_correction = true
params.t1mask2dwi_registration = true
params.masked_t1 = true
params.masked_dwi = false
params.topup_correction = true
params.eddy_correction = true
params.do_bet_mask = false
params.intensity_normalization = true
params.resample_data = true
params.register_repetitions = true
params.register_t12b0_denoised = true
params.register_syn_t12b0 = false
params.msmt_odf = false
params.seg_on_t1 = true
params.generate_segmentation = false

// T1 preprocess workflow parameters
params.denoise_t1 = true
params.nlmeans_t1 = true
// params.intensity_normalization = true
// params.resample_data = true

params.b02t1_mask_registration_config = file("$projectDir/.config/b02t1_mask_registration_config.py")
params.t1_mask_to_topup_b0_registration_config = file("$projectDir/.config/t1_mask_to_topup_b0_registration_config.py")
params.ants_transform_base_config = file("$projectDir/.config/ants_transform_base_config.py")
params.extract_mean_b0_base_config = file("$projectDir/.config/extract_mean_b0_base_config.py")
params.dwi_n4_normalization_config = file("$projectDir/.config/dwi_n4_normalization_config.py")
params.t1_n4_normalization_config = file("$projectDir/.config/t1_n4_normalization_config.py")
params.b0_to_b0_normalization_config = file("$projectDir/.config/b0_to_b0_normalization_config.py")

include { merge_channels_non_blocking; map_optional; opt_channel; replace_dwi_file; uniformize_naming; merge_repetitions; is_data } from '../modules/functions.nf'
include { extract_b0 as dwi_b0; extract_b0 as extract_topup_b0 } from '../modules/processes/preprocess.nf'
include { scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { ants_transform as ants_transform_base_t1; ants_transform as ants_transform_base_dwi; ants_transform as ants_transform_syn_t1; ants_transform as ants_transform_syn_dwi; ants_transform as ants_transform_base_wm; ants_transform as ants_transform_base_gm; ants_transform as ants_transform_base_csf; ants_transform as ants_transform_syn_wm; ants_transform as ants_transform_syn_gm; ants_transform as ants_transform_syn_csf } from '../modules/processes/register.nf'
include { convert_datatype; convert_datatype as convert_wm_segmentation; convert_datatype as convert_gm_segmentation; convert_datatype as convert_csf_segmentation; convert_datatype as topup_convert_datatype; convert_datatype as rev_convert_datatype; convert_datatype as t1_mask_convert_datatype; bet_mask; bet_mask as rev_bet_mask; bet_mask as bet_mask_topup; crop_image as crop_dwi; crop_image as crop_t1; crop_image as crop_wm; crop_image as crop_gm; crop_image as crop_csf; fit_bounding_box; average as average_t1; average as average_b0; merge_masks; apply_mask as apply_mask_to_t1_for_reg; apply_mask as apply_mask_to_b0_for_reg; dilate_mask as dilate_t1_mask; dilate_mask as dilate_b0_mask } from '../modules/processes/utils.nf'
include { gibbs_removal as dwi_gibbs_removal; gibbs_removal as rev_gibbs_removal; nlmeans_denoise; ants_gaussian_denoise; normalize_inter_b0 } from '../modules/processes/denoise.nf'
include { scilpy_resample as scilpy_resample_wm; scilpy_resample as scilpy_resample_gm; scilpy_resample as scilpy_resample_csf; scilpy_resample as scilpy_resample_t1; scilpy_resample_on_ref as scilpy_resample_t1_mask; scilpy_resample as scilpy_resample_dwi; scilpy_resample_on_ref as scilpy_resample_mask } from '../modules/processes/upsample.nf'
include { dwi_denoise_wkf; dwi_denoise_wkf as rev_denoise_wkf; squash_wkf; registration_wkf as topup_mask_registration_wkf; registration_wkf as t1_mask_registration_wkf; topup_wkf; eddy_wkf; apply_topup_wkf; n4_denoise_wkf } from "../modules/workflows/preprocess.nf"
include { cat_dwi_repetitions_wkf; cat_dwi_repetitions_wkf as cat_rev_repetitions_wkf; register_dwi_repetitions_wkf as pre_register_dwi_repetitions_wkf; register_t1_repetitions_wkf } from '../modules/workflows/repetitions.nf'
include { t12b0_registration as mask_registration_wkf; t12b0_registration as rev_mask_registration_wkf; t12b0_registration as t1_registration_wkf } from '../modules/workflows/t1_registration.nf'
include { segment_nmt } from '../modules/workflows/segment.nf'

workflow preprocess_wkf {
    take:
        dwi_channel
        rev_channel
        t1_channel
        seg_channel
        meta_channel
        rev_meta_channel
        dwi_mask_channel
        t1_mask_channel
    main:
        topup2eddy_channel = opt_channel()
        topup2eddy_b0_channel = opt_channel()

        if ( params.gaussian_noise_correction ) {
            dwi_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, dwi_mask_channel, meta_channel)
            dwi_channel = replace_dwi_file(dwi_channel, dwi_denoise_wkf.out.image)
            meta_channel = dwi_denoise_wkf.out.metadata

            if ( params.has_reverse && !params.rev_is_b0 ) {
                rev_denoise_wkf(rev_channel.map{ it.subList(0, 2) }, dwi_mask_channel, rev_meta_channel)
                rev_channel = replace_dwi_file(rev_channel, rev_denoise_wkf.out.image)
                rev_meta_channel = rev_denoise_wkf.out.metadata
            }
        }

        if ( params.gibbs_ringing_correction ) {
            dwi_gibbs_removal(dwi_channel.map{ it.subList(0, 2) }.join(meta_channel), "preprocess")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_gibbs_removal.out.image)
            meta_channel = dwi_gibbs_removal.out.metadata

            if ( params.has_reverse && !params.rev_is_b0 ) {
                rev_gibbs_removal(rev_channel.map{ it.subList(0, 2) }.join(rev_meta_channel), "preprocess")
                rev_channel = replace_dwi_file(rev_channel, rev_gibbs_removal.out.image)
                rev_meta_channel = rev_gibbs_removal.out.metadata
            }
        }

        if ( params.normalize_inter_b0 ) {
            if ( !params.has_reverse ) {
                to_normalize = dwi_channel.map{ it.subList(0, 3) + ["", ""] }.join(meta_channel.map{ it + [""] })
            } else if ( params.rev_is_b0 ) {
                to_normalize = dwi_channel.map{ it.subList(0, 3) }.join(rev_channel.map{ it + [""] }).join(meta_channel).join(rev_meta_channel)
            } else {
                to_normalize = dwi_channel.map{ it.subList(0, 3) }.join(rev_channel.map{ it.subList(0, 3) }).join(meta_channel).join(rev_meta_channel)
            }
            normalize_inter_b0(to_normalize, "preprocess", params.b0_to_b0_normalization_config)
            dwi_channel = replace_dwi_file(dwi_channel, normalize_inter_b0.out.dwi)
            meta_channel = normalize_inter_b0.out.dwi_metadata
            rev_meta_channel = normalize_inter_b0.out.rev_metadata
            if ( params.has_reverse ) {
                if (params.rev_is_b0) {
                    rev_channel = normalize_inter_b0.out.rev
                } else {
                    rev_channel = replace_dwi_file(rev_channel, normalize_inter_b0.out.rev)
                }
            }
        }

        if ( params.register_repetitions ) {
            pre_register_dwi_repetitions_wkf(
                dwi_channel,
                rev_channel,
                meta_channel,
                rev_meta_channel
            )
            register_t1_repetitions_wkf(
                t1_channel,
                params.masked_t1 ? t1_mask_channel : null
            )
            t1_channel = register_t1_repetitions_wkf.out.t1
            t1_mask_channel = params.masked_t1 ? register_t1_repetitions_wkf.out.mask : null
            dwi_channel = pre_register_dwi_repetitions_wkf.out.dwi
            rev_channel = pre_register_dwi_repetitions_wkf.out.rev
            meta_channel = pre_register_dwi_repetitions_wkf.out.metadata
            rev_meta_channel  = pre_register_dwi_repetitions_wkf.out.rev_metadata
        }

        squash_metadata = meta_channel.join(rev_meta_channel, remainder: true)
        squash_wkf(dwi_channel, rev_channel, squash_metadata)
        dwi_channel = squash_wkf.out.dwi
        rev_channel = squash_wkf.out.rev
        meta_channel = squash_wkf.out.metadata

        if ( params.has_reverse && params.topup_correction ) {
            topup_wkf(dwi_channel, rev_channel, meta_channel)

            topup2eddy_channel = topup_wkf.out.param.join(topup_wkf.out.prefix).join(topup_wkf.out.topup.map{ [it[0], it.subList(1, it.size())] })

            dwi2topup_channel = uniformize_naming(dwi_channel, "dwi_to_topup", "false", "false")
            rev2topup_channel = uniformize_naming(rev_channel, "dwi_to_topup_rev", "false", "false")
            meta2topup_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(0..<it[1].size()).step(2)] }, "dwi_to_topup_metadata", "false", "false")
            rev_meta2topup_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(1..<it[1].size()).step(2)] }, "dwi_to_topup_rev_metadata", "false", "false")
            apply_topup_wkf(dwi2topup_channel, rev2topup_channel, topup2eddy_channel, meta2topup_channel.join(rev_meta2topup_channel).map{ [it[0], it.subList(1, it.size())] })
            topup_corrected_dwi = uniformize_naming(apply_topup_wkf.out.dwi, "topup_corrected", "false", "false")
            topup_corrected_dwi_meta = uniformize_naming(apply_topup_wkf.out.metadata, "topup_corrected_metadata", "false", "false")

            extract_topup_b0(topup_corrected_dwi.map{ it.subList(0, 3) }.join(topup_corrected_dwi_meta.map{ [it[0], it.subList(1, it.size())] }), "preprocess", params.extract_mean_b0_base_config)
            b0_channel = extract_topup_b0.out.b0
            b0_metadata = extract_topup_b0.out.metadata

            if ( !params.eddy_correction ) {
                dwi_channel = topup_corrected_dwi
                meta_channel = topup_corrected_dwi_meta
            }
        }
        else {
            dwi_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess", params.extract_mean_b0_base_config)
            b0_channel = dwi_b0.out.b0
            b0_metadata = dwi_b0.out.metadata
        }

        if ( params.merge_repetitions ) {
            cat_dwi_repetitions_wkf(dwi_channel, meta_channel, "_dwi")
            dwi_channel = cat_dwi_repetitions_wkf.out.dwi
            meta_channel = cat_dwi_repetitions_wkf.out.metadata

            if ( params.has_reverse ) {
                cat_rev_repetitions_wkf(rev_channel, meta_channel, "_rev")
                rev_channel = cat_rev_repetitions_wkf.out.dwi
                meta_channel = meta_channel.join(cat_rev_repetitions_wkf.out.metadata)
            }

            t1_channel = merge_repetitions(t1_channel, false)
            average_t1(t1_channel.join(t1_channel.map{ [it[0], it[1][0].simpleName] }), "preprocess")
            t1_channel = average_t1.out.image
            if ( params.masked_t1 ) {
                t1_mask_channel = merge_repetitions(t1_mask_channel, false)
                merge_masks(t1_mask_channel.join(t1_mask_channel.map{ [it[0], it[1][0].simpleName] }), "preprocess")
                t1_mask_channel = merge_masks.out.mask
            }

            average_b0(b0_channel.join(b0_channel.map{ [it[0], it[1].simpleName] }), "preprocess")
            b0_channel = average_b0.out.image
        }

        if ( params.do_bet_mask ) {
            dwi_mask_channel = bet_mask(b0_channel, "preprocess")

            if (params.has_reverse) {
                rev_bet_mask(b0_channel, "preprocess")
                bet_merge_mask(dwi_mask_channel.join(rev_bet_mask).map { [it[0], it.subList(1, it.size()), "bet_mask"] }, "preprocess")
                dwi_mask_channel = bet_merge_mask.out.mask
            }
        }
        else if ( params.masked_t1 && params.t1mask2dwi_registration ) {
            if ( params.topup_correction ) {
                if ( !is_data(dwi_mask_channel) )
                    dwi_mask_channel = bet_mask_topup(b0_channel, "preprocess")

                apply_mask_to_b0_for_reg(b0_channel.join(dwi_mask_channel).map{ it + [""] }, "preprocess")
                topup_b0_channel = apply_mask_to_b0_for_reg.out.image

                scil_compute_dti_fa(topup_corrected_dwi.join(dwi_mask_channel), "preprocess", "preprocess")

                apply_mask_to_t1_for_reg(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess")
                topup_t1_channel = apply_mask_to_t1_for_reg.out.image

                topup_t1_mask_channel = dilate_t1_mask(t1_mask_channel, 3, "preprocess")
                topup_b0_mask_channel = dilate_b0_mask(dwi_mask_channel, 3, "preprocess")

                topup_mask_registration_wkf(
                    merge_channels_non_blocking(topup_b0_channel, scil_compute_dti_fa.out.fa),
                    topup_t1_channel.map{ [it[0], [it[1]]] },
                    t1_mask_channel,
                    topup_b0_mask_channel.join(topup_t1_mask_channel).map{ [it[0], [it[1], it[2]]] },
                    null,
                    meta_channel.map{ [it[0], it[1], ""] },
                    params.t1_mask_to_topup_b0_registration_config
                )

                topup_convert_datatype(topup_mask_registration_wkf.out.image, "int8", "preprocess")
                dwi_mask_channel = topup_convert_datatype.out.image
            }
            else {
                mask_registration_wkf(
                    dwi_channel,
                    t1_channel,
                    t1_mask_channel,
                    null,
                    meta_channel.map { [it[0], it[1]] }
                )

                convert_datatype(mask_registration_wkf.out.mask, "int8", "preprocess")
                dwi_mask_channel = convert_datatype.out.image

                if (params.has_reverse) {
                    rev_mask_registration_wkf(
                        rev_channel,
                        t1_channel,
                        t1_mask_channel,
                        null,
                        meta_channel.map{ [it[0], it[2]] }
                    )

                    rev_convert_datatype(rev_mask_registration_wkf.out.mask, "int8", "preprocess")
                    rev_mask_channel = uniformize_naming(rev_convert_datatype.out.image, "rev_mask", "false", "false")
                    reg_merge_mask(dwi_mask_channel.join(rev_mask_channel).map { [it[0], it.subList(1, it.size()), "bet_mask"] }, "preprocess")
                    dwi_mask_channel = reg_merge_mask.out.mask
                }
            }
        }
        else if ( params.masked_t1 ) {
            dwi_mask_channel = t1_mask_channel
        }

        if ( params.eddy_correction ) {
            eddy_wkf(dwi_channel, dwi_mask_channel, topup2eddy_channel, b0_channel, rev_channel, meta_channel)

            dwi_channel = eddy_wkf.out.dwi.join(eddy_wkf.out.bval).join(eddy_wkf.out.bvec)
            meta_channel = eddy_wkf.out.metadata
        }

        if ( params.intensity_normalization ) {
            n4_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, b0_channel, dwi_mask_channel, meta_channel, params.dwi_n4_normalization_config)
            dwi_channel = replace_dwi_file(dwi_channel, n4_denoise_wkf.out.image)
            meta_channel = n4_denoise_wkf.out.metadata
        }

        if ( !params.masked_t1 ) {
            t1_mask_registration_wkf(
                t1_channel.map{ [it[0], [it[1]]] },
                dwi_b0.out.b0.map{ [it[0], [it[1]]] },
                dwi_mask_channel,
                null,
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.b02t1_mask_registration_config
            )

            t1_mask_convert_datatype(t1_mask_registration_wkf.out.image, "int8", "preprocess")
            t1_mask_channel = t1_mask_convert_datatype.out.image
        }

        t1_preprocess_wkf(t1_channel.map{ it.subList(0, 2) }, t1_mask_channel)
        t1_channel = t1_preprocess_wkf.out.t1
        t1_mask_channel = t1_preprocess_wkf.out.mask

        if ( params.resample_data ) {
            scilpy_resample_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).join(meta_channel), "preprocess", "lin")
            if ( params.msmt_odf && params.seg_on_t1 ) {
                seg_to_resample = seg_channel.filter{ !it[1].isEmpty() }
                scilpy_resample_wm(seg_to_resample.map{ [it[0], it[1][0]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn")
                scilpy_resample_gm(seg_to_resample.map{ [it[0], it[1][1]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn")
                scilpy_resample_csf(seg_to_resample.map{ [it[0], it[1][2]] }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "nn")
                seg_channel = scilpy_resample_wm.out.image.join(scilpy_resample_gm.out.image).join(scilpy_resample_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(seg_channel.filter{ it[1].isEmpty() })
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
            )
            t1_channel = t1_registration_wkf.out.t1
            ants_transform_base_t1(t1_mask_channel.join(t1_channel).join(t1_registration_wkf.out.transform_base.map{ [it[0], it[2]] }).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
            ants_transform_base_dwi(t1_mask_channel.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
            t1_mask_channel = ants_transform_base_t1.out.image
            dwi_mask_channel = uniformize_naming(ants_transform_base_dwi.out.image, "dwi_mask", "false", "true")

            if ( params.seg_on_t1 ) {
                seg_to_register = seg_channel.filter{ !it[1].isEmpty() }
                ants_transform_base_wm(seg_to_register.map{ [it[0], it[1][0]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                ants_transform_base_gm(seg_to_register.map{ [it[0], it[1][1]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                ants_transform_base_csf(seg_to_register.map{ [it[0], it[1][2]] }.join(t1_registration_wkf.out.transform_base).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                seg_channel = ants_transform_base_wm.out.image.join(ants_transform_base_gm.out.image).join(ants_transform_base_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(seg_channel.filter{ it[1].isEmpty() })
            }

            if ( params.register_syn_t12b0 ) {
                ants_transform_syn_t1(t1_mask_channel.join(t1_channel).join(t1_registration_wkf.out.transform_syn.map{ [it[0], it[2]] }).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                ants_transform_syn_dwi(t1_mask_channel.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                t1_mask_channel = ants_transform_syn_t1.out.image
                dwi_mask_channel = uniformize_naming(ants_transform_syn_dwi.out.image, "dwi_mask", "false", "true")

                if ( params.seg_on_t1 ) {
                    seg_to_register = seg_channel.filter{ !it[1].isEmpty() }
                    ants_transform_syn_wm(seg_to_register.map{ [it[0], it[1][0]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                    ants_transform_syn_gm(seg_to_register.map{ [it[0], it[1][1]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                    ants_transform_syn_csf(seg_to_register.map{ [it[0], it[1][2]] }.join(t1_registration_wkf.out.transform_syn).map{ it + ["", ""] }, "preprocess", params.ants_transform_base_config)
                    seg_channel = ants_transform_syn_wm.out.image.join(ants_transform_syn_gm.out.image).join(ants_transform_syn_csf.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(seg_channel.filter{ it[1].isEmpty() })
                }
            }
        }

        crop_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).map{ it + [""] }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess")
        dwi_bbox_channel = crop_dwi.out.bbox
        fit_bounding_box(t1_channel.join(dwi_channel.map{ it.subList(0, 2) }).join(dwi_bbox_channel), "preprocess")
        dwi_bbox_channel = fit_bounding_box.out.bbox
        crop_t1(t1_channel.join(t1_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess")

        if ( params.seg_on_t1 ) {
            seg_to_crop = seg_channel.filter{ !it[1].isEmpty() }
            crop_wm(seg_to_crop.map{ [it[0], it[1][0]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess")
            crop_gm(seg_to_crop.map{ [it[0], it[1][1]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess")
            crop_csf(seg_to_crop.map{ [it[0], it[1][2]] }.join(dwi_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess")
            convert_wm_segmentation(crop_wm.out.image, "int8", "preprocess")
            convert_gm_segmentation(crop_gm.out.image, "int8", "preprocess")
            convert_csf_segmentation(crop_csf.out.image, "int8", "preprocess")
            seg_channel = convert_wm_segmentation.out.image.join(convert_gm_segmentation.out.image).join(convert_csf_segmentation.out.image).map{ [it[0], it.subList(1, it.size())] }.mix(seg_channel.filter{ it[1].isEmpty() })
        }

        dwi_channel = replace_dwi_file(dwi_channel, crop_dwi.out.image)
        meta_channel = crop_dwi.out.metadata
        dwi_mask_channel = crop_dwi.out.mask
        t1_channel = crop_t1.out.image
        t1_mask_channel = crop_dwi.out.mask

        if ( params.generate_segmentation ) {
            empty_segmentations = seg_channel.filter{ it[1].isEmpty() }.map{ [it[0]] }
            segment_nmt(empty_segmentations.join(t1_channel), empty_segmentations.join(t1_mask_channel))
            seg_channel = seg_channel.filter{ !it[1].isEmpty() }.mix(segment_nmt.out.masks.map{ [it[0], it.subList(1, it.size()).reverse()] })
        }

        dwi_channel = uniformize_naming(dwi_channel.map{ it.subList(0, 4) }, "dwi_preprocessed", "false", "false")
        meta_channel = uniformize_naming(meta_channel, "dwi_preprocessed_metadata", "false", "false")
        dwi_mask_channel = uniformize_naming(dwi_mask_channel, "mask_preprocessed", "false", "false")
    emit:
        t1 = t1_channel
        dwi = dwi_channel
        mask = dwi_mask_channel
        seg = seg_channel
        metadata = meta_channel
}

workflow t1_preprocess_wkf {
    take:
        t1_channel
        mask_channel
    main:
        if ( params.denoise_t1 ) {
            if ( params.nlmeans_t1 ) {
                nlmeans_denoise(t1_channel, "preprocess")
                t1_channel = nlmeans_denoise.out.image
            }
            else {
                ants_gaussian_denoise(t1_channel, "preprocess")
                t1_channel = ants_gaussian_denoise.out.image
            }
        }

        if ( params.intensity_normalization ) {
            n4_denoise_wkf(t1_channel, null, null, null, params.t1_n4_normalization_config)
            t1_channel = n4_denoise_wkf.out.image
        }

        if ( params.resample_data ) {
            scilpy_resample_t1(t1_channel.join(mask_channel).map{ it + [""] }, "preprocess", "lin")
            t1_channel = scilpy_resample_t1.out.image
            mask_channel = scilpy_resample_t1.out.mask
        }
    emit:
        t1 = t1_channel
        mask = mask_channel
}
