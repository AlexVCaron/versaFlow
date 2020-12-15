#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Preprocess workflow parameters
params.gaussian_noise_correction = true
params.rev_is_b0 = true
params.gibbs_ringing_correction = true
params.t1mask2dwi_registration = true
params.masked_t1 = true
params.masked_dwi = false
params.topup_correction = true
params.eddy_correction = true
params.eddy_pre_bet_mask = false
params.post_eddy_registration = true
params.intensity_normalization = true
params.resample_data = true
params.register_repetitions = true
params.register_t12b0_denoised = true
params.register_syn_t12b0 = false

// T1 preprocess workflow parameters
params.denoise_t1 = true
params.nlmeans_t1 = true
// params.intensity_normalization = true
// params.resample_data = true

params.config.workflow.preprocess.t12b0_base_registration = file("$projectDir/.config/.workflow/t12b0_base_registration.py")
params.config.workflow.preprocess.t12b0_syn_registration = file("$projectDir/.config/.workflow/t12b0_syn_registration.py")
params.config.workflow.preprocess.t12b0mask_registration = file("$projectDir/.config/.workflow/t12b0_mask_registration.py")
params.config.register.ants_transform = file("$projectDir/.config/ants_transform.py")
params.config.register.ants_motion = file("$projectDir/.config/ants_motion.py")
params.config.workflow.preprocess.b0_mean = file("$projectDir/.config/extract_b0_mean.py")
params.config.workflow.preprocess.b0_batch_mean = file("$projectDir/.config/extract_b0_batch_mean.py")
params.config.workflow.preprocess.first_b0 = file("$projectDir/.config/extract_first_b0.py")
params.config.denoise.n4_denoise = file("$projectDir/.config/n4_denoise.py")
params.config.workflow.preprocess.n4_denoise_t1 = file("$projectDir/.config/.workflow/n4_denoise_on_t1.py")

include { merge_channels_non_blocking; map_optional; opt_channel; replace_dwi_file; uniformize_naming; merge_repetitions } from '../modules/functions.nf'
include { extract_b0 as dwi_b0; extract_b0 as extract_b0_motion; extract_b0 as dwi_b0_for_t1_reg } from '../modules/processes/preprocess.nf'
include { ants_correct_motion } from '../modules/processes/register.nf'
include { scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { ants_transform } from '../modules/processes/register.nf'
include { convert_datatype; convert_datatype as t1_mask_convert_datatype; bet_mask; crop_image as crop_dwi; crop_image as crop_t1; fit_bounding_box; average; merge_masks } from '../modules/processes/utils.nf'
include { gibbs_removal as dwi_gibbs_removal; gibbs_removal as rev_gibbs_removal; nlmeans_denoise; ants_gaussian_denoise } from '../modules/processes/denoise.nf'
include { scilpy_resample as scilpy_resample_t1; scilpy_resample_on_ref as scilpy_resample_t1_mask; scilpy_resample as scilpy_resample_dwi; scilpy_resample_on_ref as scilpy_resample_mask } from '../modules/processes/upsample.nf'
include { dwi_denoise_wkf; dwi_denoise_wkf as rev_denoise_wkf; squash_wkf; registration_wkf as mask_registration_wkf; registration_wkf as t1_mask_registration_wkf; registration_wkf as t1_base_registration_wkf; registration_wkf as t1_syn_registration_wkf; topup_wkf; eddy_wkf; apply_topup_wkf; n4_denoise_wkf } from "../modules/workflows/preprocess.nf"
include { cat_dwi_repetitions_wkf; cat_dwi_repetitions_wkf as cat_rev_repetitions_wkf; register_dwi_repetitions_wkf; register_t1_repetitions_wkf } from '../modules/workflows/repetitions.nf'

workflow preprocess_wkf {
    take:
        dwi_channel
        rev_channel
        t1_channel
        meta_channel
        rev_meta_channel
    main:
        dwi_mask_channel = map_optional(dwi_channel, 4)
        t1_mask_channel = map_optional(t1_channel, 2)
        t1_channel = t1_channel.map{ it.subList(0, 2) }
        dwi_channel = dwi_channel.map{ it.subList(0, 4) }
        topup2eddy_channel = opt_channel()
        topup2eddy_b0_channel = opt_channel()

        if ( params.gaussian_noise_correction ) {
            dwi_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, dwi_mask_channel, meta_channel)
            dwi_channel = replace_dwi_file(dwi_channel, dwi_denoise_wkf.out.image)
            meta_channel = dwi_denoise_wkf.out.metadata

            if ( !params.rev_is_b0 ) {
                rev_denoise_wkf(rev_channel.map{ it.subList(0, 2) }, dwi_mask_channel, rev_meta_channel)
                rev_channel = replace_dwi_file(rev_channel, rev_denoise_wkf.out.image)
                rev_meta_channel = rev_denoise_wkf.out.metadata
            }
        }

        if ( params.gibbs_ringing_correction ) {
            dwi_gibbs_removal(dwi_channel.map{ it.subList(0, 2) }.join(meta_channel), "preprocess")
            dwi_channel = replace_dwi_file(dwi_channel, dwi_gibbs_removal.out.image)
            meta_channel = dwi_gibbs_removal.out.metadata

            if ( !params.rev_is_b0 ) {
                rev_gibbs_removal(rev_channel.map{ it.subList(0, 2) }.join(rev_meta_channel), "preprocess")
                rev_channel = replace_dwi_file(rev_channel, rev_gibbs_removal.out.image)
                rev_meta_channel = rev_gibbs_removal.out.metadata
            }
        }

        if ( params.register_repetitions ) {
            register_dwi_repetitions_wkf(
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
            dwi_channel = register_dwi_repetitions_wkf.out.dwi
            rev_channel = register_dwi_repetitions_wkf.out.rev
            meta_channel = register_dwi_repetitions_wkf.out.metadata
            rev_meta_channel  = register_dwi_repetitions_wkf.out.rev_metadata
        }

        squash_wkf(dwi_channel, rev_channel, meta_channel.join(rev_meta_channel))
        dwi_channel = squash_wkf.out.dwi
        rev_channel = squash_wkf.out.rev
        meta_channel = squash_wkf.out.metadata

        if ( params.topup_correction ) {
            topup_wkf(dwi_channel, rev_channel, meta_channel)

            topup2eddy_channel = topup_wkf.out.param.join(topup_wkf.out.prefix).join(topup_wkf.out.topup.map{ [it[0], it.subList(1, it.size())] })

            if ( !params.eddy_correction ) {
                dwi_channel = uniformize_naming(dwi_channel, "dwi_to_topup", "false")
                rev_channel = uniformize_naming(rev_channel, "dwi_to_topup_rev", "false")
                meta_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(0..<it[1].size()).step(2)] }, "dwi_to_topup_metadata", "false")
                rev_meta_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0]] + it[1][(1..<it[1].size()).step(2)] }, "dwi_to_topup_rev_metadata", "false")
                apply_topup_wkf(dwi_channel, rev_channel, topup2eddy_channel, meta_channel.join(rev_meta_channel).map{ [it[0], it.subList(1, it.size())] })
                dwi_channel = uniformize_naming(apply_topup_wkf.out.dwi, "topup_corrected", "false")
                meta_channel = uniformize_naming(apply_topup_wkf.out.metadata, "topup_corrected_metadata", "false")
            }
        }

        if ( params.merge_repetitions ) {
            if ( (!params.topup_correction || params.eddy_correction) ) {
                cat_dwi_repetitions_wkf(dwi_channel, meta_channel, "_dwi")
                cat_rev_repetitions_wkf(rev_channel, meta_channel, "_rev")
                dwi_channel = cat_dwi_repetitions_wkf.out.dwi
                rev_channel = cat_rev_repetitions_wkf.out.dwi
                meta_channel = cat_dwi_repetitions_wkf.out.metadata.join(cat_rev_repetitions_wkf.out.metadata)
            }
            t1_channel = merge_repetitions(t1_channel, false)
            average(t1_channel.join(t1_channel.map{ [it[0], it[1][0].simpleName] }), "preprocess")
            t1_channel = average.out.image
            if ( params.masked_t1 ) {
                t1_mask_channel = merge_repetitions(t1_mask_channel, false)
                merge_masks(t1_mask_channel.join(t1_mask_channel.map{ [it[0], it[1][0].simpleName] }), "preprocess")
                t1_mask_channel = merge_masks.out.mask
            }
        }

        dwi_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "", "preprocess", params.config.workflow.preprocess.b0_mean)
        b0_channel = dwi_b0.out.b0
        b0_metadata = dwi_b0.out.metadata

        if ( params.masked_t1 && params.t1mask2dwi_registration ) {
            mask_registration_wkf(
                b0_channel.map{ [it[0], [it[1]]] },
                t1_channel.map{ [it[0], [it[1]]] },
                t1_mask_channel,
                null,
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.config.workflow.preprocess.t12b0mask_registration
            )

            convert_datatype(mask_registration_wkf.out.image, "int8", "preprocess")
            dwi_mask_channel = convert_datatype.out.image
        }
        else if ( params.masked_t1 ) {
            dwi_mask_channel = t1_mask_channel
        }

        if ( (params.eddy_correction && params.eddy_pre_bet_mask) || (!( params.masked_dwi ) &&  params.masked_t1 && !params.t1mask2dwi_registration ) ) {
            dwi_mask_channel = bet_mask(b0_channel, "preprocess")
        }

        if ( params.eddy_correction ) {
            eddy_wkf(dwi_channel, dwi_mask_channel, topup2eddy_channel, b0_channel, rev_channel, meta_channel)

            dwi_channel = eddy_wkf.out.dwi.join(eddy_wkf.out.bval).join(eddy_wkf.out.bvec)
            meta_channel = eddy_wkf.out.metadata

            if ( params.post_eddy_registration ) {
                extract_b0_motion(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "_eddy", "preprocess", params.config.workflow.preprocess.b0_mean)
                ants_correct_motion(dwi_channel.map{ [it[0], [it[1]]] }.join(extract_b0_motion.out.b0.map{ [it[0], [it[1]]] }).join(meta_channel), "preprocess", params.config.register.ants_motion)
                dwi_channel = replace_dwi_file(dwi_channel, ants_correct_motion.out.image)
                meta_channel = ants_correct_motion.out.metadata
            }
        }

        if ( params.intensity_normalization ) {
            n4_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, b0_channel, dwi_mask_channel, meta_channel, params.config.denoise.n4_denoise)
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
                params.config.workflow.preprocess.t12b0mask_registration
            )

            t1_mask_convert_datatype(t1_mask_registration_wkf.out.image, "int8", "preprocess")
            t1_mask_channel = t1_mask_convert_datatype.out.image
        }

        t1_preprocess_wkf(t1_channel.map{ it.subList(0, 2) }, t1_channel.map{ [it[0], it[2]] })
        t1_channel = t1_preprocess_wkf.out.t1

        if ( params.register_t12b0_denoised ) {
            dwi_b0_for_t1_reg(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "", "preprocess", params.config.workflow.preprocess.b0_mean)
            scil_compute_dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess")
            b0_metadata = dwi_b0_for_t1_reg.out.metadata
            t1_base_registration_wkf(
                dwi_b0_for_t1_reg.out.b0.map{ [it[0], [it[1]]] },
                t1_channel.map{ [it[0], [it[1]]] },
                null,
                dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] },
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.config.workflow.preprocess.t12b0_base_registration
            )
            ants_transform(t1_mask_channel.join(t1_base_registration_wkf.out.transform).map{ it + [""] }, "preprocess", params.config.register.ants_transform)
            t1_mask_channel = ants_transform.out.image
            if ( params.register_syn_t12b0 ) {
                t1_syn_registration_wkf(
                    merge_channels_non_blocking(dwi_b0_for_t1_reg.out.b0, scil_compute_dti_fa.out.fa),
                    t1_base_registration_wkf.out.image.map{ [it[0], [it[1]]] },
                    null,
                    params.register_syn_t12b0_with_mask ? dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] } : null,
                    b0_metadata.map{ it.subList(0, 2) + [""] },
                    params.config.workflow.preprocess.t12b0_syn_registration
                )
                t1_channel = t1_syn_registration_wkf.out.image
            }
        }

        crop_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).map{ it + [""] }.join(meta_channel.map{ [it[0], it.subList(1, it.size())] }), "preprocess")
        dwi_bbox_channel = crop_dwi.out.bbox
        fit_bounding_box(t1_channel.join(dwi_channel.map{ it.subList(0, 2) }).join(dwi_bbox_channel), "preprocess")
        dwi_bbox_channel = fit_bounding_box.out.bbox
        crop_t1(t1_channel.join(t1_mask_channel).join(dwi_bbox_channel).map{ it + [""] }, "preprocess")
        dwi_channel = replace_dwi_file(dwi_channel, crop_dwi.out.image)
        dwi_mask_channel = crop_dwi.out.mask
        t1_channel = crop_t1.out.image
        t1_mask_channel = crop_t1.out.mask

        dwi_channel = uniformize_naming(dwi_channel.map{ it.subList(0, 4) }, "dwi_preprocessed", "false")
        meta_channel = uniformize_naming(meta_channel, "dwi_preprocessed_metadata", "false")
        dwi_mask_channel = uniformize_naming(dwi_mask_channel, "mask_preprocessed", "false")
    emit:
        t1 = t1_channel
        dwi = dwi_channel
        mask = dwi_mask_channel
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
            n4_denoise_wkf(t1_channel, null, null, null, params.config.workflow.preprocess.n4_denoise_t1)
            t1_channel = n4_denoise_wkf.out.image
        }

        if ( params.resample_data ) {
            scilpy_resample_t1(t1_channel.map{ it + ["", ""] }, "preprocess", "lin")
            t1_channel = scilpy_resample_t1.out.image
        }
    emit:
        t1 = t1_channel
}