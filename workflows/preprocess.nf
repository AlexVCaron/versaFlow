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
params.config.workflow.preprocess.b0_mean = file("$projectDir/.config/extract_b0_mean.py")
params.config.workflow.preprocess.b0_batch_mean = file("$projectDir/.config/extract_b0_batch_mean.py")
params.config.workflow.preprocess.first_b0 = file("$projectDir/.config/extract_first_b0.py")
params.config.workflow.preprocess.n4_denoise_t1 = file("$projectDir/.config/.workflow/n4_denoise_on_t1.py")

include { map_optional; opt_channel; replace_dwi_file; uniformize_naming } from '../modules/functions.nf'
include { extract_b0 as dwi_b0; extract_b0 as extract_b0_motion; extract_b0 as dwi_b0_for_t1_reg } from '../modules/processes/preprocess.nf'
include { ants_correct_motion } from '../modules/processes/register.nf'
include { scil_compute_dti_fa } from '../modules/processes/measure.nf'
include { ants_transform } from '../modules/processes/register.nf'
include { convert_datatype; convert_datatype as t1_mask_convert_datatype; bet_mask; crop_image as crop_dwi; crop_image as crop_t1; fit_bounding_box } from '../modules/processes/utils.nf'
include { gibbs_removal as dwi_gibbs_removal; gibbs_removal as rev_gibbs_removal; nlmeans_denoise; ants_gaussian_denoise } from '../modules/processes/denoise.nf'
include { scilpy_resample as scilpy_resample_t1; scilpy_resample_on_ref as scilpy_resample_t1_mask; scilpy_resample as scilpy_resample_dwi; scilpy_resample_on_ref as scilpy_resample_mask } from '../modules/processes/upsample.nf'
include { dwi_denoise_wkf; dwi_denoise_wkf as rev_denoise_wkf; squash_wkf; registration_wkf as mask_registration_wkf; registration_wkf as t1_mask_registration_wkf; registration_wkf as t1_base_registration_wkf; registration_wkf as t1_syn_registration_wkf; topup_wkf; eddy_wkf; apply_topup_wkf; n4_denoise_wkf } from "../modules/workflows/preprocess.nf"

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

        squash_wkf( dwi_channel, rev_channel, meta_channel.join(rev_meta_channel) )

        dwi_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "", "preprocess", "")
        b0_channel = dwi_b0.out.b0
        b0_metadata = dwi_b0.out.metadata

        if ( params.masked_t1 && params.t1mask2dwi_registration ) {
            mask_registration_wkf(
                b0_channel.groupTuple(),
                t1_channel.map{ [it[0], it[1]] }.groupTuple(),
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

        if ( params.topup_correction ) {
            topup_wkf(squash_wkf.out.dwi, squash_wkf.out.rev, squash_wkf.out.metadata)

            topup2eddy_channel = topup_wkf.out.param.join(topup_wkf.out.prefix).join(topup_wkf.out.topup.map{ [it[0], it.subList(1, it.size())] })
            b0_channel = topup_wkf.out.b0

            if ( !params.eddy_correction )
                meta_channel = topup_wkf.out.metadata.map{ [it[0], it[2]] }
        }

        if ( (params.eddy_correction && params.eddy_pre_bet_mask) || (!( params.masked_dwi ) &&  params.masked_t1 && !params.t1mask2dwi_registration ) ) {
            dwi_mask_channel = bet_mask(b0_channel, "preprocess")
        }

        if ( params.eddy_correction ) {
            eddy_wkf(squash_wkf.out.dwi, dwi_mask_channel, topup2eddy_channel, b0_channel, squash_wkf.out.rev, squash_wkf.out.metadata)

            dwi_channel = eddy_wkf.out.dwi.join(eddy_wkf.out.bval).join(eddy_wkf.out.bvec)
            meta_channel = eddy_wkf.out.metadata

            if ( params.post_eddy_registration ) {
                extract_b0_motion(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "_eddy", "preprocess", "")
                ants_correct_motion(dwi_channel.map{ it.subList(0, 2) }.groupTuple().join(extract_b0_motion.out.b0.groupTuple()).join(meta_channel), "preprocess", "")
                dwi_channel = replace_dwi_file(dwi_channel, ants_correct_motion.out.image)
                meta_channel = ants_correct_motion.out.metadata
            }
        }
        else if ( params.topup_correction ) {
            dwi_channel = uniformize_naming(squash_wkf.out.dwi, "dwi_to_topup", "false")
            rev_channel = uniformize_naming(squash_wkf.out.rev, "dwi_to_topup_rev", "false")
            meta_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0], it[1][0]] }, "dwi_to_topup_metadata", "false")
            rev_meta_channel = uniformize_naming(topup_wkf.out.in_metadata_w_topup.map{ [it[0], it[1][1]] }, "dwi_to_topup_rev_metadata", "false")
            apply_topup_wkf(dwi_channel, rev_channel, topup2eddy_channel, meta_channel.join(rev_meta_channel).map{ [it[0], it.subList(1, it.size())] })
            dwi_channel = replace_dwi_file(dwi_channel, apply_topup_wkf.out.dwi)
            dwi_channel = uniformize_naming(dwi_channel, "topup_corrected", "false")
            meta_channel = uniformize_naming(apply_topup_wkf.out.metadata, "topup_corrected_metadata", "false")
        }

        if ( params.intensity_normalization ) {
            n4_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, b0_channel, dwi_mask_channel, meta_channel, "")
            dwi_channel = replace_dwi_file(dwi_channel, n4_denoise_wkf.out.image)
            meta_channel = n4_denoise_wkf.out.metadata
        }

        if ( !params.masked_t1 ) {
            t1_mask_registration_wkf(
                t1_channel.map{ [it[0], it[1]] }.groupTuple(),
                dwi_b0.out.b0.groupTuple(),
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
            dwi_b0_for_t1_reg(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "", "preprocess", "")
            scil_compute_dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess")
            b0_metadata = dwi_b0_for_t1_reg.out.metadata
            t1_base_registration_wkf(
                dwi_b0_for_t1_reg.out.b0.groupTuple(),
                t1_channel.map{ [it[0], it[1]] }.groupTuple(),
                null,
                dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] },
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.config.workflow.preprocess.t12b0_base_registration
            )
            ants_transform(t1_mask_channel.join(t1_base_registration_wkf.out.transform).map{ it + [""] }, "preprocess")
            t1_mask_channel = ants_transform.out.image
            if ( params.register_syn_t12b0 ) {
                t1_syn_registration_wkf(
                    dwi_b0_for_t1_reg.out.b0.join(scil_compute_dti_fa.out.fa).groupTuple().map{ [it[0], it.subList(1, it.size()).inject([]){ c, t -> c + t }] },
                    t1_base_registration_wkf.out.image.groupTuple(),
                    null,
                    params.register_syn_t12b0_with_mask ? dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] } : null,
                    b0_metadata.map{ it.subList(0, 2) + [""] },
                    params.config.workflow.preprocess.t12b0_syn_registration
                )
                t1_channel = t1_syn_registration_wkf.out.image
            }
        }

        crop_dwi(dwi_channel.map{ it.subList(0, 2) }.join(dwi_mask_channel).map{ it + [""] }.join(meta_channel), "preprocess")
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
