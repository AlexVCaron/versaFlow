#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.register_syn_t12b0 = true

params.ants_transform_base_config = file("$projectDir/.config/ants_transform_base_config.py")
params.t1_registration_extract_b0_config = file("$projectDir/.config/extract_mean_b0_base_config.py")
params.t1_registration_base_registration_config = file("$projectDir/.config/t1_registration_base_registration_config.py")
params.t1_registration_syn_registration_config = file("$projectDir/.config/t1_registration_syn_registration_config.py")

include { extract_b0 } from '../processes/preprocess.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include { registration_wkf as t1_base_registration_wkf; registration_wkf as t1_syn_registration_wkf } from "./preprocess.nf"
include {
    ants_transform as ants_mask_transform_base; ants_transform as ants_mask_transform_syn;
    ants_transform as ants_t1_transform_base; ants_transform as ants_t1_transform_syn
} from '../processes/register.nf'
include { merge_channels_non_blocking; is_data } from '../functions.nf'
include {
    apply_mask as apply_mask_to_b0_for_reg; apply_mask as apply_mask_to_t1_for_reg;
    dilate_mask as dilate_t1_mask; dilate_mask as dilate_b0_mask;
    bet_mask
} from '../processes/utils.nf'

workflow t12b0_registration {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        dwi_mask_channel
        meta_channel
    main:
        extract_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "preprocess", params.t1_registration_extract_b0_config)
        b0_metadata = extract_b0.out.metadata

        if ( !is_data(dwi_mask_channel) )
            dwi_mask_channel = bet_mask(extract_b0.out.b0, "preprocess")

        apply_mask_to_b0_for_reg(extract_b0.out.b0.join(dwi_mask_channel).map{ it + [""] }, "preprocess")
        reg_b0_channel = apply_mask_to_b0_for_reg.out.image

        scil_compute_dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess")

        apply_mask_to_t1_for_reg(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess")
        reg_t1_channel = apply_mask_to_t1_for_reg.out.image

        reg_t1_mask_channel = dilate_t1_mask(t1_mask_channel, 3, "preprocess")
        reg_b0_mask_channel = dilate_b0_mask(dwi_mask_channel, 3, "preprocess")

        t1_base_registration_wkf(
            merge_channels_non_blocking(reg_b0_channel, scil_compute_dti_fa.out.fa),
            reg_t1_channel.map{ [it[0], [it[1]]] },
            null,
            reg_b0_mask_channel.join(reg_t1_mask_channel).map{ [it[0], [it[1], it[2]]] },
            null,
            b0_metadata.map{ it.subList(0, 2) + [""] },
            "",
            params.t1_registration_base_registration_config
        )

        ants_mask_transform_base(t1_mask_channel.join(t1_base_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", "", params.ants_transform_base_config)
        ants_t1_transform_base(t1_channel.join(t1_base_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", "", params.ants_transform_base_config)

        t1_mask_channel = ants_mask_transform_base.out.image
        t1_channel = ants_t1_transform_base.out.image

        if ( params.register_syn_t12b0 ) {
            if ( params.register_syn_t12b0_with_mask )
                syn_mask = dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] }
            else
                syn_mask = null

            t1_syn_registration_wkf(
                merge_channels_non_blocking(reg_b0_channel, scil_compute_dti_fa.out.fa),
                t1_base_registration_wkf.out.image.map{ [it[0], [it[1]]] },
                null,
                syn_mask,
                null,
                b0_metadata.map{ it.subList(0, 2) + [""] },
                "",
                params.t1_registration_syn_registration_config
            )

            ants_mask_transform_syn(t1_mask_channel.join(t1_syn_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", "", params.ants_transform_base_config)
            ants_t1_transform_syn(t1_channel.join(t1_syn_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", "", params.ants_transform_base_config)

            t1_mask_channel = ants_mask_transform_syn.out.image
            t1_channel = ants_t1_transform_syn.out.image
        }
    emit:
        t1 = t1_channel
        mask = t1_mask_channel
        transform_base = t1_base_registration_wkf.out.transform
        transform_syn = params.register_syn_t12b0 ? t1_syn_registration_wkf.out.transform : null
}