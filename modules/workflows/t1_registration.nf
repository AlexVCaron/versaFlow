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
include { ants_transform as ants_mask_transform_base; ants_transform as ants_mask_transform_syn; ants_transform as ants_t1_transform_base; ants_transform as ants_t1_transform_syn } from '../processes/register.nf'
include { apply_mask } from '../processes/utils.nf'
include { merge_channels_non_blocking } from '../functions.nf'

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

        if ( dwi_mask_channel ) {
            mask_channel = dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] }
            in_fa = dwi_channel.join(dwi_mask_channel)
        }
        else {
            mask_channel = t1_mask_channel.map{ [it[0], [it[1]]] }
            in_fa = dwi_channel.map{ it + [""] }
        }

        scil_compute_dti_fa(in_fa, "preprocess", "preprocess")

        apply_mask(t1_channel.join(t1_mask_channel).map{ it + [""] })

        t1_base_registration_wkf(
            merge_channels_non_blocking(extract_b0.out.b0, scil_compute_dti_fa.out.fa),
            apply_mask.out.image.map{ [it[0], [it[1]]] },
            null,
            mask_channel,
            null,
            b0_metadata.map{ it.subList(0, 2) + [""] },
            params.t1_registration_base_registration_config
        )

        ants_mask_transform_base(t1_mask_channel.join(t1_base_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", params.ants_transform_base_config)
        ants_t1_transform_base(t1_channel.join(t1_base_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", params.ants_transform_base_config)
        t1_mask_channel = ants_mask_transform_base.out.image
        t1_channel = ants_t1_transform_base.out.image

        if ( params.register_syn_t12b0 ) {
            if ( params.register_syn_t12b0_with_mask ) {
                if ( dwi_mask_channel ) {
                    syn_mask = dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] }
                }
                else {
                    syn_mask = t1_mask_channel.map{ [it[0], [it[1]]] }
                }
            }
            else {
                syn_mask = null
            }

            t1_syn_registration_wkf(
                merge_channels_non_blocking(extract_b0.out.b0, scil_compute_dti_fa.out.fa),
                t1_base_registration_wkf.out.image.map{ [it[0], [it[1]]] },
                null,
                syn_mask,
                null,
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.t1_registration_syn_registration_config
            )

            ants_mask_transform_syn(t1_mask_channel.join(t1_syn_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", params.ants_transform_base_config)
            ants_t1_transform_syn(t1_channel.join(t1_syn_registration_wkf.out.transform).map { it + ["", ""] }, "preprocess", params.ants_transform_base_config)
            t1_mask_channel = ants_mask_transform_syn.out.image
            t1_channel = ants_t1_transform_syn.out.image
        }
    emit:
        t1 = t1_channel
        mask = t1_mask_channel
        transform_base = t1_base_registration_wkf.out.transform
        transform_syn = params.register_syn_t12b0 ? t1_syn_registration_wkf.out.transform : null
}