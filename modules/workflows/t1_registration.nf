#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.register_syn_t12b0 = true

params.config.register.ants_transform = file("$projectDir/.config/ants_transform.py")
params.config.workflow.preprocess.b0_mean = file("$projectDir/.config/extract_b0_mean.py")
params.config.workflow.preprocess.t12b0_base_registration = file("$projectDir/.config/.workflow/t12b0_base_registration.py")
params.config.workflow.preprocess.t12b0_syn_registration = file("$projectDir/.config/.workflow/t12b0_syn_registration.py")

include { extract_b0 } from '../processes/preprocess.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include { registration_wkf as t1_base_registration_wkf; registration_wkf as t1_syn_registration_wkf } from "./preprocess.nf"
include { ants_transform as ants_transform_base; ants_transform as ants_transform_syn } from '../processes/register.nf'
include { merge_channels_non_blocking } from '../functions.nf'

workflow t12b0_registration {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        dwi_mask_channel
        meta_channel
    main:
        extract_b0(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "preprocess", params.config.workflow.preprocess.b0_mean)
        b0_metadata = extract_b0.out.metadata

        if ( dwi_mask_channel ) {
            mask_channel = dwi_mask_channel.join(t1_mask_channel).map{ [it[0], [it[1], it[2]]] }
        }
        else {
            mask_channel = t1_mask_channel.map{ [it[0], [it[1]]] }
        }

        t1_base_registration_wkf(
            extract_b0.out.b0.map{ [it[0], [it[1]]] },
            t1_channel.map{ [it[0], [it[1]]] },
            null,
            mask_channel,
            null,
            b0_metadata.map{ it.subList(0, 2) + [""] },
            params.config.workflow.preprocess.t12b0_base_registration
        )

        ants_transform_base(t1_mask_channel.join(t1_base_registration_wkf.out.transform).map { it + [""] }, "preprocess", params.config.register.ants_transform)
        t1_mask_channel = ants_transform_base.out.image

        if ( params.register_syn_t12b0 ) {
            if ( dwi_mask_channel ) {
                in_fa = dwi_channel.join(dwi_mask_channel)
            }
            else {
                in_fa = dwi_channel.map{ it + [""] }
            }

            if ( params.register_syn_t12b0_with_mask ) {
                syn_mask = mask_channel
            }
            else {
                syn_mask = null
            }

            scil_compute_dti_fa(in_fa, "preprocess", "preprocess")
            t1_syn_registration_wkf(
                merge_channels_non_blocking(extract_b0.out.b0, scil_compute_dti_fa.out.fa),
                t1_base_registration_wkf.out.image.map{ [it[0], [it[1]]] },
                null,
                syn_mask,
                null,
                b0_metadata.map{ it.subList(0, 2) + [""] },
                params.config.workflow.preprocess.t12b0_syn_registration
            )
            t1_channel = t1_syn_registration_wkf.out.image
            t1_syn_registration_wkf.out.transform.view()
            ants_transform_syn(t1_mask_channel.join(t1_syn_registration_wkf.out.transform).map { it + [""] }, "preprocess", params.config.register.ants_transform)
            t1_mask_channel = ants_transform_syn.out.image
        }
    emit:
        t1 = t1_channel
        mask = t1_mask_channel
        transform = params.register_syn_t12b0 ? t1_base_registration_wkf.out.transform.join(t1_syn_registration_wkf.out.transform) : t1_base_registration_wkf.out.transform
}