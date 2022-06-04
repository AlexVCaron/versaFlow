#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.register_syn_t12b0 = true
params.tissue_segmentation_root = "$moduleDir/../../.data/maccaca_mulatta/tissue_segmentation"

params.ants_transform_base_config = file("$moduleDir/../../.config/ants_transform_base_config.py")
params.t1_registration_extract_b0_config = file("$moduleDir/../../.config/extract_mean_b0_base_config.py")
params.t1_to_template_registration_config = file("$moduleDir/../../.config/t1_to_template_registration_config.py")
params.b0_to_template_registration_config = file("$moduleDir/../../.config/b0_to_template_registration_config.py")

include {
    extract_b0;
    compute_powder_average
} from '../processes/preprocess.nf'
include { scil_compute_dti_fa } from '../processes/measure.nf'
include {
    registration_wkf as t1_to_template_registration_wkf;
    registration_wkf as b0_to_template_registration_wkf
} from "./preprocess.nf"
include {
    ants_transform as ants_transform_t1_to_b0;
    ants_transform as ants_transform_t1_mask_to_b0
} from '../processes/register.nf'
include { merge_channels_non_blocking; is_data } from '../functions.nf'
include {
    apply_mask as apply_mask_to_b0_for_reg;
    apply_mask as apply_mask_to_t1_for_reg;
    prepend_sid as prepend_sid_template;
    bet_mask
} from '../processes/utils.nf'

workflow t12b0_registration {
    take:
        dwi_channel
        t1_channel
        t1_mask_channel
        dwi_mask_channel
        meta_channel
        publish_mask
        publish_t1
    main:
        extract_b0(dwi_channel.map{ it.subList(0, 3) + [""] }, "preprocess", "false", params.t1_registration_extract_b0_config)

        if ( !is_data(dwi_mask_channel) )
            dwi_mask_channel = bet_mask(extract_b0.out.b0, "preprocess", "false")

        apply_mask_to_b0_for_reg(extract_b0.out.b0.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        apply_mask_to_t1_for_reg(t1_channel.join(t1_mask_channel).map{ it + [""] }, "preprocess", "false")
        compute_powder_average(dwi_channel.map{ it.subList(0, 3) }.join(dwi_mask_channel).map{ it + [""] }, "preprocess", "false")
        scil_compute_dti_fa(dwi_channel.join(dwi_mask_channel), "preprocess", "preprocess", false)

        b0_metadata = extract_b0.out.metadata
        reg_b0_channel = apply_mask_to_b0_for_reg.out.image
        reg_t1_channel = apply_mask_to_t1_for_reg.out.image

        target_channel = apply_mask_to_b0_for_reg.out.image.join(compute_powder_average.out.image).join(scil_compute_dti_fa.out.fa)
        moving_channel = apply_mask_to_t1_for_reg.out.image.map{ [it[0], [it[1]]] }
        template_channel = prepend_sid_template(t1_channel.map{ [it[0], file("${params.tissue_segmentation_root}/tissue_segmentation_t1.nii.gz")] }).map{ [it[0], [it[1]]] }

        t1_to_template_registration_wkf(
            moving_channel,
            template_channel,
            t1_mask_channel,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            params.t1_to_template_registration_config,
            params.ants_transform_base_config
        )

        b0_to_template_registration_wkf(
            target_channel,
            template_channel,
            null,
            null,
            null,
            null,
            "",
            false,
            "",
            "",
            params.b0_to_template_registration_config,
            params.ants_transform_base_config
        )

        ants_transform_t1_to_b0(
            t1_to_template_registration_wkf.out.registration.join(b0_to_template_registration_wkf.out.inverse_transform),
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_base_config
        )
        ants_transform_t1_mask_to_b0(
            t1_to_template_registration_wkf.out.image.join(b0_to_template_registration_wkf.out.inverse_transform),
            "preprocess",
            "","$publish_mask", "mask",
            params.ants_transform_base_config
        )

        trans = t1_to_template_registration_wkf.out.transform.join(b0_to_template_registration_wkf.out.inverse_transform)

    emit:
        t1 = ants_transform_t1_to_b0.out.image
        mask = ants_transform_t1_mask_to_b0.out.image
        transform = trans.map{ [it[0], it[1] + it[2]] }
}