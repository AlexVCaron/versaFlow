#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.versaflow_inputs = false
params.versaflow_outputs = false
params.output_dir = "reports"

params.qc_screenshot_orientation = "axial"

include {
    get_id;
    load_dataset
} from '../workflows/io.nf'

include {
    qc_screenshot_parameters_wkf;
    dwi_qc_wkf;
    masks_qc_wkf;
    tissue_segmentation_qc_wkf;
    dti_qc_wkf;
    fodf_qc_wkf
} from '../modules/workflows/qc.nf'


workflow {
    if ( !params.versaflow_inputs )
        error "You must supply a versaflow input root using --versaflow_inputs"
    input_root = file(params.versaflow_inputs)

    if ( !params.versaflow_outputs )
        error "You must supply a versaflow output root using --versaflow_outputs"
    output_root = file(params.versaflow_outputs)

    input_dataloader = load_dataset(input_root)
    output_dataloader = load_versaflow_outputs(output_root)

    parameters = qc_screenshot_parameters_wkf(output_dataloader.atlases)

    // dwi_qc_wkf(
    //     input_dataloader.dwi,
    //     output_dataloader.dwi,
    //     parameters.cc_bounding_box,
    //     params.qc_screenshot_orientation
    // )

    // masks_qc_wkf(
    //     output_dataloader.b0,
    //     output_dataloader.t1,
    //     output_dataloader.dwi_mask,
    //     output_dataloader.t1_mask
    // )

    tissue_segmentation_qc_wkf(
        output_dataloader.t1,
        output_dataloader.segmentation_mask,
        output_dataloader.segmentation_3t,
        output_dataloader.segmentation_safe_mask
    )

    // dti_qc_wkf(
    //     output_dataloader.t1,
    //     output_dataloader.dti_evecs,
    //     output_dataloader.dti_scalars,
    //     output_dataloader.segmentation_mask
    // )

    fodf_qc_wkf(
        output_dataloader.t1,
        output_dataloader.fodf_peaks,
        output_dataloader.fodf_scalars,
        output_dataloader.segmentation_mask
    )
}


def join_tissues(data_channel) {
    joined_channel = data_channel
        .map{ [it[0], it[1..-1]] }
        .transpose()
        .branch{
            wm: it[1] =~ /_wm_/
            gm: it[1] =~ /_gm_/
        }

    wm_stack_channel = joined_channel.wm
        .groupTuple()
        .map{ [it[0], it[1].sort()] }
    gm_stack_channel = joined_channel.gm
        .groupTuple()
        .map{ [it[0], it[1].sort()] }

    return wm_stack_channel
        .join(gm_stack_channel)
        .transpose()
        .groupTuple()
}


workflow load_versaflow_outputs {
    take:
        versaflow_output_root
    main:

        // Preprocessed DWI and T1 images and masks

        dwi_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_dwi.{nii.gz,bval,bvec}", size: 3, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }
            .map{ sid, bval, bvec, dwi -> [sid, dwi, bval, bvec]}

        b0_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_dwi_b0.nii.gz", size: 1, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }

        t1_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_t1.nii.gz", size: 1, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }

        t1_masked_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_t1_masked.nii.gz", size: 1, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }

        dwi_mask_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_dwi_mask.nii.gz", size: 1, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }

        t1_mask_channel = Channel.fromFilePairs("$versaflow_output_root/**/*_t1_mask.nii.gz", size: 1, maxDepth: 1, flat: true)
            { get_id(it.parent, versaflow_output_root) }

        // Segmentation images

        atlases_channel = Channel.fromFilePairs("$versaflow_output_root/**/atlases/*_atlas.nii.gz", size: -1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }
            .map{ [it[0], it[1..-1]] }

        segmentation_channel = Channel.fromFilePairs("$versaflow_output_root/**/segmentation/*_segmentation.nii.gz", size:1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        segmentation_pvf_channel = Channel.fromFilePairs("$versaflow_output_root/**/segmentation/*_pvf.nii.gz", size: 1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }
            .branch{
                pvf_3t: it[1].simpleName =~ /3t_(wm|gm|csf)_pvf/
                pvf_mt: true
            }

        segmentation_pvf_3t_channel = segmentation_pvf_channel.pvf_3t
            .groupTuple()
            .map{ [it[0]] + it[1] }
            .map{ sid, csf, gm, wm -> [sid, wm, gm, csf] }

        segmentation_pvf_mt_channel = segmentation_pvf_channel.pvf_mt
            .groupTuple()
            .map{ [it[0]] + it[1] }

        segmentation_all_mask_channel = Channel.fromFilePairs("$versaflow_output_root/**/segmentation/*_{wm,gm,csf}_mask.nii.gz", maxDepth: 2, size:-1, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }
            .map{ [it[0], it[1..-1]] }
            .transpose()
            .branch{
                safe_wm: it[1].simpleName =~ /safe_wm_mask/
                masks_3t: true
            }

        segmentation_safe_mask_channel = segmentation_all_mask_channel.safe_wm
        segmentation_mask_channel = segmentation_all_mask_channel.masks_3t
            .groupTuple()
            .map{ it.flatten() }
            .map{ sid, csf, gm, wm -> [sid, wm, gm, csf] }

        // Diffusion modeling images

        dti_scalars_channel = Channel.fromFilePairs("$versaflow_output_root/**/dti/*_dti_{ad,fa,ga,md,mode,norm,rd,residuals,rgb}.nii.gz", size: 9, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }
            .map{ sid, ad, fa, ga, md, mode, norm, rd, residuals, rgb -> [sid, md, ad, rd, fa, rgb, ga, residuals, norm, mode] }

        dti_evecs_channel = Channel.fromFilePairs("$versaflow_output_root/**/dti/*_dti_evecs_v{1,2,3}.nii.gz", size: 3, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        dti_tensors_channel = Channel.fromFilePairs("$versaflow_output_root/**/dti/*_dti_dti.nii.gz", size: 1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_scalars_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_fodf_metrics_{wm,gm}*{afd,afds,afdt,nufo,rgb}.nii.gz", size: -1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_response_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_response.txt", size: 1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_ventricle_max_mask_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_ventricles_mask.nii.gz", size: 1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_coefficients_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_fodf.nii.gz", size: 1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_peaks_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_fodf_metrics_{wm,gm}*peaks.nii.gz", size: -1, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_peaks_scalars_channel = Channel.fromFilePairs("$versaflow_output_root/**/fodf/*_fodf_metrics_{wm,gm}*peaks_{indices,values}.nii.gz", size: 2, maxDepth: 2, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        fodf_scalars_channel = join_tissues(fodf_scalars_channel)
        fodf_peaks_channel = fodf_peaks_channel.map{ [it[0], it[1..-1]] }
        fodf_peaks_scalars_channel = join_tissues(fodf_peaks_scalars_channel)

        // Transforms between spaces (DWI, T1, template)

        t1_to_dwi_transforms_channel = Channel.fromFilePairs("$versaflow_output_root/**/transforms/t1_to_b0_space/composite_transforms/*_image_transform*.nii.gz", size: 2, maxDepth: 4, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        t1_to_template_transforms_channel = Channel.fromFilePairs("$versaflow_output_root/**/transforms/t1_to_template_space/composite_transforms/*_image_transform*.nii.gz", size: 2, maxDepth: 4, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

        dwi_to_template_transforms_channel = Channel.fromFilePairs("$versaflow_output_root/**/transforms/b0_to_template_space/composite_transforms/*_image_transform*.nii.gz", size: 2, maxDepth: 4, flat: true)
            { get_id(it.parent.parent, versaflow_output_root) }

    emit:
        // Preprocessed DWI and T1 images and masks
        dwi = dwi_channel
        b0 = b0_channel
        t1 = t1_channel
        t1_masked = t1_masked_channel
        dwi_mask = dwi_mask_channel
        t1_mask = t1_mask_channel
        // Segmentation images
        atlases = atlases_channel
        segmentation = segmentation_channel
        segmentation_3t = segmentation_pvf_3t_channel
        segmentation_mt = segmentation_pvf_mt_channel
        segmentation_mask = segmentation_mask_channel
        segmentation_safe_mask = segmentation_safe_mask_channel
        // Diffusion modeling images
        dti_scalars = dti_scalars_channel
        dti_evecs = dti_evecs_channel
        dti_tensors = dti_tensors_channel
        fodf_scalars = fodf_scalars_channel
        fodf_response = fodf_response_channel
        fodf_ventricle_max_mask = fodf_ventricle_max_mask_channel
        fodf_coefficients = fodf_coefficients_channel
        fodf_peaks = fodf_peaks_channel
        fodf_peaks_scalars = fodf_peaks_scalars_channel
        // Transforms between spaces (DWI, T1, template)
        t1_to_dwi_transforms = t1_to_dwi_transforms_channel
        t1_to_template_transforms = t1_to_template_transforms_channel
        dwi_to_template_transforms = dwi_to_template_transforms_channel
}