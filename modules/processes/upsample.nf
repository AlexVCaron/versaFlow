#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.resampling_resolution = 1
params.force_resampling_sequential = false

include { get_size_in_gb; remove_alg_suffixes } from '../functions.nf'

process scilpy_resample {
    label params.force_resampling_sequential ? "res_max_cpu" : "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> f.contains("metadata") ? null : f.contains("${mask.simpleName}") ? ("$publish_mask" == "true") ? mask_prefix ? "${sid}_${mask_prefix}.nii.gz" : remove_alg_suffixes(f) : null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask), file(metadata)
        val(caller_name)
        val(interpolation)
        val(publish_mask)
        val(mask_prefix)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${image.getSimpleName()}__resampled.nii.gz"), emit: image
        tuple val(sid), path("${mask.simpleName}__resampled.nii.gz"), optional: true, emit: mask
        tuple val(sid), path("${image.getSimpleName()}__resampled_metadata.py"), optional: true, emit: metadata
    script:
        def after_script = ""
        if ( !mask.empty() ) {
            after_script += "scil_resample_volume.py $mask mask_resampled.nii.gz --ref ${image.simpleName}__resampled.nii.gz --interp nn --enforce_dimensions\n"
            after_script += "scil_image_math.py round mask_resampled.nii.gz ${mask.simpleName}__resampled.nii.gz --data_type uint8 -f\n"
        }
        if ( !metadata.empty() )
            after_script += "magic-monkey metadata --in ${image.getSimpleName()}__resampled.nii.gz --update_affine --metadata $metadata\n"
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_resample_volume.py $image resampled.nii.gz --voxel_size $params.resampling_resolution --interp $interpolation
        fslmaths resampled.nii.gz -thr 0 ${image.simpleName}__resampled.nii.gz
        if [ "\$(mrinfo -datatype $image)" != "\$(mrinfo -datatype ${image.simpleName}__resampled.nii.gz)" ]
        then
            mrconvert -force -datatype "\$(mrinfo -datatype $image)" ${image.simpleName}__resampled.nii.gz ${image.simpleName}__resampled.nii.gz
        fi
        $after_script
        """
}
