#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.force_resampling_resolution = false
params.resampling_min_resolution = false
params.resampling_subdivision = 2
params.force_resampling_sequential = false

include { remove_alg_suffixes } from '../functions.nf'

process scilpy_resample {
    label "FAST"
    label params.force_resampling_sequential ? "res_max_cpu" : "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> f.contains("metadata") ? null : f.contains("${mask.simpleName}") ? ("$publish_mask" == "true") ? mask_prefix ? "${sid}_${mask_prefix}.nii.gz" : remove_alg_suffixes(f) : null : remove_alg_suffixes(f) }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(image), file(mask), file(metadata), file(reference)
        val(caller_name)
        val(interpolation)
        val(antstype)
        val(antsinterp)
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
            after_script += "scil_resample_volume.py $mask mask_resampled.nii.gz --ref ${image.simpleName}__resampled.nii.gz --interp nn\n"
            after_script += "scil_image_math.py floor mask_resampled.nii.gz ${mask.simpleName}__resampled.nii.gz --data_type uint8 -f\n"
        }
        if ( !metadata.empty() )
            after_script += "mrhardi metadata --in ${image.getSimpleName()}__resampled.nii.gz --update_affine --metadata $metadata\n"
        if ( !reference.empty() ) {
            after_script += "antsApplyTransforms -v -d 3 -e $antstype -n $antsinterp -r $reference -i ${image.simpleName}__resampled.nii.gz -o ${image.simpleName}__resampled.nii.gz -t identity\n"
        }
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_resample_volume.py $image resampled.nii.gz --voxel_size $params.force_resampling_resolution --interp $interpolation
        fslmaths resampled.nii.gz -thr 0 ${image.simpleName}__resampled.nii.gz
        if [ "\$(mrinfo -datatype $image)" != "\$(mrinfo -datatype ${image.simpleName}__resampled.nii.gz)" ]
        then
            mrconvert -force -datatype "\$(mrinfo -datatype $image)" ${image.simpleName}__resampled.nii.gz ${image.simpleName}__resampled.nii.gz
        fi
        $after_script
        """
}

process scilpy_resample_to_reference {
    label "FAST"
    label params.force_resampling_sequential ? "res_max_cpu" : "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : f.contains("${mask.simpleName}") ? ("$publish_mask" == "true") ? mask_prefix ? "${sid}_${mask_prefix}.nii.gz" : remove_alg_suffixes(f) : null : remove_alg_suffixes(f) : null }, mode: params.publish_mode, overwrite: true

    input:
        tuple val(sid), path(image), path(reference), file(mask), file(metadata)
        val(caller_name)
        val(interpolation)
        val(publish)
        val(publish_mask)
        val(mask_prefix)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${image.getSimpleName()}__resampled.nii.gz"), emit: image
        tuple val(sid), path("${mask.getSimpleName()}__resampled.nii.gz"), optional: true, emit: mask
        tuple val(sid), path("${image.getSimpleName()}__resampled_metadata.py"), optional: true, emit: metadata
    script:
        def after_script = ""
        if ( !mask.empty() ) {
            after_script += "scil_resample_volume.py $mask mask_resampled.nii.gz --ref $reference --interp nn --enforce_dimensions\n"
            after_script += "scil_image_math.py floor mask_resampled.nii.gz ${mask.getSimpleName()}__resampled.nii.gz --data_type uint8 -f\n"
        }
        if ( !metadata.empty() )
            after_script += "mrhardi metadata --in ${image.getSimpleName()}__resampled.nii.gz --update_affine --metadata $metadata\n"
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_resample_volume.py $image resampled.nii.gz --ref $reference --interp $interpolation --enforce_dimensions
        fslmaths resampled.nii.gz -thr 0 ${image.simpleName}__resampled.nii.gz
        if [ "\$(mrinfo -datatype $image)" != "\$(mrinfo -datatype ${image.simpleName}__resampled.nii.gz)" ]
        then
            mrconvert -force -datatype "\$(mrinfo -datatype $image)" ${image.simpleName}__resampled.nii.gz ${image.simpleName}__resampled.nii.gz
        fi
        $after_script
        """
}

process resampling_reference {
    label "FAST"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process.replaceAll(":", "/")}", mode: "$params.publish_all_mode", enabled: params.publish_all, overwrite: true

    input:
        tuple val(sid), path(images)
        val(caller_name)
        val(max_subdivisions)
        val(min_voxel_size)
        val(force_voxel_size)
    output:
        tuple val(sid), path("${sid}_resampling_reference.nii.gz"), emit: reference
    script:
        def args = ""
        if ( force_voxel_size )
            args += " --force_resolution $force_voxel_size"
        if ( min_voxel_size )
            args += " --min_voxel_size $min_voxel_size"
        """
        mrhardi resampling_reference \
            --subdiv $max_subdivisions \
            --in ${images.join(",")} \
            --out ${sid}_resampling_reference.nii.gz \
            $args
        """
}
