#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.config.utils.apply_mask = "$projectDir/.config/apply_mask.py"
params.config.utils.concatenate = "$projectDir/.config/cat.py"

include { get_size_in_gb } from '../functions.nf'

process apply_mask {
    memory { 2f * get_size_in_gb([img, mask]) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    beforeScript "cp $params.config.utils.apply_mask config.py"
    input:
        tuple val(sid), path(img), path(mask), path(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${img.simpleName}__masked.nii.gz"), emit: image
        tuple val(sid), path("${img.simpleName}__masked_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey apply_mask $img $mask ${img.simpleName}__masked.nii.gz --config config.py
        """
}

process bet_mask {
    memory { get_size_in_gb(img) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(img)
        val(caller_name)
    output:
        tuple val(sid), path("${img.simpleName}__bet_mask.nii.gz")
    script:
        """
        fslmaths $img -Tmean $img
        bet $img "${img.simpleName}__bet.nii.gz" -m -R -f $params.bet.f
        """
}

process cat_datasets {
    memory { 2f * get_size_in_gb(imgs) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    beforeScript "cp $params.config.utils.concatenate config.py"
    input:
        tuple val(sid), file(imgs), file(bval), file(bvec), file(metadatas)
        val(suffix)
        val(caller_name)
    output:
        tuple val(sid), file("${sid}__concatenated${suffix}.nii.gz"), emit: image
        tuple val(sid), file("${sid}__concatenated${suffix}.bval"), optional: true, emit: bval
        tuple val(sid), file("${sid}__concatenated${suffix}.bvec"), optional: true, emit: bvec
        tuple val(sid), file("${sid}__concatenated${suffix}_metadata.*"), optional: true, emit: metadata
    script:
        args = "--in ${imgs.join(',')}"

        if ( bval.size() > 0 )
            args += " --bvals ${bval.join(',')}"
        if ( bvec.size() > 0 )
            args += " --bvecs ${bvec.join(',')}"

        """
        magic-monkey concatenate $args --out ${sid}__concatenated${suffix} --config config.py
        """
}

process split_image {
    memory { get_size_in_gb(img) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(img), path(metadata)
        val(split_axis)
        val(caller_name)
    output:
        tuple val(sid), path("${img.simpleName}__splitted_ax${split_axis}_[0-9]*.nii.gz"), emit: images
        tuple val(sid), path("${img.simpleName}__splitted_ax${split_axis}_*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey split --image $img --prefix "${img.simpleName}__splitted" --axis $split_axis
        """
}

process join_images {
    memory { 2f * get_size_in_gb(imgs) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), val(prefix), path(imgs), path(metadatas)
        val(split_axis)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}__joined_ax${split_axis}.nii.gz"), emit: image
        tuple val(sid), file("${sid}__joined_ax${split_axis}_*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey split --image ${sid}__joined_ax${split_axis}.nii.gz --prefix $prefix --axis $split_axis --inverse
        """
}

process apply_topup {
    memory { 2f * (get_size_in_gb(dwis) + get_size_in_gb(revs)) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwis), path(bvals), path(bvecs), path(revs), path(topup_params), val(topup_prefix), path(topup_files), path(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}__topup_corrected*.nii.gz"), path("${sid}__topup_corrected*.bval"), path("${sid}__topup_corrected*.bvec"), emit: dwi
        tuple val(sid), file("${sid}__topup_corrected*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey apply_topup --dwi ${dwis.join(",")} --bvals ${bvals.join(",")} --bvecs ${bvecs.join(",")} --rev ${revs.join(",")} --acqp $topup_params --topup $topup_prefix --out ${sid}__topup_corrected
        """
}

process tournier2descoteaux_odf {
    memory { 2f * get_size_in_gb(odfs) }
    label params.conservative_resources ? "res_conservative" : "res_max_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(odfs)
        val(caller_name)
    output:
        tuple val(sid), path("${odfs.simpleName}__desc07_odf.nii.gz"), emit: odfs
    script:
        """
        scil_convert_sh_basis.py $odfs ${odfs.simpleName}__desc07_odf.nii.gz tournier07
        """
}

process convert_datatype {
    memory { 2f * get_size_in_gb(image) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(image)
        val(datatype)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}_dt_${datatype}.nii.gz"), emit: image
    script:
        """
        magic-monkey convert --in $image --out "${image.simpleName}_dt_${datatype}.nii.gz" --dt $datatype
        """
}

process replicate_image {
    memory { 2f * get_size_in_gb([img, ref_img]) }
    label "res_single_cpu"
    errorStrategy "finish"

    input:
        tuple val(sid), path(img), path(ref_img)
        val(idx_to_rep)
    output:
        tuple val(sid), path("${img.simpleName}__replicated.nii.gz"), emit: image
    script:
        args = ""
        if ( "$idx_to_rep" )
            args += "--idx $idx_to_rep"
        """
        magic-monkey replicate --in $img --ref $ref_img --out ${img.simpleName}__replicated.nii.gz $args
        """
}

process check_dwi_conformity {
    label "res_single_cpu"
    errorStrategy "finish"

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(metadata)
        val(error_strategy)
    output:
        tuple val(sid), path("${dwi.simpleName}_checked.nii.gz"), path("${dwi.simpleName}_checked.bval"), path("${dwi.simpleName}_checked.bvec"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}_checked_metadata.*"), emit: metadata, optional: true
    script:
        """
        magic-monkey check --in $dwi --bvals $bval --bvecs $bvec --strat $error_strategy --out ${dwi.simpleName}_checked
        """
}

process crop_image {
    memory { 2f * get_size_in_gb([image, mask]) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("cropped.nii.gz") ? f : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask), file(bounding_box), file(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}_cropped.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}_bbox.pkl"), emit: bbox, optional: true
        tuple val(sid), path("${mask.simpleName}_cropped.nii.gz"), emit: mask, optional: true
    script:
        args = ""
        after_script = []

        if ( !bounding_box.empty() ) {
            args += "--input_bbox $bounding_box"
            after_script += ["magic-monkey fit2box --in ${image.simpleName}_cropped.nii.gz --out ${image.simpleName}_cropped.nii.gz --pbox $bounding_box"]
        }
        else
            args += "--output_bbox ${image.simpleName}_bbox.pkl"

        if ( !mask.empty() ) {
            mask_script = "magic-monkey fit2box --in $mask --out ${mask.simpleName}_cropped.nii.gz"
            if ( !bounding_box.empty() )
                mask_script += " --pbox $bounding_box"
            else
                mask_script += " --pbox ${image.simpleName}_bbox.pkl"
            after_script += [mask_script]
            after_script += ["scil_image_math.py convert ${mask.simpleName}_cropped.nii.gz ${mask.simpleName}_cropped.nii.gz --data_type uint8 -f"]
        }

        if ( metadata instanceof nextflow.util.BlankSeparatedList ? !metadata.isEmpty() : !metadata.empty() )
            after_script += ["magic-monkey metadata --in ${image.getSimpleName()}_cropped.nii.gz --update_affine --metadata $metadata"]

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_crop_volume.py $image ${image.simpleName}_cropped.nii.gz $args
        if [ "\$(mrinfo -datatype $image)" != "\$(mrinfo -datatype ${image.simpleName}_cropped.nii.gz)" ]
        then
            mrconvert -force -datatype "\$(mrinfo -datatype $image)" ${image.simpleName}_cropped.nii.gz ${image.simpleName}_cropped.nii.gz
        fi
        ${after_script.join('\n')}
        """
}

process fit_bounding_box {
    memory { get_size_in_gb([image]) }
    label "res_single_cpu"
    errorStrategy "finish"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.process}_${task.index}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("cropped.nii.gz") ? f : null }, mode: params.publish_mode

    input:
        tuple val(sid), file(image), file(reference), file(bounding_box)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}_bbox.pkl"), emit: bbox, optional: true
    script:
    """
    magic-monkey fitbox --in $image --ref $reference --pbox $bounding_box --out ${image.simpleName}_bbox
    """
}

process average {
    input:
        tuple val(sid), path(images), val(base_name)
        val(caller_name)
    output:
        tuple val(sid), path("${base_name}_averaged.nii.gz"), emit: image
    script:
    """
    magic-monkey concatenate --in ${images.join(",")} --out cat_images --ts
    fslmaths cat_images.nii.gz -Tmean ${base_name}_averaged.nii.gz
    """
}

process merge_masks {
    input:
    tuple val(sid), path(masks), val(base_name)
    val(caller_name)
    output:
    tuple val(sid), path("${base_name}_merged.nii.gz"), emit: mask
    script:
    """
    magic-monkey concatenate --in ${masks.join(",")} --out cat_images --ts
    fslmaths cat_images.nii.gz -Tmax ${base_name}_merged.nii.gz
    """
}