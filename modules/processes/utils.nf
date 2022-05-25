#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.add_odd_dimension = false
params.bet_f = 0.5
params.duplicates_merge_method = "mean"

include { remove_alg_suffixes; add_suffix } from '../functions.nf'

process apply_mask {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_masked") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(img), path(mask), file(metadata)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${img.simpleName}__masked.nii.gz"), emit: image
        tuple val(sid), path("${img.simpleName}__masked_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey apply_mask --in $img --mask $mask --out ${img.simpleName}__masked.nii.gz
        """
}

process bet_mask {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_bet_mask") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(img)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${img.simpleName}_bet_mask.nii.gz")
    script:
        """
        fslmaths $img -Tmean mean_image.nii.gz
        bet mean_image.nii.gz "${img.simpleName}_bet.nii.gz" -m -R -f $params.bet_f
        """
}

process cat_datasets {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(imgs), file(bval), file(bvec), file(metadatas)
        val(prefix)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${sid}_${prefix}__concatenated.nii.gz"), emit: image
        tuple val(sid), path("${sid}_${prefix}__concatenated.bval"), optional: true, emit: bval
        tuple val(sid), path("${sid}_${prefix}__concatenated.bvec"), optional: true, emit: bvec
        tuple val(sid), path("${sid}_${prefix}__concatenated_metadata.*"), optional: true, emit: metadata
    script:
        def args = "--in ${imgs.join(',')}"

        if ( bval.size() > 0 )
            args += " --bvals ${bval.join(',')}"
        if ( bvec.size() > 0 )
            args += " --bvecs ${bvec.join(',')}"

        """
        magic-monkey concatenate $args --out ${sid}_${prefix}__concatenated --config $config
        """
}

process split_image {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), path(img), path(metadata)
        val(split_axis)
        val(caller_name)
    output:
        tuple val(sid), path("${img.simpleName}_splitted_ax${split_axis}_[0-9]*.nii.gz"), emit: images
        tuple val(sid), path("${img.simpleName}_splitted_ax${split_axis}_*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey split --image $img --prefix "${img.simpleName}_splitted" --axis $split_axis
        """
}

process join_images {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/$caller_name", saveAs: { f -> f.contains("metadata") ? null : f }, mode: params.publish_mode

    input:
        tuple val(sid), val(prefix), path(imgs), path(metadatas)
        val(split_axis)
        val(caller_name)
    output:
        tuple val(sid), path("${sid}__joined_ax${split_axis}.nii.gz"), emit: image
        tuple val(sid), path("${sid}__joined_ax${split_axis}_*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey split --image ${sid}__joined_ax${split_axis}.nii.gz --prefix $prefix --axis $split_axis --inverse
        """
}

process apply_topup {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwis), path(bvals), path(bvecs), path(revs), path(topup_params), val(topup_prefix), path(topup_files), path(metadata)
        val(caller_name)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${sid}_dwi__topup_corrected_*.nii.gz"), path("${sid}_dwi__topup_corrected_*.bval"), path("${sid}_dwi__topup_corrected_*.bvec"), emit: dwi
        tuple val(sid), path("${sid}_dwi__topup_corrected_*_metadata.*"), optional: true, emit: metadata
    script:
        """
        magic-monkey apply_topup --dwi ${dwis.join(",")} --bvals ${bvals.join(",")} --bvecs ${bvecs.join(",")} --rev ${revs.join(",")} --acqp $topup_params --topup $topup_prefix --out ${sid}_dwi__topup_corrected
        """
}

process tournier2descoteaux_odf {
    label params.conservative_resources ? "res_conservative_cpu" : "res_max_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}/fodf", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(odfs)
        val(caller_name)
    output:
        tuple val(sid), path("${odfs.simpleName}_desc07_odf.nii.gz"), emit: odfs
    script:
        """
        scil_convert_sh_basis.py $odfs ${odfs.simpleName}_desc07_odf.nii.gz tournier07
        """
}

process convert_float_to_integer {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : publish_suffix ? "${sid}_${publish_suffix}.nii.gz" : remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image)
        val(datatype)
        val(caller_name)
        val(publish)
        val(publish_suffix)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${image.simpleName}__uint8.nii.gz"), emit: image
    script:
        """
        scil_image_math.py round $image ${image.simpleName}__${datatype}.nii.gz -f --data_type $datatype
        """
}

process replicate_image {
    label "res_single_cpu"

    input:
        tuple val(sid), path(img), path(ref_img)
        val(idx_to_rep)
    output:
        tuple val(sid), path("${img.simpleName}__replicated.nii.gz"), emit: image
    script:
        def args = ""
        if ( "$idx_to_rep" )
            args += "--idx $idx_to_rep"
        """
        magic-monkey replicate --in $img --ref $ref_img --out ${img.simpleName}__replicated.nii.gz $args
        """
}

process check_dwi_conformity {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode
    
    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(metadata)
        val(error_strategy)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__checked.nii.gz"), path("${dwi.simpleName}__checked.bval"), path("${dwi.simpleName}__checked.bvec"), emit: dwi
        tuple val(sid), path("${dwi.simpleName}__checked_metadata.*"), emit: metadata, optional: true
    script:
        """
        magic-monkey check --in $dwi --bvals $bval --bvecs $bvec --strat $error_strategy --out ${dwi.simpleName}__checked
        """
}

process pvf_to_mask {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(wm_pvf), path(gm_pvf), path(csf_pvf), path(brain_mask)
        val(caller_name)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${sid}_wm_mask.nii.gz"), emit: wm_mask
        tuple val(sid), path("${sid}_gm_mask.nii.gz"), emit: gm_mask
        tuple val(sid), path("${sid}_csf_mask.nii.gz"), emit: csf_mask
        tuple val(sid), path("${sid}_safe_wm_mask.nii.gz"), emit: safe_wm_mask
    script:
        """
        scil_image_math.py round $wm_pvf ${sid}_wm_mask.nii.gz --data_type uint8
        scil_image_math.py round $gm_pvf ${sid}_gm_mask.nii.gz --data_type uint8
        scil_image_math.py round $csf_pvf ${sid}_csf_mask.nii.gz --data_type uint8

        scil_image_math.py lower_threshold_eq $csf_pvf 0.001 csf_map.nii.gz
        scil_image_math.py dilation csf_map.nii.gz 1 csf_map.nii.gz -f --data_type uint8
        scil_image_math.py difference ${sid}_wm_mask.nii.gz csf_map.nii.gz ${sid}_safe_wm_mask.nii.gz
        scil_image_math.py difference ${sid}_safe_wm_mask.nii.gz ${sid}_gm_mask.nii.gz ${sid}_safe_wm_mask.nii.gz -f
        scil_image_math.py intersection ${sid}_safe_wm_mask.nii.gz $brain_mask ${sid}_safe_wm_mask.nii.gz -f
        """
}

process crop_image {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${["${params.output_root}/${sid}", additional_publish_path].findAll({ it }).join("/")}", saveAs: { f -> f.contains("${mask.simpleName}") ? ("$publish_mask" == "true") ? mask_prefix ? "${sid}_${mask_prefix}.nii.gz" : remove_alg_suffixes(f) : null : f.contains("cropped.nii.gz") ? remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), file(mask), file(bounding_box), file(metadata)
        val(caller_name)
        val(publish_mask)
        val(mask_prefix)
        val(additional_publish_path)
    output:
        tuple val(sid), path("${image.simpleName}__cropped.nii.gz"), emit: image
        tuple val(sid), path("${image.simpleName}__bbox.pkl"), emit: bbox, optional: true
        tuple val(sid), path("${mask.simpleName}__cropped.nii.gz"), emit: mask, optional: true
        tuple val(sid), path("${image.simpleName}__cropped_metadata.py"), emit: metadata, optional: true
    script:
        def args = ""
        def img = "$image"
        def before_script = []
        def after_script = []

        if ( !bounding_box.empty() ) {
            args += "--input_bbox $bounding_box"
            after_script += ["magic-monkey fit2box --in ${image.simpleName}__cropped.nii.gz --out ${image.simpleName}__cropped.nii.gz --pbox $bounding_box"]
        }
        else
            args += "--output_bbox ${image.simpleName}__bbox.pkl"

        if ( !mask.empty() ) {
            before_script = "magic-monkey apply_mask --in $image --mask $mask --out masked_image.nii.gz"
            img = "masked_image.nii.gz"
            mask_script = "magic-monkey fit2box --in $mask --out ${mask.simpleName}__cropped.nii.gz"
            img_script = "magic-monkey fit2box --in $image --out ${image.simpleName}__cropped.nii.gz"
            if ( !bounding_box.empty() ) {
                mask_script += " --pbox $bounding_box"
                img_script += " --pbox $bounding_box"
            }
            else {
                mask_script += " --pbox ${image.simpleName}__bbox.pkl"
                img_script += " --pbox ${image.simpleName}__bbox.pkl"
            }
            after_script += [img_script]
            after_script += [mask_script]
            after_script += ["scil_image_math.py round ${mask.simpleName}__cropped.nii.gz ${mask.simpleName}__cropped.nii.gz --data_type uint8 -f"]
        }

        if ( metadata instanceof nextflow.util.BlankSeparatedList ? !metadata.isEmpty() : !metadata.empty() )
            after_script += ["magic-monkey metadata --in ${image.getSimpleName()}__cropped.nii.gz --update_affine --metadata $metadata"]

        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        $before_script
        scil_crop_volume.py $img ${image.simpleName}__cropped.nii.gz $args
        ${after_script.join('\n')}
        if [ "\$(mrinfo -datatype $image)" != "\$(mrinfo -datatype ${image.simpleName}__cropped.nii.gz)" ]
        then
            mrconvert -force -datatype "\$(mrinfo -datatype $image)" ${image.simpleName}__cropped.nii.gz ${image.simpleName}__cropped.nii.gz
        fi
        """
}

process fit_bounding_box {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("cropped.nii.gz") ? remove_alg_suffixes(f) : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(image), path(reference), path(bounding_box)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}__bbox.pkl"), emit: bbox, optional: true
    script:
        """
        magic-monkey fitbox --in $image --ref $reference --pbox $bounding_box --out ${image.simpleName}__bbox
        """
}

process average {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(images), val(base_name)
        val(caller_name)
    output:
        tuple val(sid), path("${base_name}__averaged.nii.gz"), emit: image
    script:
        """
        magic-monkey concatenate --in ${images.join(",")} --out cat_images --ts
        fslmaths cat_images.nii.gz -Tmean ${base_name}__averaged.nii.gz
        """
}

process merge_masks {
    input:
        tuple val(sid), path(masks), val(base_name)
        val(caller_name)
    output:
        tuple val(sid), path("${base_name}__merged.nii.gz"), emit: mask
    script:
        """
        magic-monkey concatenate --in ${masks.join(",")} --out cat_images --ts
        fslmaths cat_images.nii.gz -Tmax ${base_name}__merged.nii.gz
        """
}

process timeseries_mean {
    input:
        tuple val(sid), path(image)
        val(caller_name)
    output:
        tuple val(sid), path("${image.simpleName}__mean.nii.gz"), emit: image
    script:
        """
        fslmaths $image -Tmean ${image.simpleName}__mean.nii.gz
        """
}

process extract_shells {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__extracted_shells.nii.gz"), path("${dwi.simpleName}__extracted_shells.bval"), path("${dwi.simpleName}__extracted_shells.bvec"), emit: dwi
    script:
        """
        magic-monkey shells --in $dwi --bvals $bval --bvecs $bvec --out ${dwi.simpleName}__extracted_shells --config $config
        """
}

process dilate_mask {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(mask)
        val(dilation_factor)
        val(caller_name)
    output:
        tuple val(sid), path("${mask.simpleName}__dilated.nii.gz")
    script:
        """
        scil_image_math.py dilation $mask $dilation_factor ${mask.simpleName}__dilated.nii.gz --data_type uint8
        """
}

process segmentation_to_binary {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(segmentation)
        val(caller_name)
    output:
        tuple val(sid), path("${segmentation.simpleName}_wm.nii.gz"), emit: wm_seg
        tuple val(sid), path("${segmentation.simpleName}_gm.nii.gz"), emit: gm_seg
        tuple val(sid), path("${segmentation.simpleName}_csf.nii.gz"), emit: csf_seg
        tuple val(sid), path("${segmentation.simpleName}_dgm.nii.gz"), emit: dgm_seg
        tuple val(sid), path("${segmentation.simpleName}_all_gm.nii.gz"), emit: all_gm_seg
    script:
        """
        magic-monkey seg2mask --in $segmentation --values 1,2,3,4 --labels csf,gm,dgm,wm --out ${segmentation.simpleName}
        scil_image_math.py addition ${segmentation.simpleName}_gm.nii.gz ${segmentation.simpleName}_dgm.nii.gz ${segmentation.simpleName}_all_gm.nii.gz --data_type uint8 -f
        """
}

process prepend_sid {
    label "res_single_cpu"

    input:
        tuple val(sid), path(file)
    output:
        tuple val(sid), path("${sid}_${file.getName()}")
    script:
        """
        ln -s $file ${sid}_${file.getName()}
        """
}

process generate_b0_bval {
    label "res_single_cpu"

    input:
        tuple val(sid), path(b0_image)
        val(with_bvec)
    output:
        tuple val(sid), path("${b0_image.simpleName}.bval"), emit: bval
        tuple val(sid), path("${b0_image.simpleName}.bvec"), optional: true, emit: bvec
    script:
        if (with_bvec == "true") {
            """
            echo "0" >> ${b0_image.simpleName}.bval
            echo "0\n0\n0" >> ${b0_image.simpleName}.bvec
            """
        }
        else {
            """
            echo "0" >> ${b0_image.simpleName}.bval
            """
        }
}

process check_odd_dimensions {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("${reverse.simpleName}") ? null : f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(reverse), file(rval), file(rvec), file(mask), file(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__even_dims.nii.gz"), path("${dwi.simpleName}__even_dims.bval"), path("${dwi.simpleName}__even_dims.bvec"), emit: dwi
        tuple val(sid), path("${reverse.simpleName}__even_dims.nii.gz"), optional: true, emit: rev
        tuple val(sid), path("${reverse.simpleName}__even_dims.bval"), path("${reverse.simpleName}__even_dims.bvec"), optional: true, emit: rev_bval_bvec
        tuple val(sid), path("${mask.simpleName}__even_dims.nii.gz"), optional: true, emit: mask
        tuple val(sid), path("*__even_dims_metadata.*"), optional: true, emit: metadata
    script:
        def args = "--strat ${params.add_odd_dimension ? "add" : "sub"}"
        def after_script = ""
        def assoc = []
        if ( !reverse.empty() ) assoc += ["$reverse"]
        if ( !mask.empty() ) assoc += ["$mask"]
        args += " --assoc ${assoc.join(",")}"
        if ( !rval.empty() ) after_script += "cp $rval ${reverse.simpleName}__even_dims.bval\n"
        if ( !rvec.empty() ) after_script += "cp $rvec ${reverse.simpleName}__even_dims.bvec\n"
        """
        magic-monkey even_dimensions --in $dwi --suffix __even_dims $args
        cp $bval ${dwi.simpleName}__even_dims.bval
        cp $bvec ${dwi.simpleName}__even_dims.bvec
        $after_script
        """
}

process check_for_duplicates {
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("${reverse.simpleName}") ? null : f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.nii.gz"), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.bval"), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.bvec"), emit: dwi
        tuple val(sid), path("*__${params.duplicates_merge_method}_duplicates_metadata.*"), optional: true, emit: metadata   
    script:
        """
        magic-monkey duplicates --in $dwi --bvals $bval --bvecs $bvec --merge $params.duplicates_merge_method --out ${dwi.simpleName}__${params.duplicates_merge_method}_duplicates
        """
}