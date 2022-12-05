#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.add_odd_dimension = false
params.b0_threshold = false
params.shell_threshold = false
params.bet_f = 0.5
params.min_pvf_threshold = 0.001
params.max_safe_csf_pvf_threshold = 0.01
params.max_safe_gm_pvf_threshold = 0.01
params.safe_csf_mask_dilation = 1
params.safe_gm_mask_dilation = 1
params.duplicates_merge_method = "mean"
params.validate_bvecs_fa_thr = 0.2

include { remove_alg_suffixes; add_suffix } from '../functions.nf'

process apply_mask {
    label "FAST"
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
        mrhardi apply_mask --in $img --mask $mask --out ${img.simpleName}__masked.nii.gz
        """
}

process bet_mask {
    label "BET"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> ("$publish" == "true") ? f.contains("metadata") ? null : add_suffix(remove_alg_suffixes(f), "_bet_mask") : null }, mode: params.publish_mode

    input:
        tuple val(sid), path(img)
        val(caller_name)
        val(publish)
    output:
        tuple val(sid), path("${img.simpleName}_bet_mask.nii.gz"), emit: mask
    script:
        """
        fslmaths $img -Tmean mean_image.nii.gz
        bet mean_image.nii.gz "${img.simpleName}_bet.nii.gz" -m -R -f $params.bet_f
        """
}

process cat_datasets {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(imgs), file(bval), file(bvec), file(metadatas)
        val(axis)
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
        def single_copy_and_exit = ""
        if ( bval.size() > 0 )
            args += " --bvals ${bval.join(',')}"
        if ( bvec.size() > 0 )
            args += " --bvecs ${bvec.join(',')}"
        if ( !("$axis" == "") )
            args += " --axis $axis"
        if ( imgs.size() == 1 ) {
            single_copy_and_exit = "cp $imgs ${sid}_${prefix}__concatenated.nii.gz\n"
            if ( bval.size() > 0 ) single_copy_and_exit += "cp $bval ${sid}_${prefix}__concatenated.bval\n"
            if ( bvec.size() > 0 ) single_copy_and_exit += "cp $bvec ${sid}_${prefix}__concatenated.bvec\n"
            if ( metadata.size() > 0 ) single_copy_and_exit += "cp $metadata ${sid}_${prefix}__concatenated_metadata.py\n"
            single_copy_and_exit += "exit 0\n"
        }
        """
        $single_copy_and_exit
        mrhardi concatenate $args --out ${sid}_${prefix}__concatenated --config $config
        """
}

process split_image {
    label "MEDIUM"
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
        if [[ \$(( $split_axis + 1 )) -gt \$(mrinfo -ndim ${img}) ]]
        then
            cp $img ${img.simpleName}_splitted_ax${split_axis}_0.nii.gz
            cp $metadata ${img.simpleName}_splitted_ax${split_axis}_0_metadata.py
        else
            mrhardi split --image $img --prefix "${img.simpleName}_splitted" --axis $split_axis
        fi
        """
}

process join_images {
    label "MEDIUM"
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
        mrhardi split --image ${sid}__joined_ax${split_axis}.nii.gz --prefix $prefix --axis $split_axis --inverse
        """
}

process apply_topup {
    label "FAST"
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
        mrhardi apply_topup --dwi ${dwis.join(",")} --bvals ${bvals.join(",")} --bvecs ${bvecs.join(",")} --rev ${revs.join(",")} --acqp $topup_params --topup $topup_prefix --out ${sid}_dwi__topup_corrected
        """
}

process tournier2descoteaux_odf {
    label "MEDIUM"
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
    label "LIGHTSPEED"
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
        scil_image_math.py floor $image ${image.simpleName}__${datatype}.nii.gz -f --data_type $datatype
        """
}

process replicate_image {
    label "FAST"
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
        mrhardi replicate --in $img --ref $ref_img --out ${img.simpleName}__replicated.nii.gz $args
        """
}

process check_dwi_conformity {
    label "FAST"
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
        def args = ""
        if (params.b0_threshold)
            args += " --ceil ${params.b0_threshold}"
        """
        mrhardi check --in $dwi --bvals $bval --bvecs $bvec --strat $error_strategy --out ${dwi.simpleName}__checked $args
        """
}

process pvf_to_mask {
    label "FAST"
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
        python3 - <<'END_SCRIPT'
        import nibabel as nib
        import numpy as np
        wm_pvf = nib.load("$wm_pvf").get_fdata()
        gm_pvf = nib.load("$gm_pvf").get_fdata()
        csf_pvf = nib.load("$csf_pvf")
        affine = csf_pvf.affine
        csf_pvf = csf_pvf.get_fdata()

        wm_background = wm_pvf < $params.min_pvf_threshold
        wm_pvf[wm_background] = 0.
        gm_background = gm_pvf < $params.min_pvf_threshold
        gm_pvf[gm_background] = 0.
        csf_background = csf_pvf < $params.min_pvf_threshold
        csf_pvf[csf_background] = 0.

        data = np.concatenate(
            (wm_pvf[..., None], gm_pvf[..., None], csf_pvf[..., None]),
            axis=-1
        )
        maxes = np.argmax(data, axis=-1)

        wm_mask = maxes == 0
        wm_mask[wm_background] = False
        gm_mask = maxes == 1
        gm_mask[gm_background] = False
        csf_mask = maxes == 2
        csf_mask[csf_background] = False

        nib.save(nib.Nifti1Image(csf_mask.astype(np.uint8), affine), "${sid}_csf_mask.nii.gz")
        nib.save(nib.Nifti1Image(gm_mask.astype(np.uint8), affine), "${sid}_gm_mask.nii.gz")
        nib.save(nib.Nifti1Image(wm_mask.astype(np.uint8), affine), "${sid}_wm_mask.nii.gz")
        END_SCRIPT

        scil_image_math.py lower_threshold_eq $csf_pvf $params.max_safe_csf_pvf_threshold csf_map.nii.gz --data_type uint8
        if [[ $params.safe_csf_mask_dilation > 0 ]]
        then
            scil_image_math.py dilation csf_map.nii.gz $params.safe_csf_mask_dilation csf_map.nii.gz -f --data_type uint8
        fi
        scil_image_math.py lower_threshold_eq $gm_pvf $params.max_safe_gm_pvf_threshold gm_map.nii.gz --data_type uint8
        if [[ $params.safe_gm_mask_dilation > 0 ]]
        then
            scil_image_math.py dilation gm_map.nii.gz $params.safe_gm_mask_dilation gm_map.nii.gz -f --data_type uint8
        fi

        scil_image_math.py difference ${sid}_wm_mask.nii.gz csf_map.nii.gz ${sid}_safe_wm_mask.nii.gz
        scil_image_math.py difference ${sid}_safe_wm_mask.nii.gz gm_map.nii.gz ${sid}_safe_wm_mask.nii.gz -f
        scil_image_math.py intersection ${sid}_safe_wm_mask.nii.gz $brain_mask ${sid}_safe_wm_mask.nii.gz -f
        """
}

process crop_image {
    label "FAST"
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
            after_script += ["mrhardi fit2box --in ${image.simpleName}__cropped.nii.gz --out ${image.simpleName}__cropped.nii.gz --pbox $bounding_box"]
        }
        else
            args += "--output_bbox ${image.simpleName}__bbox.pkl"

        if ( !mask.empty() ) {
            before_script = "mrhardi apply_mask --in $image --mask $mask --out masked_image.nii.gz"
            img = "masked_image.nii.gz"
            mask_script = "mrhardi fit2box --in $mask --out ${mask.simpleName}__cropped.nii.gz"
            img_script = "mrhardi fit2box --in $image --out ${image.simpleName}__cropped.nii.gz"
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
            after_script += ["scil_image_math.py floor ${mask.simpleName}__cropped.nii.gz ${mask.simpleName}__cropped.nii.gz --data_type uint8 -f"]
        }

        if ( metadata instanceof nextflow.util.BlankSeparatedList ? !metadata.isEmpty() : !metadata.empty() )
            after_script += ["mrhardi metadata --in ${image.getSimpleName()}__cropped.nii.gz --update_affine --metadata $metadata"]

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
    label "LIGHTSPEED"
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
        mrhardi fitbox --in $image --ref $reference --pbox $bounding_box --out ${image.simpleName}__bbox
        """
}

process average {
    label "FAST"
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
        mrhardi concatenate --in ${images.join(",")} --out cat_images --ts
        fslmaths cat_images.nii.gz -Tmean ${base_name}__averaged.nii.gz
        """
}

process merge_masks {
    label "LIGHTSPEED"
    input:
        tuple val(sid), path(masks), val(base_name)
        val(caller_name)
    output:
        tuple val(sid), path("${base_name}__merged.nii.gz"), emit: mask
    script:
        """
        mrhardi concatenate --in ${masks.join(",")} --out cat_images --ts
        fslmaths cat_images.nii.gz -Tmax ${base_name}__merged.nii.gz
        """
}

process timeseries_mean {
    label "FAST"
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
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
        val(caller_name)
        path(config)
    output:
        tuple val(sid), path("${dwi.simpleName}__extracted_shells.nii.gz"), path("${dwi.simpleName}__extracted_shells.bval"), path("${dwi.simpleName}__extracted_shells.bvec"), emit: dwi
    script:
        def args = ""
        if (params.b0_threshold)
            args += " --ceil ${params.b0_threshold}"
        if (params.shell_threshold)
            args += " --gap ${params.shell_threshold}"
        """
        mrhardi shells --in $dwi --bvals $bval --bvecs $bvec --out ${dwi.simpleName}__extracted_shells --config $config $args
        """
}

process dilate_mask {
    label "LIGHTSPEED"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(mask)
        val(dilation_factor)
        val(caller_name)
    output:
        tuple val(sid), path("${mask.simpleName}__dilated.nii.gz"), emit: mask
    script:
        """
        scil_image_math.py dilation $mask $dilation_factor ${mask.simpleName}__dilated.nii.gz --data_type uint8
        """
}

process erode_mask {
    label "LIGHTSPEED"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all

    input:
        tuple val(sid), path(mask)
        val(erosion_factor)
        val(caller_name)
    output:
        tuple val(sid), path("${mask.simpleName}__dilated.nii.gz"), emit: mask
    script:
        """
        scil_image_math.py erosion $mask $erosion_factor ${mask.simpleName}__dilated.nii.gz --data_type uint8
        """
}

process clean_mask_borders {
    label "LIGHTSPEED"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(mask)
        val(factor)
        val(caller_name)
    output:
        tuple val(sid), path("${mask.simpleName}__clean_borders.nii.gz"), emit: mask
    script:
        """
        scil_image_math.py opening $mask $factor ${mask.simpleName}__clean_borders.nii.gz --data_type uint8 -f
        scil_image_math.py closing ${mask.simpleName}__clean_borders.nii.gz $factor ${mask.simpleName}__clean_borders.nii.gz --data_type uint8 -f
        """
}

process segmentation_to_binary {
    label "LIGHTSPEED"
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
        mrhardi seg2mask --in $segmentation --values 1,2,3,4 --labels csf,gm,dgm,wm --out ${segmentation.simpleName}
        scil_image_math.py addition ${segmentation.simpleName}_gm.nii.gz ${segmentation.simpleName}_dgm.nii.gz ${segmentation.simpleName}_all_gm.nii.gz --data_type uint8 -f
        """
}

process prepend_sid {
    label "LIGHTSPEED"
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
    label "LIGHTSPEED"
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
    label "FAST"
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
        mrhardi even_dimensions --in $dwi --suffix __even_dims $args
        cp $bval ${dwi.simpleName}__even_dims.bval
        cp $bvec ${dwi.simpleName}__even_dims.bvec
        $after_script
        """
}

process check_for_duplicates {
    label "FAST"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), file(metadata)
        val(caller_name)
    output:
        tuple val(sid), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.nii.gz"), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.bval"), path("${dwi.simpleName}__${params.duplicates_merge_method}_duplicates.bvec"), emit: dwi
        tuple val(sid), path("*__${params.duplicates_merge_method}_duplicates_metadata.*"), optional: true, emit: metadata   
    script:
        def args = ""
        if (params.b0_threshold)
            args += " --ceil ${params.b0_threshold}"
        """
        mrhardi duplicates --in $dwi --bvals $bval --bvecs $bvec --merge $params.duplicates_merge_method --out ${dwi.simpleName}__${params.duplicates_merge_method}_duplicates $args
        """
}

process validate_gradients {
    label "MEDIUM"
    label "res_single_cpu"

    publishDir "${params.output_root}/all/${sid}/$caller_name/${task.index}_${task.process.replaceAll(":", "_")}", mode: params.publish_mode, enabled: params.publish_all
    publishDir "${params.output_root}/${sid}", saveAs: { f -> f.contains("metadata") ? null : remove_alg_suffixes(f) }, mode: params.publish_mode

    input:
        tuple val(sid), path(bvec), path(peaks), path(fa), file(mask), file(peaks_vals)
        val(caller_name)
    output:
        tuple val(sid), path("${bvec.simpleName}__validated.bvec"), emit: bvecs
    script:
        def args = ""
        if ( !mask.empty() ) args += " --mask $mask"
        if ( !peaks_vals.empty() ) args += " --peaks_vals $peaks_vals"
        """
        scil_validate_and_correct_bvecs.py $bvec $peaks $fa ${bvec.simpleName}__validated.bvec --fa_th $params.validate_bvecs_fa_thr -f $args
        """

}