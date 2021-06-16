#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.merge_repetitions = false
params.eddy_on_rev = true
params.has_reverse =true
params.use_cuda = false

params.ants_registration_base_config = file("$projectDir/.config/ants_registration_base_config.py")
params.ants_transform_base_config = file("$projectDir/.config/ants_transform_base_config.py")
params.preproc_extract_b0_topup_config = file("$projectDir/.config/extract_mean_b0_base_config.py")
params.preproc_squash_b0_config = file("$projectDir/.config/preproc_squash_b0_config.py")
params.prepare_topup_base_config = file("$projectDir/.config/prepare_topup_base_config.py")
params.prepare_eddy_base_config = file("$projectDir/.config/prepare_eddy_base_config.py")
params.prepare_eddy_cuda_base_config = file("$projectDir/.config/prepare_eddy_cuda_base_config.py")
params.concatenate_base_config = file("$projectDir/.config/concatenate_base_config.py")

include { merge_channels_non_blocking; group_subject_reps; join_optional; map_optional; opt_channel; replace_dwi_file; uniformize_naming; sort_as_with_name; merge_repetitions; interleave; is_data } from '../functions.nf'
include { extract_b0 as b0_topup; extract_b0 as b0_topup_rev; squash_b0 as squash_dwi; squash_b0 as squash_rev } from '../processes/preprocess.nf'
include { n4_denoise; dwi_denoise; prepare_topup; topup; prepare_eddy; eddy } from '../processes/denoise.nf'
include { ants_register; ants_transform } from '../processes/register.nf'
include { cat_datasets; cat_datasets as cat_topup; cat_datasets as cat_eddy_on_rev; bet_mask; split_image; apply_topup; convert_datatype; replicate_image; check_dwi_conformity } from '../processes/utils.nf'


workflow registration_wkf {
    take:
        target_channel
        moving_channel
        trans_channel
        mask_channel
        bvecs_channel
        metadata_channel
        parameters
    main:
        reg_metadata = null
        trans_metadata = null
        if ( is_data(metadata_channel) ) {
            reg_metadata = metadata_channel.map { it.subList(0, it.size() - 1) }
            trans_metadata = metadata_channel.map { [it[0], it[-1]] }
        }
        into_register = moving_channel.join(target_channel).join(target_channel.map{ [it[0], it[1][0]] })
        into_register = join_optional(into_register, mask_channel)
        ants_register(join_optional(into_register, reg_metadata), "preprocess", parameters ? parameters : params.ants_registration_base_config)

        if ( is_data(trans_channel) ) {
            in_ants_trans = trans_channel.join(ants_register.out.reference).join(ants_register.out.transformation)
            in_ants_trans = join_optional(in_ants_trans, bvecs_channel)
            ants_transform(join_optional(in_ants_trans, trans_metadata), "preprocess", params.ants_transform_base_config)
            img = ants_transform.out.image
        }
        else {
            img = ants_register.out.image
        }
    emit:
        image = img
        registration = ants_register.out.image
        transform = ants_register.out.reference.join(ants_register.out.transformation)
}

// TODO : Here there is probably some metadatas from squashed process being tangled in i/o. The
// i/o bridge should be removed and the squashed metadatas should be passed directly to
// the delegated workflows/processes
workflow topup_wkf {
    take:
        dwi_channel
        rev_channel
        metadata_channel
    main:
        meta_channel = metadata_channel.map{ [it[0], [it[1]]] }.join(
            metadata_channel.map{ [it[0], [it[2]]] }
        ).map{ [it[0], it[1] + it[2]] }

        b0_topup(dwi_channel.map{ it.subList(0, 3) }.join(meta_channel), "preprocess", params.preproc_extract_b0_topup_config)
        b0_topup_rev(rev_channel.map{ it.subList(0, 3) }.join(meta_channel), "preprocess", params.preproc_extract_b0_topup_config)

        b0_channel = b0_topup.out.b0
        b0_rev_channel = b0_topup_rev.out.b0
        b0_metadata_channel = b0_topup.out.metadata.join(b0_topup_rev.out.metadata).map{ [it[0], it.subList(1, it.size())] }

        if ( params.merge_repetitions ) {
            b0_channel = merge_repetitions(b0_channel, false) //.map{ [it[0], it[1].inject([]){ c, t -> c + t }] }
            b0_rev_channel = merge_repetitions(b0_rev_channel, false) //.map{ [it[0], it[1].inject([]){ c, t -> c + t }] }
            b0_metadata_channel = merge_repetitions(b0_metadata_channel, false).map{ [it[0], it[1].inject([]){ c, t -> c + t }] }
            meta_channel = merge_repetitions(meta_channel, false).map{ [it[0], it[1].inject([]){ c, t -> c + t }] }
            dwi_channel = merge_repetitions(dwi_channel, false)
            rev_channel = merge_repetitions(rev_channel, false)
        }
        else {
            b0_channel = b0_channel.map{ [it[0], it.subList(1, it.size())] }
            b0_rev_channel = b0_rev_channel.map{ [it[0], it.subList(1, it.size())] }
        }

        acq_channel = dwi_channel.map{ [it[0], [it[2]]] }.join(rev_channel.map{ [it[0], [it[2]]] })
        b0_data_channel = b0_channel.join(b0_rev_channel).map{ [it[0], it.subList(1, it.size()).inject([]){ c, t -> c + t }] }
        b0_data_channel = b0_data_channel.map{ it + [[], []] }.join(b0_metadata_channel)

        cat_topup(b0_data_channel, "b0", "preprocess", params.concatenate_base_config)

        metadata_channel = cat_topup.out.metadata.join(meta_channel).map{ [it[0], [it[1]] + it[2]] }

        prepare_topup(cat_topup.out.image.join(dwi_channel.map{ [it[0], it[2]] }).join(rev_channel.map{ [it[0], it[2]] }).join(metadata_channel), params.prepare_topup_base_config)
        data_channel = prepare_topup.out.config.map{ it.subList(0, 4) }.join(cat_topup.out.image)

        topup(data_channel.join(prepare_topup.out.metadata), "preprocess")
    emit:
        b0 = topup.out.image
        field = topup.out.field
        movpar = topup.out.transfo.map{ [it[0], it[1]] }
        coeff = topup.out.transfo.map{ [it[0], it[2]] }
        param = prepare_topup.out.config.map{ [it[0], it[2]] }
        prefix = prepare_topup.out.config.map{ [it[0], it[-1]] }
        topup = topup.out.pkg
        metadata = prepare_topup.out.metadata
        in_metadata_w_topup = sort_as_with_name(prepare_topup.out.in_metadata_w_topup, acq_channel.map{ it.flatten() })
}

workflow apply_topup_wkf {
    take:
        dwi_channel
        rev_channel
        topup_channel
        meta_channel
    main:
        if ( params.merge_repetitions ) {
            // meta_channel = merge_repetitions(meta_channel, false).map{ [it[0], it[1].inject([]){ c, t -> c + t }] }
            dwi_channel = merge_repetitions(dwi_channel, false)
            rev_channel = merge_repetitions(rev_channel, false)
        }
        data_channel = dwi_channel.join(rev_channel.map{ it.subList(0, 2) })
        apply_topup(data_channel.join(topup_channel).join(meta_channel), "preprocess")
        dwi = apply_topup.out.dwi
        metadata = apply_topup.out.metadata
        if ( params.merge_repetitions ) {
            cat_datasets(dwi.join(metadata), "dwi", "preprocess", params.concatenate_base_config)
            dwi = cat_datasets.out.image.join(cat_datasets.out.bval).join(cat_datasets.out.bvec)
            metadata = cat_datasets.out.metadata
        }
    emit:
        dwi = dwi
        metadata = metadata
}

workflow squash_wkf {
    take:
        dwi_channel
        rev_channel
        metadata_channel
    main:
        dwi_meta_channel = metadata_channel.map{ it.subList(0, 2) }

        squash_dwi(dwi_channel.join(dwi_meta_channel), "preprocess", params.preproc_squash_b0_config)
        meta_channel = squash_dwi.out.metadata
        if ( params.has_reverse ) {
            rev_meta_channel = metadata_channel.map{ [it[0], it[2]] }
            squash_rev(rev_channel.join(rev_meta_channel), "preprocess", params.preproc_squash_b0_config)
            rev_channel = squash_rev.out.dwi
            meta_channel =meta_channel.join(squash_rev.out.metadata)
        }
    emit:
        dwi = squash_dwi.out.dwi
        rev = rev_channel
        metadata = meta_channel
}

workflow eddy_wkf {
    take:
        dwi_channel
        mask_channel
        topup_channel
        topup_b0_channel
        rev_channel
        metadata_channel
    main:

        bval_channel = dwi_channel.map{ [it[0], "${it[2].getName()}".tokenize(".")[0]] }
        if ( is_data(topup_channel) ) {
            bval_channel = join_optional(bval_channel, topup_channel.map { [it[0], it[1]] })
        }
        else {
            bval_channel = bval_channel.map{ it + [""] }
        }
        if ( params.has_reverse ) {
            bval_channel = join_optional(bval_channel, rev_channel.map { [it[0], "${it[2].getName()}".tokenize(".")[0]] })
        }
        else {
            bval_channel = bval_channel.map{ it + [""] }
        }

        metadata_channel = metadata_channel.map{ [it[0], it.subList(1, it.size())] }

        prep_eddy_channel = dwi_channel
        if ( params.has_reverse ) {
            prep_eddy_channel = dwi_channel.join(rev_channel)
        }

        in_prep = bval_channel.join(prep_eddy_channel.map{ [it[0], it.subList(1, it.size())] })
        prepare_eddy(
            bval_channel.join(prep_eddy_channel.map{ [it[0], it.subList(1, it.size())] }).join(metadata_channel),
            params.use_cuda ? params.prepare_eddy_cuda_base_config : params.prepare_eddy_base_config
        )

        if ( params.has_reverse ) {
            dwi_channel = dwi_channel.map{ it.subList(0, 3)}.join(prepare_eddy.out.bvec.map{ [it[0], it[1].find{ f -> f.simpleName.indexOf("_rev") == -1 }] })
            rev_channel = rev_channel.map { it.subList(0, 3) }.join(prepare_eddy.out.bvec.map { [it[0], it[1].find { f -> f.simpleName.indexOf("_rev") >= 0 }] })
        }
        else {
            dwi_channel = dwi_channel.map{ it.subList(0, 3)}.join(prepare_eddy.out.bvec)
        }

        if ( params.has_reverse && params.eddy_on_rev ) {
            cat_eddy_on_rev(merge_channels_non_blocking(dwi_channel, rev_channel).join(metadata_channel), "dwi", "preprocess", params.concatenate_base_config)
            dwi_channel = cat_eddy_on_rev.out.image.join( cat_eddy_on_rev.out.bval).join(cat_eddy_on_rev.out.bvec)
            metadata_channel = cat_eddy_on_rev.out.metadata.map{ [it[0], it.subList(1, it.size())] }
        }

        // if ( params.eddy_pre_denoise ) {
        //    dwi_denoise_wkf(dwi_channel.map{ it.subList(0, 2) }, null, metadata_channel)
        //    dwi_channel = replace_dwi_file(dwi_channel, dwi_denoise_wkf.out.image)
        //    metadata_channel = dwi_denoise_wkf.out.metadata.map{ [it[0], it.subList(1, it.size())] }
        // }

        dwi_channel = join_optional(dwi_channel, mask_channel)
        if ( params.has_reverse ) {
            dwi_channel = join_optional(dwi_channel, topup_channel.map { [it[0], it[2], it[3]] })
        }
        else {
            dwi_channel = dwi_channel.map{ it + ["", ""] }
        }

        eddy_in = prepare_eddy.out.config

        if ( params.use_cuda )
            eddy_in = eddy_in.join(prepare_eddy.out.slspec)
        else
            eddy_in = eddy_in.map{ it + [""] }

        eddy(eddy_in.join(dwi_channel).join(metadata_channel), "preprocess")
        check_dwi_conformity(eddy.out.dwi.join(eddy.out.bval).join(eddy.out.bvec).join(eddy.out.metadata), "fix", "preprocess")
    emit:
        dwi = check_dwi_conformity.out.dwi.map{ [it[0], it[1]] }
        bval = check_dwi_conformity.out.dwi.map{ [it[0], it[2]] }
        bvec = check_dwi_conformity.out.dwi.map{ [it[0], it[3]] }
        metadata = check_dwi_conformity.out.metadata
}

workflow dwi_denoise_wkf {
    take:
        dwi_channel
        mask_channel
        metadata_channel
    main:
        dwi_channel = join_optional(dwi_channel, mask_channel)
        dwi_denoise(dwi_channel.join(metadata_channel), "preprocess")
    emit:
        image = dwi_denoise.out.image
        metadata = dwi_denoise.out.metadata
}

workflow n4_denoise_wkf {
    take:
        dwi_channel
        b0_channel
        mask_channel
        metadata_channel
        config
    main:
        dwi_channel = join_optional(dwi_channel, b0_channel)
        dwi_channel = join_optional(dwi_channel, mask_channel)
        dwi_channel = join_optional(dwi_channel, metadata_channel)
        n4_denoise(dwi_channel, "preprocess", config)
    emit:
        image = n4_denoise.out.image
        metadata = n4_denoise.out.metadata
}