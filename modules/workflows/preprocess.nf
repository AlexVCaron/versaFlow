#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.eddy_with_reverse = true
params.use_cuda = false

params.ants_registration_base_config = file("$projectDir/.config/ants_registration_base_config.py")
params.ants_transform_base_config = file("$projectDir/.config/ants_transform_base_config.py")
params.preproc_extract_b0_topup_config = file("$projectDir/.config/extract_mean_b0_base_config.py")
params.preproc_squash_b0_config = file("$projectDir/.config/preproc_squash_b0_config.py")
params.prepare_topup_base_config = file("$projectDir/.config/prepare_topup_base_config.py")
params.prepare_eddy_base_config = file("$projectDir/.config/prepare_eddy_base_config.py")
params.prepare_eddy_cuda_base_config = file("$projectDir/.config/prepare_eddy_cuda_base_config.py")
params.concatenate_base_config = file("$projectDir/.config/concatenate_base_config.py")

include {
    filter_datapoints; separate_b0_from_dwi; exclude_missing_datapoints; fill_missing_datapoints;
    merge_channels_non_blocking; join_optional; sort_as_with_name; is_data
} from '../functions.nf'
include { extract_b0 as b0_topup; extract_b0 as b0_topup_rev; squash_b0 as squash_dwi; squash_b0 as squash_rev } from '../processes/preprocess.nf'
include { n4_denoise; dwi_denoise; nlmeans_denoise; prepare_topup; topup; prepare_eddy; eddy } from '../processes/denoise.nf'
include { ants_register; ants_transform } from '../processes/register.nf'
include {
    cat_datasets; cat_datasets as cat_topup; cat_datasets as cat_eddy_on_rev;
    apply_topup; check_dwi_conformity; generate_b0_bval
} from '../processes/utils.nf'


workflow registration_wkf {
    take:
        target_channel
        moving_channel
        trans_channel
        mask_channel
        bvecs_channel
        metadata_channel
        additional_publish_path
        publish
        publish_suffix
        trans_publish_suffix
        registration_parameters
        transformation_parameters
    main:
        reg_metadata = null
        trans_metadata = null
        if ( is_data(metadata_channel) ) {
            reg_metadata = metadata_channel.map { it.subList(0, it.size() - 1) }
            trans_metadata = metadata_channel.map { [it[0], it[-1]] }
        }
        into_register = moving_channel.join(target_channel).join(target_channel.map{ [it[0], it[1][0]] })
        into_register = join_optional(into_register, mask_channel)
        ants_register(
            join_optional(into_register, reg_metadata),
            "preprocess",
            additional_publish_path,
            publish, publish_suffix,
            registration_parameters ? registration_parameters : params.ants_registration_base_config
        )

        if ( is_data(trans_channel) ) {
            in_ants_trans = trans_channel.join(ants_register.out.reference).join(ants_register.out.transformation)
            in_ants_trans = join_optional(in_ants_trans, bvecs_channel)
            ants_transform(
                join_optional(in_ants_trans, trans_metadata),
                "preprocess",
                additional_publish_path,
                publish, trans_publish_suffix,
                transformation_parameters ? transformation_parameters : params.ants_transform_base_config
            )
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
        existing_rev = exclude_missing_datapoints(rev_channel, 1, "")

        (dwi_rev, b0_rev) = separate_b0_from_dwi(existing_rev)
        topupable_indexes = existing_rev.map{ [it[0]] }

        topupable_meta_channel = topupable_indexes.join(metadata_channel.map{ [it[0], [it[1], it[2]]] })
        topupable_dwi_channel = topupable_indexes.join(dwi_channel)

        b0_topup(topupable_dwi_channel.map{ it.subList(0, 3) }.join(topupable_meta_channel), "preprocess", "false", params.preproc_extract_b0_topup_config)
        b0_topup_rev(dwi_rev.map{ it.subList(0, 3) }.join(topupable_meta_channel), "preprocess", "false", params.preproc_extract_b0_topup_config)

        b0_channel = b0_topup.out.b0
        b0_rev_channel = b0_topup_rev.out.b0.mix(b0_rev.map{ it.subList(0, 2) })
        b0_metadata_channel = b0_topup.out.metadata.join(
            b0_topup_rev.out.metadata.mix(b0_rev.map{ [it[0]] }.join(topupable_meta_channel).map{ [it[0], it[1][1]] })
        ).map{ [it[0], it.subList(1, it.size())] }

        b0_channel = b0_channel.map{ [it[0], it.subList(1, it.size())] }
        b0_rev_channel = b0_rev_channel.map{ [it[0], it.subList(1, it.size())] }

        acq_channel = topupable_dwi_channel.map{ [it[0], [it[2]]] }.join(existing_rev.map{ [it[0], [it[2]]] })
        b0_data_channel = b0_channel.join(b0_rev_channel).map{ [it[0], it.subList(1, it.size()).inject([]){ c, t -> c + t }] }
        b0_data_channel = b0_data_channel.map{ it + [[], []] }.join(b0_metadata_channel)

        cat_topup(b0_data_channel, "b0", "preprocess", params.concatenate_base_config)

        metadata_channel = cat_topup.out.metadata.join(topupable_meta_channel).map{ [it[0], [it[1]] + it[2]] }

        generate_b0_bval(b0_rev.map{ it.subList(0, 2) }, "false")

        prepare_topup(
            cat_topup.out.image.join(
                topupable_dwi_channel.map{ [it[0], it[2]] }
            ).join(
                dwi_rev.map{ [it[0], it[2]] }.mix(generate_b0_bval.out.bval)
            ).join(metadata_channel),
            params.prepare_topup_base_config
        )
        data_channel = prepare_topup.out.config.map{ it.subList(0, 4) }.join(cat_topup.out.image)

        topup(data_channel.join(prepare_topup.out.metadata), "preprocess")
        excluded_indexes = filter_datapoints(rev_channel, { it[1] == "" })
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
        topupable_indexes = topupable_indexes
        excluded_indexes = excluded_indexes
        excluded_dwi = excluded_indexes.join(dwi_channel)
        excluded_dwi_metadata = excluded_indexes.join(metadata_channel).map{ it.subList(0, 2) }
}

workflow apply_topup_wkf {
    take:
        dwi_channel
        rev_channel
        topup_channel
        meta_channel
        additional_publish_path
    main:
        data_channel = dwi_channel.join(rev_channel.map{ it.subList(0, 2) })
        apply_topup(data_channel.join(topup_channel).join(meta_channel), "preprocess", additional_publish_path)
        dwi = apply_topup.out.dwi
        metadata = apply_topup.out.metadata
    emit:
        dwi = dwi
        metadata = metadata
}

workflow squash_wkf {
    take:
        dwi_channel
        rev_channel
        metadata_channel
        additional_publish_path
    main:

        squash_dwi(dwi_channel.join(metadata_channel.map{ it.subList(0, 2) }), "preprocess", "true", params.preproc_squash_b0_config, additional_publish_path)

        (dwi_rev, b0_rev) = separate_b0_from_dwi(exclude_missing_datapoints(rev_channel.join(metadata_channel.map{ [it[0], it[2]] }), 1, ""))

        squash_rev(dwi_rev, "preprocess", "false", params.preproc_squash_b0_config, additional_publish_path)

        ref_id_channel = dwi_channel.map{ [it[0]] }
    emit:
        dwi = squash_dwi.out.dwi
        rev = fill_missing_datapoints(squash_rev.out.dwi.mix(b0_rev.map{ it.subList(0, it.size() - 1) }), ref_id_channel, 1, ["", "", ""])
        metadata = squash_dwi.out.metadata.join(fill_missing_datapoints(squash_rev.out.metadata.mix(b0_rev.map{ [it[0], it[-1]] }), ref_id_channel, 1, [""] ))
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
        ref_id_channel = dwi_channel.map{ [it[0]] }

        bval_channel = dwi_channel.map{ [it[0], "${it[2].getName()}".tokenize(".")[0]] }
        bval_channel = bval_channel.join(topup_channel.map { [it[0], it[1]] })

        absent_reverse_ids = filter_datapoints(rev_channel, { it[1] == "" }).map{ [it[0]] }
        rev_channel = exclude_missing_datapoints(rev_channel, 1, "")

        (dwi_rev, b0_rev) = separate_b0_from_dwi(rev_channel)

        generate_b0_bval(b0_rev.map{ it.subList(0, 2) }, "true")
        b0_bval = generate_b0_bval.out.bval
        b0_bvec = generate_b0_bval.out.bvec

        rev_channel = fill_missing_datapoints(
            dwi_rev.mix(b0_rev.map{ it.subList(0, 2) }.join(b0_bval).join(b0_bvec)),
            ref_id_channel,
            1, ["", "", ""]
        )

        bval_channel = bval_channel.join(
            fill_missing_datapoints(
                dwi_rev.mix(b0_bval).map{ [it[0], "${it[1].getName()}".tokenize(".")[0]] },
                ref_id_channel,
                1, [""]
            )
        )

        metadata_channel = metadata_channel.map{ [it[0], it.subList(1, it.size())] }

        prepare_eddy(
            bval_channel.join(dwi_channel.join(rev_channel).map{ [it[0], it.subList(1, it.size())] }).join(metadata_channel),
            params.use_cuda ? params.prepare_eddy_cuda_base_config : params.prepare_eddy_base_config
        )

        dwi_channel = dwi_channel.map{ it.subList(0, 3) }.join(prepare_eddy.out.bvec.map{ [it[0], it[1].find{ f -> f.simpleName.indexOf("_rev") == -1 }] })
        rev_channel = rev_channel.map{ it.subList(0, 3) }.join(prepare_eddy.out.bvec.map{ [it[0], it[1].find { f -> f.simpleName.indexOf("_rev") >= 0 }] })

        if ( params.eddy_with_reverse ) {
            cat_eddy_on_rev(merge_channels_non_blocking(dwi_channel, rev_channel).join(metadata_channel), "dwi", "preprocess", params.concatenate_base_config)
            dwi_channel = cat_eddy_on_rev.out.image.join( cat_eddy_on_rev.out.bval ).join(cat_eddy_on_rev.out.bvec).mix( absent_reverse_ids.join(dwi_channel) )
            metadata_channel = cat_eddy_on_rev.out.metadata.map{ [it[0], it.subList(1, it.size())] }.mix( absent_reverse_ids.join(metadata_channel) )
        }

        dwi_channel = dwi_channel.join(mask_channel).join(topup_channel.map { [it[0], it[2], it[3]] })
        eddy_in = prepare_eddy.out.config.join(prepare_eddy.out.slspec)

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
        publish
    main:
        (dwi_into_denoise, b0_into_denoise) = separate_b0_from_dwi(exclude_missing_datapoints(dwi_channel.join(mask_channel).join(metadata_channel), 1, ""))
        dwi_denoise(dwi_into_denoise.map{ it.subList(0, 2) + it.subList(4, it.size()) }, "preprocess", "$publish")
        nlmeans_denoise(b0_into_denoise.map{ it.subList(0, 2) + it.subList(4, it.size()) }, "preprocess", "$publish")

        ref_id_channel = dwi_channel.map{ [it[0]] }
    emit:
        image = fill_missing_datapoints(dwi_denoise.out.image.mix(nlmeans_denoise.out.image), ref_id_channel, 1, [""])
        metadata = fill_missing_datapoints(dwi_denoise.out.metadata.mix(nlmeans_denoise.out.metadata), ref_id_channel, 1, [""])
}

workflow n4_denoise_wkf {
    take:
        image_channel
        ref_anat_channel
        mask_channel
        metadata_channel
        config
    main:
        ref_id_channel = image_channel.map{ [it[0]] }
        n4_denoise(
            image_channel.join(
                fill_missing_datapoints(ref_anat_channel, ref_id_channel, 1, [""])
            ).join(
                fill_missing_datapoints(mask_channel, ref_id_channel, 1, [""])
            ).join(
                fill_missing_datapoints(metadata_channel, ref_id_channel, 1, [""])
            ),
            "preprocess",
            config
        )
    emit:
        image = n4_denoise.out.image
        metadata = n4_denoise.out.metadata
}