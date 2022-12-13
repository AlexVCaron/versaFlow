#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    filter_datapoints;
    separate_b0_from_dwi;
    exclude_missing_datapoints;
    fill_missing_datapoints;
    merge_channels_non_blocking;
    join_optional;
    sort_as_with_name;
    is_data;
    get_config_path;
    collect_paths;
    is_path_list
} from '../functions.nf'
include {
    extract_b0 as b0_topup;
    extract_b0 as b0_topup_rev;
    squash_b0 as squash_dwi;
    squash_b0 as squash_rev
} from '../processes/preprocess.nf'
include {
    n4_denoise;
    apply_n4_bias_field;
    dwi_denoise;
    nlmeans_denoise;
    nlmeans_denoise as nlmeans_denoise_b0_from_fwd_dwi;
    nlmeans_denoise as nlmeans_denoise_b0_from_rev_dwi;
    prepare_topup;
    topup;
    prepare_eddy;
    eddy
} from '../processes/denoise.nf'
include {
    ants_register;
    ants_transform;
    align_to_closest as b0_align_to_closest;
    align_to_closest as rev_align_to_closest;
    align_to_average as b0_align_to_average
} from '../processes/register.nf'
include {
    cat_datasets as cat_topup;
    cat_datasets as cat_eddy_on_rev;
    cat_datasets as concatenate_for_average;
    cat_datasets as concatenate_b0;
    cat_datasets as concatenate_rev_b0;
    timeseries_mean as get_average;
    apply_topup;
    check_dwi_conformity;
    generate_b0_bval;
    split_image as split_b0;
    split_image as split_rev
} from '../processes/utils.nf'

params.eddy_with_reverse = true
params.use_cuda = false

params.ants_registration_base_config = file("${get_config_path()}/ants_registration_base_config.py")
params.ants_transform_base_config = file("${get_config_path()}/ants_transform_base_config.py")
params.preproc_extract_b0_topup_config = file("${get_config_path()}/extract_mean_b0_base_config.py")
params.preproc_squash_b0_config = file("${get_config_path()}/preproc_squash_b0_config.py")
params.prepare_topup_base_config = file("${get_config_path()}/prepare_topup_base_config.py")
params.prepare_eddy_base_config = file("${get_config_path()}/prepare_eddy_base_config.py")
params.prepare_eddy_cuda_base_config = file("${get_config_path()}/prepare_eddy_cuda_base_config.py")
params.concatenate_base_config = file("${get_config_path()}/concatenate_base_config.py")


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
            in_ants_trans = trans_channel
                .join(ants_register.out.reference)
                .join(ants_register.out.transformation.map{ [
                    it[0],
                    it[1] instanceof ArrayList ? it[1].reverse() : [it[1]]
                ] })
                .map{ it + [""] }

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
        reference = ants_register.out.reference
        transform = ants_register.out.transformation.map{ [
            it[0], it[1] instanceof ArrayList ? it[1].reverse() : [it[1]]
        ] }
        inverse_transform = ants_register.out.inverse_transformation
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

        align_b0_for_topup(
            b0_channel,
            b0_rev_channel,
            b0_topup.out.metadata,
            b0_topup_rev.out.metadata.mix(b0_rev.map{ [it[0]] }.join(topupable_meta_channel).map{ [it[0], it[1][1]] })
        )

        b0_channel = align_b0_for_topup.out.b0.map{ [it[0], it.subList(1, it.size())] }
        b0_rev_channel = align_b0_for_topup.out.rev_b0.map{ [it[0], it.subList(1, it.size())] }

        acq_channel = topupable_dwi_channel.map{ [it[0], [it[2]]] }.join(existing_rev.map{ [it[0], [it[2]]] })
        b0_data_channel = b0_channel
            .join(b0_rev_channel)
            .map{ [it[0], it.subList(1, it.size()).inject([]){ c, t -> c + t }] }
            .map{ it + [[], []] }
            .join(align_b0_for_topup.out.metadata)

        cat_topup(b0_data_channel, 3, "b0", "preprocess", params.concatenate_base_config)

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

workflow align_b0_for_topup {
    take:
        b0_channel
        b0_rev_channel
        b0_meta_channel
        b0_rev_meta_channel
    main:
        meta_channel = b0_meta_channel
            .join(b0_rev_meta_channel)
            .map{ [it[0], it.subList(1, it.size())] }

        split_b0(b0_channel.join(b0_meta_channel), 3, "preprocess")
        split_rev(b0_rev_channel.join(b0_rev_meta_channel), 3, "preprocess")

        b0_align_to_closest(
            split_b0.out.images
                .map{ it[1] instanceof Path ? [it[0], [it[1]]] : it  }
                .join(split_b0.out.metadata),
            1,
            false,
            "preprocess",
            "",
            false,
            ""
        )
        rev_align_to_closest(
            split_rev.out.images
                .map{ it[1] instanceof Path ? [it[0], [it[1]]] : it  }
                .join(split_rev.out.metadata),
            1,
            false,
            "preprocess",
            "",
            false,
            ""
        )

        b0_data_channel = b0_align_to_closest.out.images
            .join(rev_align_to_closest.out.images)
            .map{ [it[0], it.subList(1, it.size()).inject([]){ c, t -> t instanceof Path ? c + [t] : c + t }] }

        concatenate_for_average(b0_data_channel.map{ it + [[], []] }.join(meta_channel), 3, "b0", "preprocess", params.concatenate_base_config)
        get_average(concatenate_for_average.out.image, "preprocess")

        metadata_channel = b0_align_to_closest.out.metadata.join(rev_align_to_closest.out.metadata).map{ [it[0], [it[1], it[2]]] }
        b0_align_to_average(b0_data_channel.join(get_average.out.image).join(metadata_channel), 1, false, "preprocess", "", false, "")

        b0_map = b0_align_to_average.out.images
            .join(b0_align_to_closest.out.images)
            .multiMap{ it ->
                forward: [it[0], it[2] instanceof Path ? [it[1][0]] : it[1].subList(0, it[2].size())]
                reverse: [it[0], it[2] instanceof Path ? [it[1][1]] : it[1].subList(it[2].size(), it[1].size())] 
            }

        concatenate_b0(b0_map.forward.map{ it + [[], []] }.join(b0_align_to_average.out.metadata), 3, "b0__aligned", "preprocess", params.concatenate_base_config)
        concatenate_rev_b0(b0_map.reverse.map{ it + [[], []] }.join(b0_align_to_average.out.metadata), 3, "b0_rev__aligned", "preprocess", params.concatenate_base_config)
    emit:
        b0 = concatenate_b0.out.image
        rev_b0 = concatenate_rev_b0.out.image
        metadata = concatenate_b0.out.metadata
            .join(concatenate_rev_b0.out.metadata)
            .map{ [it[0], it[1..-1]] }
}

workflow apply_topup_wkf {
    take:
        dwi_channel
        rev_channel
        topup_channel
        meta_channel
        additional_publish_path
    main:
        data_channel = dwi_channel.join(rev_channel.map{ it[0..1] })
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
        absent_reverse_ids = filter_datapoints(rev_channel, { it[1] == "" })
            .map{ [it[0]] }
        reverse_ids = filter_datapoints(rev_channel, { it[1] })
            .map{ [it[0]] }

        rev_channel = exclude_missing_datapoints(rev_channel, 1, "")
        (dwi_rev, b0_rev) = separate_b0_from_dwi(rev_channel)

        reverse_dwi_ids = dwi_rev.map{ [it[0]] }
        reverse_b0_ids = b0_rev.map{ [it[0]] }

        generate_b0_bval(b0_rev.map{ it[0..1] }, "true")
        b0_bval = generate_b0_bval.out.bval
        b0_bvec = generate_b0_bval.out.bvec

        rev_channel = fill_missing_datapoints(
            b0_rev.map{ it[0..1] }
                .join(b0_bval)
                .join(b0_bvec)
                .mix(dwi_rev),
            ref_id_channel,
            1, ["", "", ""]
        )

        rev_prefix_channel = rev_channel
            .map{ [it[0], it[1] ? it[1].simpleName.tokenize(".")[0] : ""] }

        prepare_eddy(
            dwi_channel
                .map{ [it[0], it[1].simpleName.tokenize(".")[0]] }
                .join(topup_channel.map{ it[0..1] })
                .join(rev_prefix_channel)
                .join(collect_paths(dwi_channel.join(rev_channel)))
                .join(collect_paths(metadata_channel)),
            params.use_cuda ? params.prepare_eddy_cuda_base_config : params.prepare_eddy_base_config
        )

        non_zero_bvec = collect_paths(prepare_eddy.out.bvec.map{ it.flatten() })

        dwi_channel = dwi_channel
            .map{ it[0..2] }
            .join(non_zero_bvec.map{ [it[0], it[1].find{ f -> !f.simpleName.contains("_rev")}] })
        rev_channel = rev_channel
            .map{ it[0..2] }
            .join(non_zero_bvec.map{ [it[0], it[1].find{ f -> f.simpleName.contains("_rev")}] })
            .mix(absent_reverse_ids.mix(reverse_b0_ids).join(rev_channel))
        dwi_metadata_channel = collect_paths(metadata_channel)
            .map{ [it[0], it[1].find{ f -> !f.simpleName.contains("_rev") }] }

        if ( params.eddy_with_reverse ) {
            cat_eddy_on_rev(
                reverse_ids
                    .join(collect_paths(dwi_channel))
                    .join(collect_paths(rev_channel))
                    .map{ [it[0]] + it[1..-1].transpose() }
                    .join(collect_paths(metadata_channel)),
                3, "dwi", "preprocess",
                params.concatenate_base_config
            )

            dwi_channel = cat_eddy_on_rev.out.image
                .join(cat_eddy_on_rev.out.bval)
                .join(cat_eddy_on_rev.out.bvec)
                .mix(absent_reverse_ids.join(dwi_channel))

            dwi_metadata_channel = cat_eddy_on_rev.out.metadata
                .mix(absent_reverse_ids.join(dwi_metadata_channel))
        }

        slspec_channel = fill_missing_datapoints(
            prepare_eddy.out.slspec,
            ref_id_channel,
            1, ""
        )

        eddy(
            prepare_eddy.out.config
                .join(slspec_channel)
                .join(dwi_channel)
                .join(mask_channel)
                .join(topup_channel.map{ [it[0], it[2], it[3]] })
                .join(dwi_metadata_channel),
            "preprocess"
        )

        check_dwi_conformity(
            eddy.out.dwi
                .join(eddy.out.bval)
                .join(eddy.out.bvec)
                .join(eddy.out.metadata),
            "fix",
            "preprocess"
        )

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
        publish
    main:
        ref_id_channel = image_channel.map{ [it[0]] }
        absent_ref_anat_id_channel = filter_datapoints(ref_anat_channel, { it[1] == "" })
            .map{ [it[0]] }

        ref_anat_channel = exclude_missing_datapoints(ref_anat_channel, 1, "")
        mask_channel = fill_missing_datapoints(mask_channel, ref_id_channel, 1, [""])

        n4_denoise(
            absent_ref_anat_id_channel
                .join(image_channel)
                .mix(ref_anat_channel)
                .map{ it + [""] }
                .join(mask_channel)
                .join(fill_missing_datapoints(metadata_channel, ref_id_channel, 1, [""])),
            "preprocess",
            publish
            config
        )
        apply_n4_bias_field(
            ref_anat_channel
                .map{ [it[0]] }
                .join(image_channel)
                .join(n4_denoise.out.bias_field)
                .join(mask_channel)
                .join(n4_denoise.out.metadata),
            "preprocess",
            publish
        )
    emit:
        reference = n4_denoise.out.image
        image = absent_ref_anat_id_channel
            .join(n4_denoise.out.image)
            .mix(apply_n4_bias_field.out.image)
        metadata = absent_ref_anat_id_channel
            .join(n4_denoise.out.metadata)
            .mix(apply_n4_bias_field.out.metadata)
}