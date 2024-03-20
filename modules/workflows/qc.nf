#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import groovy.json.*

params.qc_gif_slice_gap = 5
params.qc_image_size = [2000, 1000]
params.qc_screenshot_orientation = "axial"
params.qc_extra_overlay_opacity = 0.5 // overlay is no the primary QC target
params.qc_primary_overlay_opacity = 0.7 // overlay is the primary QC target
params.qc_masks_as_contours = true

params.qc_dwi_gif_framerate = 0.6 // display an image each 2/3 of a second
params.qc_masking_gif_framerate = 0.6
params.qc_dti_gif_framerate = 0.6

params.qc_pvf_cmap_name = "jet"
params.qc_pvf_opacity = 0.5

params.qc_gradients_sphere = false
params.qc_gradients_disable_sphere = false
params.qc_gradients_disable_supershell = false
params.qc_gradients_split_shells = false
params.qc_gradients_sphere_opacity = 0.5
params.qc_gradients_shells_same_color = false

params.qc_b0_colormap_name = false
params.qc_dwi_colormap_name = false
params.qc_t1_colormap_name = false
params.qc_masks_colormap_name = false
params.qc_pvf_colormap_name = false
params.qc_d99_colormap_name = false
params.qc_inia19_colormap_name = false
params.qc_charm_colormap_name = false
params.qc_sarm_colormap_name = false

screenshot_axes = [
    "axial": 2,
    "coronal": 1,
    "sagittal": 0
]

include {
    extract_labels_from_atlas as extract_mask_from_atlas;
    image_bounding_box;
    transform_bounding_box;
    extract_slice_from_volume as extract_from_raw_volume;
    extract_slice_from_volume as extract_from_processed_volume;
    average_to_3d_volume
} from "../processes/utils.nf"

include {
    screenshot_3d_nifti as screenshot_dwi_volumes;
    screenshot_3d_nifti as screenshot_dwi_mask;
    screenshot_3d_nifti as screenshot_t1_mask;
    screenshot_3d_nifti as screenshot_tissues_masks;
    screenshot_3d_nifti as screenshot_tissues_pvfs;
    screenshot_3d_nifti as screenshot_dti_metrics;
    screenshot_3d_nifti as screenshot_fodf_metrics;
    screenshot_3d_nifti as screenshot_fodf_peaks_coronal;
    screenshot_3d_nifti as screenshot_fodf_peaks_sagittal;
    screenshot_gradient_sampling_fsl;
    stitch_images as stitch_dwi_raw_vs_processed;
    stitch_images as stitch_gradients_raw_vs_processed;
    stitch_images as stitch_tissues_screenshots;
    stitch_images as stitch_peaks_screenshots;
    compose_gif_from_images as compose_dwi_qc_gifs;
    compose_gif_from_images as compose_mask_qc_gifs;
    compose_gif_from_images as compose_tissues_qc_gifs;
    compose_gif_from_images as compose_dti_qc_gifs;
    compose_gif_from_images as compose_fodf_qc_gifs;
    generate_report_from_screenshots as generate_dwi_processing_report;
    generate_report_from_screenshots as generate_gradients_processing_report;
    generate_report_from_screenshots as generate_dwi_masking_report;
    generate_report_from_screenshots as generate_t1_masking_report;
    generate_report_from_screenshots as generate_tissues_masks_report;
    generate_report_from_screenshots as generate_tissues_pvf_report;
    generate_report_from_screenshots as generate_safe_wm_masks_report;
    generate_report_from_screenshots as generate_dti_report;
    generate_report_from_screenshots as generate_fodf_report
} from "../processes/qc.nf"

include {
    get_natural_sorter
} from "../functions.nf"


workflow qc_screenshot_parameters_wkf {
    take:
        segmentation_atlases_channel
    main:
        atlases_channel = segmentation_atlases_channel
            .map{ [it[0], it[1]] }
            .transpose()
            .branch{
                inia19: it[1].simpleName =~ /INIA19_atlas/
                d99: it[1].simpleName =~ /D99_atlas/
                charm: it[1].simpleName =~ /CHARM_atlas/
                sarm: it[1].simpleName =~ /SARM_atlas/
            }

        extract_mask_from_atlas(
            atlases_channel.inia19,
            [124, 1124],
            false,
            ""
        )

        image_bounding_box(
            extract_mask_from_atlas.out.mask,
            false,
            ""
        )

    emit:
        cc_bounding_box = image_bounding_box.out.box
}


workflow dti_qc_wkf {
    take:
        t1_channel
        dti_peaks_channel
        dti_metrics_channel
        masks_3t_channel
    main:
        masks_3t_splitting_channel = masks_3t_channel
            .map{ [it[0], it[1..-1]] }
            .transpose()
            .branch{
                wm_mask: it[1].simpleName =~ /(?<!safe)_wm_mask/
                gm_mask: it[1].simpleName =~ /_gm_mask/
                csf_mask: it[1].simpleName =~ /_csf_mask/
            }

        masks_3t_to_ovelays_channel = masks_3t_splitting_channel.wm_mask
            .join(masks_3t_splitting_channel.gm_mask)
            .join(masks_3t_splitting_channel.csf_mask)
            .map{ [it[0], it[1..-1]] }

        dti_metrics_splitting_channel = dti_metrics_channel
            .map{ [it[0], it[1..-1]] }
            .transpose()
            .branch{
                md: it[1].simpleName =~ /dti_md/
                    return [it[0], "dti_md"] + it[1..-1]
                ad: it[1].simpleName =~ /dti_ad/
                    return [it[0], "dti_ad"] + it[1..-1]
                rd: it[1].simpleName =~ /dti_rd/
                    return [it[0], "dti_rd"] + it[1..-1]
                fa: it[1].simpleName =~ /dti_fa/
                    return [it[0], "dti_fa"] + it[1..-1]
                rgb: it[1].simpleName =~ /dti_rgb/
                    return [it[0], "dti_rgb"] + it[1..-1]
                residuals: it[1].simpleName =~ /dti_residuals/
                    return [it[0],  "dti_residuals"] + it[1..-1]
            }

        dti_peaks_splitting_channel = t1_channel
            .join(masks_3t_splitting_channel.wm_mask)
            .join(dti_peaks_channel)
            .map{ [it[0], it[1], it[2], it[3..-1]] }
            .transpose()
            .branch{
                peaks: it[3].simpleName =~ /dti_evecs_v1/
                    return [it[0], "dti_peaks", it[1]] + it[2..-1]
            }

        md_to_shot_channel = dti_metrics_splitting_channel.md
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        ad_to_shot_channel = dti_metrics_splitting_channel.ad
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        rd_to_shot_channel = dti_metrics_splitting_channel.rd
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        fa_to_shot_channel = dti_metrics_splitting_channel.fa
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        rgb_to_shot_channel = dti_metrics_splitting_channel.rgb
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        residuals_to_shot_channel = dti_metrics_splitting_channel.residuals
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        peaks_to_shot_channel = dti_peaks_splitting_channel.peaks
            //.join(masks_3t_to_ovelays_channel)
            .map{ it[0..-3] + [[it[-2]], [], [], it[-1]] }

        screenshot_dti_metrics(
            md_to_shot_channel
                .mix(ad_to_shot_channel)
                .mix(rd_to_shot_channel)
                .mix(fa_to_shot_channel)
                .mix(rgb_to_shot_channel)
                .mix(residuals_to_shot_channel)
                .mix(peaks_to_shot_channel),
            [params.qc_screenshot_orientation, params.qc_extra_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [true, true],
            params.qc_image_size,
            [false, false, false]
        )

        compose_dti_qc_gifs(
            screenshot_dti_metrics.out.screenshots
                .map{ it[0..1] + [it[2].sort(get_natural_sorter())] },
            params.qc_image_size,
            [params.qc_dti_gif_framerate, 0, false, false]
        )

        generate_dti_report(
            compose_dti_qc_gifs.out.gif
                .collect{ it[1] }
                .map{ [it, []] },
            "dti_metrics",
            [false, true]
        )

}


workflow fodf_qc_wkf {
    take:
        t1_channel
        fodf_peaks_channel
        fodf_metrics_channel
        masks_3t_channel
    main:
        masks_3t_splitting_channel = masks_3t_channel
            .map{ [it[0], it[1..-1]] }
            .transpose()
            .branch{
                wm_mask: it[1].simpleName =~ /(?<!safe)_wm_mask/
                gm_mask: it[1].simpleName =~ /_gm_mask/
                csf_mask: it[1].simpleName =~ /_csf_mask/
            }

        masks_3t_to_ovelays_channel = masks_3t_splitting_channel.wm_mask
            .join(masks_3t_splitting_channel.gm_mask)
            .join(masks_3t_splitting_channel.csf_mask)
            .map{ [it[0], it[1..-1]] }

        fodf_metrics_splitting_channel = fodf_metrics_channel
            .transpose()
            .map{ [it[0], it[1..-1]] }
            .transpose()
            .branch{
                afd: it[1].simpleName =~ /fodf_metrics.*_afd(?![st])/
                    return [it[0], "fodf_afd", it[1]]
                afds: it[1].simpleName =~ /fodf_metrics.*_afds/
                    return [it[0], "fodf_afds", it[1]]
                afdt: it[1].simpleName =~ /fodf_metrics.*_afdt/
                    return [it[0], "fodf_afdt", it[1]]
                nufo: it[1].simpleName =~ /fodf_metrics.*_nufo/
                    return [it[0], "fodf_nufo", it[1]]
                rgb: it[1].simpleName =~ /fodf_metrics.*_rgb/
                    return [it[0], "fodf_rgb", it[1]]
            }

        fodf_peaks_splitting_channel = t1_channel
            .join(fodf_peaks_channel)

        afd_to_shot_channel = fodf_metrics_splitting_channel.afd
            .groupTuple(by: [0, 1])
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        afds_to_shot_channel = fodf_metrics_splitting_channel.afds
            .groupTuple(by: [0, 1])
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        afdt_to_shot_channel = fodf_metrics_splitting_channel.afdt
            .groupTuple(by: [0, 1])
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        nufo_to_shot_channel = fodf_metrics_splitting_channel.nufo
            .groupTuple(by: [0, 1])
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        rgb_to_shot_channel = fodf_metrics_splitting_channel.rgb
            .groupTuple(by: [0, 1])
            //.join(masks_3t_to_ovelays_channel)
            .map{ it + [[], [], [], []] }

        peaks_to_shot_channel = fodf_peaks_splitting_channel
            .join(masks_3t_splitting_channel.wm_mask.map{ [it[0], [it[1]]] })
            .map{ [it[0], "fodf_peaks", it[1], [], [], [], it[2]] }

        screenshot_fodf_metrics(
            afd_to_shot_channel
                .mix(afds_to_shot_channel)
                .mix(afdt_to_shot_channel)
                .mix(nufo_to_shot_channel)
                .mix(rgb_to_shot_channel),
            [params.qc_screenshot_orientation, params.qc_extra_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [params.qc_masks_as_contours, true],
            params.qc_image_size,
            [false, false, false]
        )

        screenshot_fodf_peaks_coronal(
            peaks_to_shot_channel,
            ["coronal", 1.0],
            [params.qc_gif_slice_gap, false],
            [true, false],
            [params.qc_image_size[0].intdiv(3), params.qc_image_size[1]],
            [false, false, false]
        )

        screenshot_fodf_peaks_sagittal(
            peaks_to_shot_channel,
            ["sagittal", 1.0],
            [params.qc_gif_slice_gap, false],
            [true, false],
            [2 * params.qc_image_size[0].intdiv(3), params.qc_image_size[1]],
            [false, false, false]
        )

        peaks_sagittal_shots_channel = screenshot_fodf_peaks_sagittal.out.screenshots
            .transpose()
            .map{ it[0..1] + [(it[2].simpleName =~ /_slice_([0-9]+)/)[0][1], it[2]] }

        peaks_coronal_shots_channel = screenshot_fodf_peaks_coronal.out.screenshots
            .transpose()
            .map{ it[0..1] + [(it[2].simpleName =~ /_slice_([0-9]+)/)[0][1], it[2]] }

        peaks_to_stitch_channel = peaks_coronal_shots_channel
            .combine(peaks_sagittal_shots_channel, by: [0, 2])
            .map{ [it[0], "peaks_${it[1]}", [it[3], it[5]]] }

        stitch_peaks_screenshots(peaks_to_stitch_channel, false)

        peaks_stitch_to_gif_channel = stitch_peaks_screenshots.out.image
            .map{ [it[0], it[1].simpleName.tokenize('_')[0..-2].join("_"), it[1]] }
            .groupTuple(by: [0, 1])
            .map{ [it[0], it[1].tokenize('_')[1..-1].join("_"), it[2].sort(get_natural_sorter())] }

        compose_fodf_qc_gifs(
            screenshot_fodf_metrics.out.screenshots
                .map{ it[0..1] + [it[2].sort(get_natural_sorter())] }
                .mix(peaks_stitch_to_gif_channel),
            params.qc_image_size,
            [params.qc_dti_gif_framerate, 0, false, false]
        )

        generate_fodf_report(
            compose_fodf_qc_gifs.out.gif
                .collect{ it[1] }
                .map{ [it, []] },
            "fodf_metrics",
            [false, true]
        )
}


workflow tissue_segmentation_qc_wkf {
    take:
        t1_channel
        masks_3t_channel
        pvf_3t_channel
        safe_wm_channel
    main:
        masks_3t_to_shot_channel = t1_channel
            .join(masks_3t_channel)
            .map{ it[0..1] + [it[2..-1]] }
            .transpose()
            .branch{
                wm_mask: it[2].simpleName =~ /(?<!safe)_wm_mask/
                    return [it[0], "wm_mask"] + it[1..-1] + [[], [], []]
                gm_mask: it[2].simpleName =~ /_gm_mask/
                    return [it[0], "gm_mask"] + it[1..-1] + [[], [], []]
                csf_mask: it[2].simpleName =~ /_csf_mask/
                    return [it[0], "csf_mask"] + it[1..-1] + [[], [], []]
            }

        pvf_3t_to_shot_channel = t1_channel
            .join(pvf_3t_channel)
            .map{ it[0..1] + [it[2..-1]] }
            .transpose()
            .branch{
                wm_pvf: it[2].simpleName =~ /_wm_pvf/
                    return [it[0], "wm_pvf", it[1], [], it[2], [], []]
                gm_pvf: it[2].simpleName =~ /_gm_pvf/
                    return [it[0], "gm_pvf", it[1], [], it[2], [], []]
                csf_pvf: it[2].simpleName =~ /_csf_pvf/
                    return [it[0], "csf_pvf", it[1], [], it[2], [], []]
            }

        safe_wm_mask_to_shot_channel = t1_channel
            .join(safe_wm_channel)
            .map{ it[0..1] + [it[2..-1]] }
            .map{ [it[0], "safe_wm_mask"] + it[1..-1] + [[], [], []] }

        screenshot_tissues_masks(
            masks_3t_to_shot_channel.wm_mask
                .mix(masks_3t_to_shot_channel.gm_mask)
                .mix(masks_3t_to_shot_channel.csf_mask)
                .mix(safe_wm_mask_to_shot_channel),
            [params.qc_screenshot_orientation, params.qc_primary_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [params.qc_masks_as_contours, true],
            [params.qc_image_size[0].intdiv(3), params.qc_image_size[1]],
            [false, false, false]
        )

        screenshot_tissues_pvfs(
            pvf_3t_to_shot_channel.wm_pvf
                .mix(pvf_3t_to_shot_channel.gm_pvf)
                .mix(pvf_3t_to_shot_channel.csf_pvf),
            [params.qc_screenshot_orientation, params.qc_primary_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [params.qc_masks_as_contours, true],
            [params.qc_image_size[0].intdiv(3), params.qc_image_size[1]],
            [false, params.qc_pvf_cmap_name, params.qc_pvf_opacity]
        )

        tissues_shots_channel = screenshot_tissues_masks.out.screenshots
            .mix(screenshot_tissues_pvfs.out.screenshots)
            .transpose()
            .map{ it[0..1] + [(it[2].simpleName =~ /_slice_([0-9]+)/)[0][1], it[2]] }
            .branch{
                wm_mask: it[1] == "wm_mask"
                gm_mask: it[1] == "gm_mask"
                csf_mask: it[1] == "csf_mask"
                safe_wm_mask: it[1] == "safe_wm_mask"
                wm_pvf: it[1] == "wm_pvf"
                gm_pvf: it[1] == "gm_pvf"
                csf_pvf: it[1] == "csf_pvf"
            }

        mask_3t_to_stitch_channel = tissues_shots_channel.wm_mask
            .combine(tissues_shots_channel.gm_mask, by: [0, 2])
            .map{ [it[0], it[2], it[1]] + it[3..-1] }
            .combine(tissues_shots_channel.csf_mask, by: [0, 2])
            .map{ [it[0], "3t_masks_${it[1]}", [it[3], it[5], it[7]]] }

        pvf_3t_to_stitch_channel = tissues_shots_channel.wm_pvf
            .combine(tissues_shots_channel.gm_pvf, by: [0, 2])
            .map{ [it[0], it[2], it[1]] + it[3..-1] }
            .combine(tissues_shots_channel.csf_pvf, by: [0, 2])
            .map{ [it[0], "3t_pvfs_${it[1]}", [it[3], it[5], it[7]]] }

        safe_wm_mask_to_stitch_channel = tissues_shots_channel.safe_wm_mask
            .combine(tissues_shots_channel.wm_mask, by: [0, 2])
            .map{ [it[0], "safe_wm_mask_${it[1]}", [it[3], it[5]]] }

        stitch_tissues_screenshots(
            mask_3t_to_stitch_channel
                .mix(pvf_3t_to_stitch_channel)
                .mix(safe_wm_mask_to_stitch_channel),
            false
        )

        compose_tissues_qc_gifs(
            stitch_tissues_screenshots.out.image
                .map{ [it[0], it[1].simpleName.tokenize('_')[0..-2].join("_"), it[1]] }
                .groupTuple(by: [0, 1])
                .map{ [it[0], it[1].tokenize('_')[1..-1].join("_"), it[2].sort(get_natural_sorter())] },
            params.qc_image_size,
            [params.qc_masking_gif_framerate, 0, false, false]
        )

        tissues_gifs_channel = compose_tissues_qc_gifs.out.gif
            .branch{
                tissues_mask: it[1].simpleName.contains("3t_masks")
                tissues_pvfs: it[1].simpleName.contains("3t_pvfs")
                safe_wm_mask: true
            }

        generate_tissues_masks_report(
            tissues_gifs_channel.tissues_mask
                .collect{ it[1] }
                .map{ [it, []] },
            "tissues_masking",
            [false, true]
        )

        generate_tissues_pvf_report(
            tissues_gifs_channel.tissues_pvfs
                .collect{ it[1] }
                .map{ [it, []] },
            "tissues_fractions",
            [false, true]
        )

        generate_safe_wm_masks_report(
            tissues_gifs_channel.safe_wm_mask
                .collect{ it[1] }
                .map{ [it, []] },
            "safe_wm_masking",
            [false, true]
        )

}

workflow masks_qc_wkf {
    take:
        b0_channel
        t1_channel
        dwi_mask_channel
        t1_mask_channel
    main:

        average_to_3d_volume(b0_channel)

        screenshot_dwi_mask(
            average_to_3d_volume.out.image
                .join(dwi_mask_channel)
                .map{ [it[0], "dwi_mask"] + it[1..-1] + [[], [], []] },
            [params.qc_screenshot_orientation, params.qc_extra_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [params.qc_masks_as_contours, true],
            params.qc_image_size,
            [params.qc_b0_colormap_name, false, false]
        )

        screenshot_t1_mask(
            t1_channel
                .join(t1_mask_channel)
                .map{ [it[0], "t1_mask"] + it[1..-1] + [[], [], []] },
            [params.qc_screenshot_orientation, params.qc_extra_overlay_opacity],
            [params.qc_gif_slice_gap, false],
            [params.qc_masks_as_contours, true],
            params.qc_image_size,
            [params.qc_t1_colormap_name, false, false]
        )

        compose_mask_qc_gifs(
            screenshot_dwi_mask.out.screenshots
                .mix(screenshot_t1_mask.out.screenshots)
                .map{ it[0..1] + [it[2].sort(get_natural_sorter())] },
            params.qc_image_size,
            [params.qc_masking_gif_framerate, 0, false, false]
        )

        masks_gifs_channel = compose_mask_qc_gifs.out.gif
            .branch{
                dwi_mask: it[1].simpleName.contains("dwi_mask")
                t1_mask: it[1].simpleName.contains("t1_mask")
            }

        generate_dwi_masking_report(
            masks_gifs_channel.dwi_mask
                .collect{ it[1] }
                .map{ [it, []] },
            "dwi_masking",
            [false, true]
        )

        generate_t1_masking_report(
            masks_gifs_channel.t1_mask
                .collect{ it[1] }
                .map{ [it, []] },
            "t1_masking",
            [false, true]
        )

}

workflow dwi_qc_wkf {
    take:
        raw_dwi_channel
        processed_dwi_channel
        processed_screenshot_positions_channel
        screenshot_axis
    main:

        transform_bounding_box(
            processed_screenshot_positions_channel
                .join(processed_dwi_channel.map{ it[0..1] })
                .join(raw_dwi_channel.map{ it[0..1] }),
            false,
            ""
        )

        raw_slice_extract_id_channel = transform_bounding_box.out.box
            .map{ sid, bbox_json -> 
                def bbox = (new JsonSlurper()).parseText(bbox_json.getText())
                [sid, screenshot_axis, bbox.lower[screenshot_axes[screenshot_axis]]]
            }

        processed_slice_extract_id_channel = processed_screenshot_positions_channel
            .map{ sid, bbox_json -> 
                def bbox = (new JsonSlurper()).parseText(bbox_json.getText())
                [sid, screenshot_axis, bbox.lower[screenshot_axes[screenshot_axis]]]
            }

        extract_from_raw_volume(
            raw_dwi_channel
                .map{ it[0..1] }
                .join(raw_slice_extract_id_channel),
            false,
            ""
        )
        extract_from_processed_volume(
            processed_dwi_channel
                .map{ it[0..1] }
                .join(processed_slice_extract_id_channel),
            false,
            ""
        )

        raw_for_screenshot_channel = extract_from_raw_volume.out.image
            .map{ [it[0], "raw_dwi", it[1], [], [], [], []] }
        processed_for_screenshot_channel = extract_from_processed_volume.out.image
            .map{ [it[0], "processed_dwi", it[1], [], [], [], []] }

        screenshot_dwi_volumes(
            raw_for_screenshot_channel
                .mix(processed_for_screenshot_channel),
            [screenshot_axes[screenshot_axis], 0.7],
            [false, false],
            [params.qc_masks_as_contours, true],
            [params.qc_image_size[0].intdiv(2), params.qc_image_size[1]],
            [params.qc_dwi_colormap_name, false, false]
        )

        dwi_screenshots_channel = screenshot_dwi_volumes.out.screenshots
            .transpose()
            .map{ it[0..1] + [(it[2].simpleName =~ /_slice_([0-9]+)/)[0][1], it[2]] }
            .branch{
                raw_dwi: it[1] == "raw_dwi"
                processed_dwi: it[1] == "processed_dwi"
            }

        stitch_dwi_raw_vs_processed(
            dwi_screenshots_channel.raw_dwi
                .combine(dwi_screenshots_channel.processed_dwi, by: [0, 2])
                .map{ [it[0], "dwi_raw_vs_processed_${it[1]}", [it[3], it[5]]] },
            false
        )

        compose_dwi_qc_gifs(
            stitch_dwi_raw_vs_processed.out.image
                .map{ [it[0].tokenize('-')[0]] + it[1..-1] }
                .groupTuple()
                .map{ [it[0], it[1].sort(get_natural_sorter())] }
                .map{ [it[0], "dwi_raw_vs_processed", it[1]] },
            [params.qc_image_size[0], params.qc_image_size[1]],
            [params.qc_dwi_gif_framerate, 0, false, false]
        )

        generate_dwi_processing_report(
            compose_dwi_qc_gifs.out.gif
                .collect{ it[1] }
                .map{ [it, []] },
            "dwi_raw_vs_processed",
            [false, true]
        )

        raw_gradients_channel = raw_dwi_channel
            .map{ [it[0], "raw", it[2], it[3]] }

        processed_gradients_channel = processed_dwi_channel
            .map{ [it[0], "processed", it[2], it[3]] }

        screenshot_gradient_sampling_fsl(
            raw_gradients_channel
                .mix(processed_gradients_channel),
            params.qc_image_size[1],
            [params.qc_gradients_sphere, false],
            [params.qc_gradients_disable_sphere, params.qc_gradients_disable_supershell, params.qc_gradients_split_shells],
            [params.qc_gradients_sphere_opacity, params.qc_gradients_shells_same_color]
        )

        gradients_screenshots_channel = screenshot_gradient_sampling_fsl.out.screenshots
            .branch{
                raw: it[1] == "raw"
                processed: it[1] == "processed"
            }

        stitch_gradients_raw_vs_processed(
            gradients_screenshots_channel.raw
                .join(gradients_screenshots_channel.processed)
                .map{ [it[0], "gradients_raw_vs_processed", [it[2], it[4]]] },
            false
        )

        generate_gradients_processing_report(
            stitch_gradients_raw_vs_processed.out.image
                .collect{ it[1] }
                .map{ [it, []] },
            "gradients_raw_vs_processed",
            [false, true]
        )
}

/*workflow t1_qc_wkf {
    take:
    main:
}*/