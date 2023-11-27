#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def prepare_screenshot_parameters(
    orientation = "axial",
    overlay_opacity = 0.7,
    slice_gap = 5,
    slice_ids = false,
    overlay_as_contours = true,
    annotate = true,
    image_size = [4000, 4000],
    volume_colormap_name = false,
    labelmap_colormap_name = false,
    labelmap_opacity = 0.7,
    gif_framerate = 0.6,
    gif_n_loops = -1,
    gif_last_frame_dispose = false,
    gif_background_color = false
) {
    return [
        "orientation": orientation,
        "overlay_opacity": overlay_opacity,
        "slice_gap": slice_gap,
        "slice_ids": slice_ids,
        "overlay_as_contours": overlay_as_contours,
        "annotate": annotate,
        "image_size": image_size,
        "volume_colormap_name": volume_colormap_name,
        "labelmap_colormap_name": labelmap_colormap_name,
        "labelmap_opacity": labelmap_opacity,
        "gif_framerate": gif_framerate,
        "gif_n_loops": gif_n_loops,
        "gif_last_frame_dispose": gif_last_frame_dispose,
        "gif_background_color": gif_background_color
    ]
}


workflow stitching_wkf {
    take:
        images_to_screenshot_channel
        screenshot_parameters
        data_to_stitcher_closure
    main:
        stitch_screenshotter(
            images_to_screenshot_channel,
            [screenshot_parameters.orientation, screenshot_parameters.overlay_opacity],
            [screenshot_parameters.slice_gap, screenshot_parameters.slice_ids],
            [screenshot_parameters.overlay_as_contours, screenshot_parameters.annotate],
            screenshot_parameters.image_size,
            screenshot_parameters.colormap_names + [screenshot_parameters.labelmap_opacity]
        )

        images_to_stitch_channel = data_to_stitcher_closure(
            stitch_screenshotter.out.screenshots
        )

        stitch_stitcher(images_to_stitch_channel, false)

    emit:
        images = stitch_tissues_masks.out.image
        screenshots = stitch_screenshotter.out.screenshots
}


workflow gif_composing_wkf {
    take:
        images_to_screenshot_channel
        screenshot_parameters
        data_to_stitcher_closure
    main:
        screenshot_split_channel = images_to_screenshot_channel
            .groupTuple()
            .branch{
                requires_stitching: it[1].size() > 1
                    return it[0]
                direct_screenshotting: true
                    return it[0]
            }

        stitching_wkf(
            screenshot_split_channel.requires_stitching
                .join(images_to_screenshot_channel),
            screenshot_parameters,
            data_to_stitcher_closure
        )

        composing_screenshotter(
            screenshot_split_channel.direct_screenshotting
                .join(images_to_screenshot_channel),
            [screenshot_parameters.orientation, screenshot_parameters.overlay_opacity],
            [screenshot_parameters.slice_gap, screenshot_parameters.slice_ids],
            [screenshot_parameters.overlay_as_contours, screenshot_parameters.annotate],
            screenshot_parameters.image_size,
            screenshot_parameters.colormap_names + [screenshot_parameters.labelmap_opacity]
        )

        screenshots_channel = stitching_wkf.out.image
            .map{ [it[0], it[1].simpleName.tokenize('_')[0..-2].join("_"), it[1]] }
            .mix(composing_screenshotter.out.screenshots)
            .groupTuple(by: [0, 1])
            .map{ [it[0], it[2].sort(get_natural_sorter())] }
            .map{ [it[0], it[1].simpleName.tokenize('_')[1..-1].join("_"), it[1]] }

        compose_tissues_masks_qc_gifs(
            screenshots_channel,
            screenshot_parameters.image_size,
            [
                screenshot_parameters.gif_framerate,
                screenshot_parameters.gif_n_loops,
                screenshot_parameters.gif_last_frame_dispose,
                screenshot_parameters.gif_background_color
            ]
        )

    emit:
        gifs = compose_tissues_masks_qc_gifs.out.gif
        screenshots = screenshots_channel
}
