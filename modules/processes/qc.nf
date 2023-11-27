#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output_dir = false

include {
    is_path_list
} from "../functions.nf"


def unpack_optional ( files ) {
    if ( files instanceof Path ) {
        if ( files.empty() ) return []
        if ( is_path_list(files) ) {
            return files
        } else {
            return [files]
        }
    }
    else if ( files instanceof ArrayList) {
        if ( files.isEmpty() ) return []
    }
    return files
}


process screenshot_3d_nifti {

input:
    tuple val(sid), val(directory_name), path(image), path(overlays, stageAs: 'ovl*.nii.gz') /* optional : can be [] */, path(labelmap, stageAs: 'lbl*.nii.gz') /* optional : can be [] */, path(transparency, stageAs: 'trans*.nii.gz') /* optional : can be [] */, path(peaks, stageAs: 'peaks*.nii.gz') /* optional : can be [] */
    tuple val(axis), val(overlay_alpha) // required
    tuple val(slice_gap), val(slice_ids) // optional, values can be false
    tuple val(masks_as_contours), val(annotate) // optional, values can be false
    tuple val(window_width), val(window_height) // optional, values can be false
    tuple val(volume_cmap), val(labelmap_cmap), val(labelmap_alpha) // optional, values can be false

output:
    tuple val(sid), val(directory_name), path("${directory_name}/${sid}_${axis}*.png"), emit: screenshots

script:

def axes = ["axial": 2, "coronal": 1, "sagittal": 0]
def overlays = unpack_optional(overlays)
def labelmap = unpack_optional(labelmap)
def peaks = unpack_optional(peaks)
def transparency = unpack_optional(transparency)
def optional_args = []
if ( overlays ) optional_args += ["--in_masks $overlays"]
if ( labelmap ) optional_args += ["--in_labelmap $labelmap"]
if ( transparency ) optional_args += ["--in_transparency_mask $transparency"]
if ( peaks ) optional_args += ["--in_peaks $peaks"]
if ( masks_as_contours ) optional_args += ["--masks_as_contours"]
if ( annotate ) optional_args += ["--display_slice_number", "--display_lr"]
if ( volume_cmap ) optional_args += ["--volume_cmap_name $volume_cmap"]
if ( labelmap_cmap ) optional_args += ["--labelmap_cmap_name $labelmap_cmap"]
if ( labelmap_alpha ) optional_args += ["--labelmap_alpha $labelmap_alpha"]
if ( window_width && window_height ) optional_args += ["--win_dims $window_width $window_height"]
if ( slice_gap && slice_ids ) error "Cannot specify both slice_gap and slice_ids for $task.name"
if ( slice_gap ) slice_ids = "\$(mrinfo -size ${image[0]} | cut -d' ' -f${axes[axis] + 1} --output-delimiter=NUL | tr -d '\n' | xargs printf '%s-1\n' | bc | tr -d '\n' | xargs -0 seq -s ' ' 0 $slice_gap)"
if ( slice_ids ) optional_args += ["--slice_ids $slice_ids"]

"""
mkdir -p $directory_name
scil_screenshot_volume.py $image $directory_name/${sid}_${axis}.png \
    ${optional_args.join(' \\\n    ')} \
    --axis_name $axis \
    --masks_alpha $overlay_alpha
"""

}

process screenshot_mosaic_nifti {

input:
    tuple val(sid), val(output_prefix), path(image), path(transparency), path(overlays, stageAs: 'ovl*.nii.gz') /* optional : can be [] */, path(labelmap, stageAs: 'lbl*.nii.gz') /* optional : can be [] */
    tuple val(mosaic_rows), val(mosaic_cols), val(axis), val(overlay_alpha) // required
    tuple val(slice_gap), val(slice_ids) // mutually exclusive, either values can be false
    tuple val(width_overlap), val(height_overlap) // optional, values can be false
    tuple val(masks_as_contours), val(annotate) // optional, values can be false
    tuple val(window_width), val(window_height) // optional, values can be false
    tuple val(volume_cmap), val(labelmap_cmap), val(labelmap_alpha) // optional, values can be false

output:
    tuple val(sid), path("${sid}_${output_prefix}_${axis}_mosaic.png"), emit: mosaic

script:

def axes = ["axial": 2, "coronal": 1, "sagittal": 0]
def overlays = unpack_optional(overlays)
def labelmap = unpack_optional(labelmap)
def optional_args = []
if ( overlays ) optional_args += ["--in_masks $overlays"]
if ( labelmap ) optional_args += ["--in_labelmap $labelmap"]
if ( masks_as_contours ) optional_args += ["--masks_as_contours"]
if ( annotate ) optional_args += ["--display_slice_number", "--display_lr"]
if ( volume_cmap ) optional_args += ["--volume_cmap_name $volume_cmap"]
if ( labelmap_cmap ) optional_args += ["--labelmap_cmap_name $labelmap_cmap"]
if ( labelmap_alpha ) optional_args += ["--labelmap_alpha $labelmap_alpha"]
if ( window_width && window_height ) optional_args += ["--win_dims $window_width $window_height"]
if ( width_overlap && height_overlap ) optional_args += ["--overlap_factor $width_overlap $height_overlap"]

if ( slice_gap && slice_ids ) error "Cannot specify both slice_gap and slice_ids for $task.name"
if ( slice_gap ) slice_ids = "\$(mrinfo -size $image | cut -d' ' -f${axes[axis] + 1} | tr -d '\n' | xargs printf '%s-1\n' | bc | tr -d '\n' | xargs seq -s ' ' 0 $slice_gap)"

"""
scil_screenshot_volume_mosaic_overlap.py $image ${sid}_${output_prefix}_${axis}_mosaic.png $transparency \
    $slice_ids $mosaic_rows $mosaic_cols \
    ${optional_args.join(' \\\n    ')} \
    --axis_name $axis \
    --masks_alpha $overlay_alpha
"""

}

process screenshot_gradient_sampling_fsl {

input:
    tuple val(sid), val(output_prefix), path(bval), path(bvec)
    val(resolution) // required
    tuple val(sphere), val(disable_symmetry) // optional, values can be false
    tuple val(disable_sphere), val(disable_supershell), val(individual_shells) // optional, values can be false
    tuple val(opacity), val(uniform_color) // optional, values can be false

output:
    tuple val(sid), val(output_prefix), path("${sid}_${output_prefix}_gradient_sampling.png"), emit: screenshots

script:

def optional_args = []
if ( sphere ) optional_args += ["--dipy_sphere $sphere"]
if ( disable_symmetry ) optional_args += ["--dis-sym"]
if ( disable_sphere ) optional_args += ["--dis-sphere"]
if ( disable_supershell ) optional_args += ["--dis-proj"]
if ( individual_shells ) optional_args += ["--plot_shells"]
if ( uniform_color ) optional_args += ["--same-color"] 
if ( opacity ) optional_args += ["--opacity $opacity"]

"""
scil_visualize_gradients.py \
    ${optional_args.join(' \\\n    ')} \
    --in_gradient_scheme $bval $bvec \
    --out_basename ${sid}_${output_prefix}_gradient_sampling \
    --res $resolution
"""

}

process compose_gif_from_images {

input:
    tuple val(sid), val(output_prefix), path(images, stageAs: 'img*.png')
    tuple val(width), val(height) // optional, values can be false
    tuple val(framerate), val(loop), val(dispose), val(bg_color) // optional, values can be false

output:
    tuple val(sid), path("${sid}_${output_prefix}.gif"), emit: gif

script:

def optional_args = []
if ( framerate ) optional_args += ["-delay \$(echo \"60.0/${framerate}\" | bc -l)"]
if ( loop ) optional_args += ["-loop $loop"]
if ( dispose ) optional_args += ["-dispose $dispose"]
if ( width && height ) optional_args += ["-size ${width}x${height}"]
if ( bg_color ) optional_args += ["$bg_color"]

"""
convert ${optional_args.join(' \\\n    ')} \
    \$(echo "$images" | tr ' ' '\n' | sort -V | xargs) \
    ${sid}_${output_prefix}.gif
"""

}

process stitch_images {

input:
    tuple val(sid), val(output_prefix), path(images, stageAs: 'img*.nii.gz')
    val(append_bottom)
output:
    tuple val(sid), path("${sid}_${output_prefix}.*"), emit: image

script:

def append_method = append_bottom ? "-append" : "+append"

"""
convert $images $append_method ${sid}_${output_prefix}.png
"""

}

process generate_report_from_screenshots {

scratch true

publishDir "${params.output_dir}/${output_name}", mode: 'copy'

input:
    tuple path(gifs, stageAs: "images/*"), path(stats) /* optional : can be [] */
    val(output_name)
    tuple val(symlink_files), val(use_cdn) // optional, values can be false

output:
    path("${output_name}.html"), emit: report
    tuple path("data/"), path("libs/"), emit: report_data

script:

def stats = unpack_optional(stats)
def optional_args = []
if ( stats ) optional_args += ["--stats $stats"]
if ( symlink_files ) optional_args += ["--sym_link"]
if ( use_cdn ) optional_args += ["--online"]

"""
mkdir -p screenshots
for img in images/*
do
    sid="\$(basename \$img | cut -d"." -f1 | cut -d"_" -f1)"
    ext="\$(basename \$img | cut -d"." -f2)"
    tab="\$(basename \$img | cut -d"." -f1 | cut --complement -d"_" -f1)"
    tab="\${tab// /_}"
    mkdir -p screenshots/\$tab/
    cp \$img screenshots/\$tab/\${sid}_\$tab.\$ext
done

dmriqc_from_screenshot.py ${output_name}.html \
    ${optional_args.join(' \\\n    ')} \
    --data screenshots/*
"""

}
