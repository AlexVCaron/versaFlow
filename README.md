[![Release build](https://github.com/AlexVCaron/mrHARDI/actions/workflows/build-release.yml/badge.svg)](https://github.com/AlexVCaron/mrHARDI/actions/workflows/build-release.yml)
[![Docker containers](https://img.shields.io/badge/Docker%20images-dockerhub-blue?style=plaflat&logo=docker&labelColor=2e343b)](https://hub.docker.com/repository/docker/avcaron/versa)
[![DOI](https://zenodo.org/badge/430873937.svg)](https://zenodo.org/badge/latestdoi/430873937)

# versaFlow

versaFlow is a diffusion MRI processing pipeline configured to handle high 
angular and spatial resolution data. By default it is configured to process 
Maccaca Mulatta brain images; profiles to handle other types of primates, as 
well as the human brain are being developped. The pipeline's pre-processing 
workflow contains a collection of steps focussed at increasing SNR and CNR of 
noisy diffusion data and improving the alignment of longitudinal acquisitions, 
using state-of-the-art algorithms such as MP-PCA denoising, Eddy and Topup, and 
libraries like Scipy, Dipy, Mrtrix, FSL and ANTs.

If you use this tool for your research, **please cite the following**

```
Valcourt Caron A., Shmuel A., Hao Z., Descoteaux M.,
“versaFlow: a versatile pipeline for resolution adapted diffusion MRI processing 
and its application to studying the variability of the PRIME-DE database”,
Frontiers in Neuroinformatics, 10.3389/fninf.2023.1191200, doi.org/10.3389/fninf.2023.1191200.
```

## Requirements

To run the pipeline, the following tools must be installed :

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (we recommend installing the last available release of version 21)
- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) 3.7.1 or higher, or [Docker](https://docs.docker.com/engine/install/)

## Dependencies

All dependencies required to run the pipeline have been packaged in singularity 
and docker images. The `latest` versions come pre-packaged with the CUDA runtime 
and require a Nvidia GPU to execute. For usage on a machine without a Nvidia GPU, 
use the images tagged `nogpu`.

- Docker : [Docker Hub](https://hub.docker.com/r/avcaron/versa)
  - `docker pull avcaron/versa:latest`
- Singularity :
  - Singularity images are no longer produced in house. To build your singularity, 
    use a docker tag and the following command : `singularity build <image.sif> docker://avcaron/versa:<tag>`.
  - Old versions of the singularity images can still be access through the ORAS 
    container registry :
    - Nvidia GPU  : `singularity pull oras://mrhardi.azurecr.io/mrHARDI/mrhardi:latest`
    - Without GPU : `singularity pull oras://mrhardi.azurecr.io/mrHARDI/mrhardi:nogpu`

# Data input format

**\*\* All images must be in the compressed Nifti format (.nii.gz) \*\***

A digestible subject must be composed of a minimum of the following images :

- A DWI image and its b-values (.bval) and b-vectors (.bvec) files, in the FSL format
- A T1 anatomical image
- A T1 or DWI mask (or both)

In addition to those, a subject can contain :

- A reverse phase acquired b0 or DWI image

A single subject is identified by a **unique key**, prefixing the name of all images 
associated to it. This [**key**](#) must be composed of the directory names leading to the 
files, separated by a underscore symbol `"_"`. Thus, all images must abide to the 
following naming convention :

- \<subject folder\>
    - ...
      - \<key\>_dwi.nii.gz
      - \<key\>_dwi.bval
      - \<key\>_dwi.bvec
      - \<key\>_dwi.json *(see [here](#alternative-specification-of-the-diffusion-acquisition-parameters))*
      - \<key\>_t1.nii.gz
      - \<key\>_t1_mask.nii.gz *(optional, if DWI mask is supplied)*
      - \<key\>_dwi_mask.nii.gz *(optional, if T1 mask is supplied)*
      - \<key\>_rev.nii.gz *(optional)*
      - \<key\>_rev.bval *(optional, unrequired if reverse acquisition is only a b0)*
      - \<key\>_rev.bvec *(optional, unrequired if reverse acquisition is only a b0)*
      - \<key\>_rev.json *(see [here](#alternative-specification-of-the-diffusion-acquisition-parameters))*
      - \<key\>_wm_pvf.nii.gz *(optional, will also be used to generate tissues mask)*
      - \<key\>_gm_pvf.nii.gz *(optional, will also be used to generate tissues mask)*
      - \<key\>_csf_pvf.nii.gz *(optional, will also be used to generate tissues mask)*

____

## Examples of name formatting :


- For a subject in a directory named **sub1**, all files associated to it must have 
  the **key** prefix **sub1**. A diffusion image should thus be named **sub1_dwi.nii.gz**.
- For a session in a subdirectory **ses1** of the main directory of a subject **sub1**, 
  the key would be **sub1_ses1**. A diffusion image would then be named **sub1_ses1_dwi.nii.gz**.

___

## Content of diffusion acquisition specification file (.json) :

The json file included with a dwi or rev image describes the parameters related to 
their specific acquisition. They are listed as a json object, in the fashion displayed 
below :

```
{
  "direction" : <phase encoding direction>,
  "slice_direction" : <slicing direction>,
  "readout" : <readout>,
  "interleaved" : <false|true>,
  "multiband_factor" : <multiband factor>
}
```

Optionally, one can replace the "interleaved" and "multiband_factor" fields by the 
field "slice_indexes", a list of list of indexes defining the sequence of acquisition 
of the image's slices. Two examples :

- An image with 10 slices, acquired interleaved, multiband factor of 1 :
  ```
  [[0], [2], [4], [6], [8], [1], [3], [5], [7], [9]]
  ```
- An image with interleaved slices, multiband factor of 2, 10 slices :
  ```
  [[0, 5], [2, 7], [4, 9], [1, 6], [3, 8]]
  ```

It is also possible instead of providing `.json` files to specify those parameters 
using default values passed through command line when calling the pipeline. For more 
information, refer to [this section](#alternative-specification-of-the-diffusion-acquisition-parameters) 
below.

The best way to extract the information to fill the .json configuration is to use 
`dcm2niix` with the output option `--json` when unpacking the base DICOM images. 
It will produce a json file of metadata.

As of now, this json is not directly consumable by the pipeline and must be parsed 
and transformed by the user to the structure presented above. This will be 
automatically handled in a future version of the pipeline.

# Running the pipeline

To run the pipeline, we recommend creating a directory per run that contains log and 
cache of execution from Nextflow. In that directory, then run the following command

```
nextflow run \
    -resume \
    -w cache \
    -dsl2 \
    <root of versaFlow>/main.nf \
    --data_root <data root>
```

Additional parameters to this command can be supplied, such as :

- Running with docker image : `-with-docker <docker image>`
- Running with singularity image : `-with-singularity <singularity image>`
- Changing output directory : `--output_root <output folder>`
- Display help and usage : `--help`

___

## Alternative specification of the diffusion acquisition parameters

By construction, the `.json` files provided with the dwi volumes are **required**. 
However, it is possible to define **default values** for all the fields in the `.json` 
specification presented [above](#content-of-diffusion-acquisition-specification-file-json), 
either via **command line** or through the `Nextflow.config` file, in which case 
the `.json` files will be **optional**. The command line parameters are :

```
--default_readout <readout value>
--default_multiband_factor <multiband factor>
--default_slicing_direction <slicing direction>
--default_phase_direction <phase direction>
--default_is_interleaved
```

___

## Nextflow.config file

If a **nextflow.config** file is present in the execution directory, it will be used 
as configuration for the run. A user can then modify the effective configuration for 
the run and turn on and off the execution of different processes.

Some configuration capabilities have been extracted to separate configuration files. 
We refer the user to the header of the **nextflow.config** file for more information. 
Those configuration files can be found in the *.config* directory.

___
