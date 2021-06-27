MAGIC-MONKEY FLOW
=================

Magic-Monkey is a diffusion MRI processing pipeline configured to handle high 
angular and spatial resolution data. It's pre-processing workflow contains a 
collection of steps focussed at increasing SNR and CNR of noisy diffusion data 
and improving the alignment of longitudinal acquisitions, using stat-of-the-art 
algorithms such as MP-PCA denoising, Eddy and Topup, and libraries like Scipy, 
Dipy, Mrtrix, FSL and ANTs.

Requirements
------------

To run the pipeline, the following tools must be installed :

- Nextflow (we recommend installing the last available release of version 21)
- Singularity 3.7.1 or higher, or Docker

Dependencies
------------

All dependencies required to run the pipeline have been packaged in a singularity 
image, as well as a docker image. They can be found at :

- Singularity : [Singularity Cloud](https://cloud.sylabs.io/library/avcaron/default/magic-monkey)
  - Pull command : `singularity pull library://avcaron/default/magic-monkey:v2.0`

- Docker : [Docker Hub](https://hub.docker.com/r/avcaron/magic-monkey)
  - Pull command : `docker pull avcaron/magic-monkey:latest`

Data input format
-----------------

A digestible subject must be composed of a minimum of the following images :

- A DWI image and its b-values (.bval) and b-vectors (.bvec) files, in the FSL format
- A T1 anatomical image
- A T1 or DWI mask (or both)

In addition to those, a subject can contain :

- A reverse phase acquired b0 or DWI image

**All images must be in the compressed Nifti format (.nii.gz)**

A single subject is identified by a **unique key**, prefixing the name of all images 
associated to it. This **key** must be composed of the directory names leading to the files, 
separated by a underscore symbol **\_**. Thus, all images must abide to the following naming convention :

- \<subject folder\>
    - \<key\>_dwi.nii.gz
    - \<key\>_dwi.bval
    - \<key\>_dwi.bvec
    - \<key\>_dwi.json *(optional)*
    - \<key\>_t1.nii.gz
    - \<key\>_t1_mask.nii.gz *(optional, if DWI mask is supplied)*
    - \<key\>_dwi_mask.nii.gz *(optional, if T1 mask is supplied)*
    - \<key\>_rev.nii.gz *(optional)*
    - \<key\>_rev.bval *(optional, unrequired if reverse acquisition is only a b0)*
    - \<key\>_rev.bvec *(optional, unrequired if reverse acquisition is only a b0)*
    - \<key\>_rev.json *(optional)*

Examples of name formatting :

- For a subject in a directory named **sub1**, all files associated to it must have the **key**
  prefix **sub1**. A diffusion image should thus be named **sub1_dwi.nii.gz**.
- For a session in a subdirectory **ses1** of the main directory of a subject **sub1**, the key would 
  be **sub1_ses1**. A diffusion image would then be named **sub1_ses1_dwi.nii.gz**.

Content of diffusion acquisition specification file (.json) :

The json file included with a dwi or rev image describes the parameters related to their specific acquisition. They 
are listed as a json object, in the fashion displayed below :

{
  "direction" : <phase encoding direction>,
  "slice_direction" : <slicing direction>,
  "readout" : <readout>,
  "interleaved" : <false|true>,
  "multiband_factor" : <multiband factor>
}

Optionally, one can omit the last two items and replace them by the parameter "slice_indexes", which is a list of list containing the indexes of
the different slices acquired in order. Two examples :

- An image with interleaved slices, multiband factor of 1, 10 slices : \[\[0\], \[2\], \[4\], \[6\], \[8\], \[1\], \[3\], \[5\], \[7\], \[9\]\]
- An image with interleaved slices, multiband factor of 2, 10 slices : \[\[0, 5\], \[2, 7\], \[4, 9\], \[1, 6\], \[3, 8\]\]

Running the pipeline
--------------------

To run the pipeline, we recommend creating a directory per run that contains log and 
cache of execution from Nextflow. In that directory, then run the following command

`nextflow run -resume -w cache -dsl2 <root_of_magic-monkey_nextflow_library>/main.nf --data_root <root_of_datasets> --resampling_resolution <final_image_resolution>`

Additional parameters to this command can be supplied. All parameters prefixed with a 
single hyphen ( - ) must be placed before the pipeline name, and with double hypen 
( -- ) after. Those parameters are :

- Running with docker image : `-with-docker avcaron/magic-monkey:latest`
- Running with singularity image : `-with-singularity <singularity_image>`
- Changing output directory : `--output_root <output_folder>`

If a **nextflow.config** file is present in the execution directory, it will be used 
as configuration for the run. A user can then modify the effective configuration for 
the run and turn on and off the execution of different processes.

Some configuration capabilities have been extracted to separate configuration files. 
We refer the user to the header of the **nextflow.config** file for more information. 
Those configuration files can be found in the *.config* directory.