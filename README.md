# Expanded signal model for image reconstruction
This MATLAB code reconstructs MR images using the expanded signal model image reconstruction [(Wilm et al., 2011)](https://doi.org/10.1002/mrm.22767). Field measurements and a static B<sub>0</sub> non-uniformity map are incorporated into the image reconstruction to provide images with minimal distortions and artifacts.

## Requirements
* A GPU is required to run the code.
* The MATLAB ISMRMRD reader is required to read the raw MR data: [github.com/ismrmrd/ismrmrd](https://github.com/ismrmrd/ismrmrd)


## ISMRMRD format
Using ISMRMRD files that have raw k-space and trajectory data is the most convenient way to use the code. Information in ISMRMRD header files is used during image reconstruction.

## How to use
Follow the steps in ***example_cg_recon.m*** to learn how to use the code.

The maximum size of a matrix defined in ***cg_encoding_loops_gpu.m*** needs to be determined depending on the model of the GPU. The default value is for RTX 3090 and Tesla v100. For best performance, tune this value until reaching the shortest image reconstruction duration, or maximum power drawn by the GPU.

## This version
*	is only for single-band excitation sequences;
* does NOT use gridding for image reconstruction;
* does NOT use multi-frequency interpolation (MFI) for correction of B<sub>0</sub> non-uniformity.

**If you have any question regarding the code, please contact: [sajjad.feizollah@mail.mcgill.ca](mailto:sajjad.feizollah@mail.mcgill.ca)**
