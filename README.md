# Expanded signal model for image reconstruction
This MATLAB code reconstructs MR images using expanded signal model image reconstruction (Wilm et al., 2011). Field measurements and B<sub>0</sub> nonuniformity map are incorporated into image reconstruction to provide images with minimal distortions and artifacts.

## Requirements
* GPU is required to use the code.
* ISMRMRD reader for MATLAB is needed to read data. It can be found here: [github.com/ismrmrd/ismrmrd](https://github.com/ismrmrd/ismrmrd)


## ISMRMRD format
Using ISMRMRD files that have raw k-space and trajectory data is the most convenient way to use the code. Also, some information in ISMRMRD header files are used during image reconstruction.

## How to use
Follow steps in ***example_cg_recon.m*** to learn how to use the code.

In order to efficiently use GPU, maximum variable size defined in ***cg_encoding_loops_gpu.m*** needs to be determined depending on the model of the GPU. The default value is tested for RTX 3090 and Tesla v100. For the best performance, tune this value until reaching the shortest image reconstruction duration, or maximum power drawn by the GPU.

## Limitations
* This implementation is for single-band excitation sequences
* Does NOT use gridding for image reconstruction
* Multifrequency interpolation (MFI) for correction of B<sub>0</sub> nonuniformity is NOT implemented.

**If you have any question regarding the code, please contact: [sajjad.feizollah@mail.mcgill.ca](mailto:sajjad.feizollah@mail.mcgill.ca)**
