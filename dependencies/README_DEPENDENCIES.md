# Dependencies

This folder contains third-party scripts and libraries required by the medical imaging pipeline, as well as standalone components like Deep Learning architectures.

## Third-Party Scripts

- **NIfTI Tools**: The NIfTI scripts included in this directory belong to their respective authors and retain their original open-source licenses. They are used to read and write imaging arrays.
- **IVIM Scripts**: The IVIM (Intravoxel Incoherent Motion) scripts are the property of their respective authors and are distributed under their original open-source licenses. They handle Bayesian and segmented fitting approaches.
- **Deep Learning Components**:
  - `apply_dncnn_symmetric.m`: Handles the application of symmetric DnCNN image spatial denoising.
- **Miscellaneous Medical Imaging Tools**:
  - `dvh.m` and `sample_rtdose_on_image.m` are utilized for sampling dose on structures and producing fractional Dose-Volume Histograms.

Please refer to the individual files or subdirectories for specific license details, copyright notices, and usage restrictions regarding these third-party dependencies.
