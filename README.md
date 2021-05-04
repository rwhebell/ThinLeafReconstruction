# ThinLeafReconstruction
MATLAB code for the article "Implicit reconstructions of thin leaf surfaces from large, noisy point clouds".

## Citing this work
If you use this code in your work, please cite the following paper:

Whebell, R. M., Moroney, T. J., Turner, I. W., Pethiyagoda, R., & McCue, S. W. (2020). Implicit reconstructions of thin leaf surfaces from large, noisy point clouds. _arXiv preprint arXiv:2009.10286_.

## Dependencies
This source code is dependent on the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox,
- Partial Differential Equation Toolbox, and
- Computer Vision Toolbox.

## Parameters
The main parameters of interest and their effects are:
- Downsampling grid size: `downsampleParam`. Higher values will reduce runtime but also reduce reconstruction resolution.
- Smoothing parameter: `rho`. A higher value will yield a smoother (as defined by the thin plate penalty) surface reconstruction.
- Sampling resolution: `Ngrid` for marching cubes, or `Hmax` for marching tetrahedra. 
    - A higher value for `Ngrid` will result in a finer grid for marching cubes isosurfacing. 
    - A lower value for `Hmax` will result in a finer tetrahedral mesh for marching tetrahedra isosurfacing.
