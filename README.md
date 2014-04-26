Shadow Free Segmentation
========================

This software package provides a MATLAB implementation of the foreground background segmentation approach that robust to the effects of shadows described in the following paper:

A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in Still Images using Local Density Measure", ICCP 2014

Required 3rd part packages
==========================

There is a number of 3rd party software packages that are used by this software:

- My package of usefull MATLAB functions [[link](https://github.com/aecins/utilities)]
- Berkeley segmentation code v1.2.0 [[link](https://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/)]
- Maxflow v3.0.1 [[link](http://vision.csd.uwo.ca/code/)]
- Maxflow wrapper for MATLAB [[link](http://www.mathworks.com/matlabcentral/fileexchange/21310-maxflow)]
- EDISSON wrapper for MATLAB [[link](http://www.wisdom.weizmann.ac.il/~bagon/matlab.html)]
- libsvm v3.18 [[link](https://github.com/cjlin1/libsvm)]

    NOTE: I have renamed functions `svmtrain` and `svmpredict` to `libsvmtrain` and `libsvmpredict` in my installation to avoid confusion with the same functions available in MATLAB Statistics Toolbox.
    
- Local Binary Patterns (already included in this package) [[link](http://www.cse.oulu.fi/MVG/Downloads/LBPMatlab)]


Usage
=====

The package contains 4 important pieces of code. Each of them can be used standalone for their respective purpose:

- `fixation_segmentation` - fixation based foreground background segmentation
- `shadow_boundary_recognition` - shadow boundary classfier training
- `density_map` - image transformation that is robust to illumination changes
- `shadow_free_seg` - shadow free segmentation

Each package had a demo file included. To try them out install the required 3rd party packages, open the demo file (e.g. `shadow_free_seg_demo.m`) set the required paths to the isntalled packages and run the demo file.

All software was tested in Ubuntu 12.04, MATLAB R2014a. If you have any questions or comments feel free to contact me.
