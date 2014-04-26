%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo of the image density map transformation presented in:
%   A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in
%   Still Images using Local Density Measure", ICCP 2014
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear and close
close all
clear all

%% Read image
im = imread('./im_gray.png');
im = im2double(im);

%% Compute density
imDensity = density(im, 2:5);

%% Show results
imDensityVis = histeq(mat2gray(imDensity),20);             % Rescale density to [0, 1] range and contrast normalize
figure;
imshow(im);
title('Original grayscale image');

figure;
imshow(imDensityVis);
title('Corresponding density map');