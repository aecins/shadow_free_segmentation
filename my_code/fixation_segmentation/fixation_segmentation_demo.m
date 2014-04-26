%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo of the fixation based foreground background segmentation
% algorithm presented in:
%   A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in
%   Still Images using Local Density Measure", ICCP 2014
%
% The original fixation based segmentation algorithm is by Ajay Misra:
%   Mishra, Ajay, Yiannis Aloimonos, and Cornelia Fermuller. "Active
%   segmentation for robotics." Intelligent Robots and Systems, 2009. IROS 
%   2009. IEEE/RSJ International Conference on. IEEE, 2009.
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths

% My code
addpath(genpath('../../../utilities/'));                                 % Utility functions

% 3rd party
otherCodeDirPath = '../../../../';                                          % Directory containing the rest of the 3rd party code
addpath(fullfile(otherCodeDirPath, 'maxflow'));                             % Maxflow
addpath(fullfile(otherCodeDirPath, 'segbench/lib/matlab/'));                % Berkeley edges and filterbank

%% Clear and close
close all
clearAllPreserveBreakpoints

%% Prepare image

% Read image
imName = 'dog_n_elephant';
im = imread(['Images/', imName, '.jpg']);
im = im2double(im);

% Choose fixation point
fig_select_fixation = figure;
imshow(im);
title('Select a fixation point inside the object');
[fix_x, fix_y] = ginput(1);
fix_x = round(fix_x);
fix_y = round(fix_y);
close(fig_select_fixation);

%% Segment
% [edgeGrad, ~]   = pbCGTG(im);                                               % Compute Berkeley edges
load(['Images/', imName '_grad.mat']);                                      % Use precomputed Berkeley edges
[labels, ~, ~] = segment(im, edgeGrad, fix_x, fix_y);                       % Segment

%% Show segmentation results
segBoundary = bwmorph(labels==1, 'remove');
segBoundary = bwmorph(segBoundary, 'dilate');
segOrig     = imRegionHighlight(im, segBoundary == 1, 'g');

figure;
imshow(segOrig);
title('Segmented image');