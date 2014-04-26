%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo of the shadow free segmentation algorithm presented in the
% paper:
%   A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in
%   Still Images using Local Density Measure", ICCP 2014
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add paths

% My code
addpath(genpath('../../../utilities/'));                                    % Utility functions
addpath('../density_map/');                                                 % Density map
addpath('../fixation_segmentation/');                                       % Fixation segmentation
addpath(genpath('../shadow_boundary_recognition/'));                        % Feature extraction functions

% 3rd party
addpath('../../3rd_party/lbp/');                                            % Local Binary Patterns
otherCodeDirPath = '../../../../';                                     % Directory containing the rest of the 3rd party code
addpath(fullfile(otherCodeDirPath, 'maxflow'));                             % Maxflow
addpath(fullfile(otherCodeDirPath, '/segbench/lib/matlab/'));               % Berkeley edges and filterbank
addpath(fullfile(otherCodeDirPath, 'libsvm/matlab/'));                      % libsvm
addpath(fullfile(otherCodeDirPath, 'edison_wrapper/'));                     % Mean-shift

%% Clear and close
close all
clearAllPreserveBreakpoints

%% Get shadow recognition model and extract parameters
modelDirPath = '../shadow_boundary_recognition/models/';
modelFileName = '1_1_0_0_1_tex_texton128_malik_0x5_3_imagespec_density5_8.mat';

% Get scale
nameSplit = regexp(modelFileName(1:end-4),'_','split');
scales = str2double(nameSplit{end});

for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));
end

% Get texture mode
startIdx    = regexp(modelFileName, 'tex_');
endIdx      = regexp(modelFileName, '_');
endIdx      = endIdx(end)-1;
texMode     = {modelFileName(startIdx:endIdx)};

% Get feature use
featureUse.intRatio     = str2double(nameSplit(1));
featureUse.RGBratio     = str2double(nameSplit(2));
featureUse.chromAlign   = str2double(nameSplit(3));
featureUse.labRatio     = str2double(nameSplit(4));
featureUse.texture      = str2double(nameSplit(5));

%% Prepare image

% Read image
im              = imread('./input_image.png');                              % Read image
im              = im2double(im);
[yRes, xRes, ~] = size(im);

% Choose fixation point
fig_select_fixation = figure;
imshow(im);
title('Select a fixation point inside the object');
[fix_x, fix_y] = ginput(1);
fix_x = round(fix_x);
fix_y = round(fix_y);
close(fig_select_fixation);

% Compute Berkeley edges
% NOTE: here we are loading precomputed Berkeley edges for the sake of
% speed but if you have your own image just compute the edges by
% uncommenting the line below.
fprintf('Computing Berkeley edges...\n');
% load('input_image_grad.mat');                                               % Load Berkeley edges
[edgeGrad, ~]   = pbCGTG(im);                                               % Compute Berkeley edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running initial segmentation...\n');
[labelsOrig, unaryFG, unaryBG] = segment(im, edgeGrad, fix_x, fix_y);

% Visualize segmentation
segBoundary = bwmorph(labelsOrig==1, 'remove');
segOrig     = imRegionHighlight(im, segBoundary == 1, 'g');

figure;
imshow(segOrig);
title('Original segmentation without shadow removal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix shadows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Split the segmentation boundary
fprintf('Splitting the segmentation boundary...\n');

%% Extract features
fprintf('Extracting features...\n');
[boundaries, internalBoundaries] = segmentBoundary(labelsOrig, im);         % Extract internal boundaries and boundary segments
[features, reject] = getFeatures(im, boundaries, scales, texMode);          % Extract features

if ~isempty(reject)
    error(['For now we have no means of dealing with cases when some of the',... 
            'segmentation boundary is too close to the image boundary and thus',...
            'cannot be classified as shadow/non-shadow.']);
end

curFeatures = getFeatureVector(features, scaleNames, featureUse, texMode);  % Get feature vector

%% Classify shadow edges
fprintf('Classifying...\n');
load(fullfile(modelDirPath, modelFileName));                                % Load model    
[~, ~, probValues] = libsvmpredict(ones(numel(boundaries), 1), curFeatures, curModel, '-b 1'); % Classify

% Visualize classification  results
nbColors = 32;
cmap = jet(nbColors);
[h,probQuant] = histc(probValues(:,1), linspace(0,1+eps,nbColors));
imVis = im;
for j=1:numel(boundaries)

    tmpMask = false([yRes, xRes]);
    tmpMask(boundaries{j}.pixelIdx) = true;
    tmpMask = bwmorph(tmpMask, 'dilate');
    imVis = imRegionHighlight(imVis, tmpMask, cmap(probQuant(j), :));
end

figure;
imshow(imVis);
title('Shadow recognition results');

%% Modify segmentation energy function
fprintf('Modfying the segmentation energy function\n');
strongShadowThresh = 0.8;
strongNonShadowThresh = 0.2;
if any(probValues(:,1) > strongShadowThresh)


    % Get shadow probability mask
    shadowProbMask          = zeros([yRes, xRes]);
    segBnd                  = zeros([yRes, xRes]);
    strongShadowSegments    = zeros([yRes, xRes]);
    strongNonShadowSegments = zeros([yRes, xRes]);
    for j=1:numel(boundaries)
        shadowProbMask(boundaries{j}.pixelIdx) = probValues(j,1);
        segBnd(boundaries{j}.pixelIdx) = 1;

    end

    % Find strong shadow edges
    shadowEdges = im2bw(shadowProbMask, strongShadowThresh);

    % Find strong non-shadow edges
    nonShadowEdges = shadowProbMask;
    nonShadowEdges(~segBnd) = 1;
    nonShadowEdges = imcomplement(im2bw(nonShadowEdges, strongNonShadowThresh));        

    % Modify binary weights

    % Reduce shadow edge
    tmp = bwmorph(shadowEdges, 'dilate', 2);
    edgeGrad(tmp) = 0.0001;        

    % Increase all internal segmentation area and internal boundaries
    edgeGrad(labelsOrig) = max(edgeGrad(labelsOrig), 0.00000000001);
    edgeGrad(internalBoundaries) = max(edgeGrad(internalBoundaries), 0.1);

    % Modify unary weights
    [y, x] = find(nonShadowEdges);
    fix_x = [fix_x; x];
    fix_y = [fix_y; y];

    %% Resegment
    fprintf('Resegmenting...\n');
    [labelsMod,  unaryFG, unaryBG] = segment(im, edgeGrad, fix_x, fix_y);

else
    segMod = segOrig;
    labelsMod = labelsOrig;
end

fprintf('Done!\n');
% Visualize segmentation
segBoundary = bwmorph(labelsMod==1, 'remove');
segMod     = imRegionHighlight(im, segBoundary == 1, 'g');

figure;
imshow(segMod);
title('Shadow free segmentation');