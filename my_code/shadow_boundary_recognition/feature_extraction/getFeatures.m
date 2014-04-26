function [feat, reject] = getFeatures(im, boundaries, scales, texMode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute boundary features for input image
%
% Input:
%   im,                 input image
%   boundaries,         input boundaries on which features are computed
%   scales,             scales of the halfdisc for feature extraction
%   texMode,            texture description mode used
%
% Output:
%   feat,               data structure with features for current boundaries
%   reject,             indices of boudnaries that were rejected
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Feature collection parameters
nOrient = 16;
radius = scales*2;

% Filterbank parameters
fbParams.no = 8;
fbParams.ns = 1;
fbParams.sc = sqrt(2);

for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));
end

%% Preprocess
im = im2double(im);                                                        % Convert to double
imGray = rgb2gray(im);                                                     % Get grayscale image
[yRes, xRes, ~] = size(im);                                                % Get the size of the image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare edge pixel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filter out boudnaries that contain pixels that are too close to the image
%% boundary (closer than the radius of the halfdisc).

% Limits for the boundary pixels
leftX   = max(radius)+1;
rightX  = xRes - (max(radius));
topY    = max(radius)+1;
bottomY = yRes - (max(radius));

% List of boundaries that need to be filtered out
reject = [];

% Loop through all boundaries
for i=1:numel(boundaries)
    
    % Get boundary pixel coordinates
    pixelCrd = [];
    [pixelCrd(:,1), pixelCrd(:,2)] = ind2sub([yRes xRes], boundaries{i}.pixelIdx);
    boundaries{i}.pixelCrd = pixelCrd;
    
    % Check if any pixels are too close to image boundaries
    outliers =  pixelCrd(:,2) < leftX |...
                pixelCrd(:,2) > rightX |...
                pixelCrd(:,1) < topY |...
                pixelCrd(:,1) > bottomY;
  
    % Add current boundary to filter list
    if any(outliers)
        reject(end+1) = i;
%     else
%         boundaries{i}.pixelIdx = boundaries{i}.pixelIdx(~outliers);
%         boundaries{i}.pixelCrd = boundaries{i}.pixelCrd(~outliers, :);
%         boundaries{i}.orient   = boundaries{i}.orient(~outliers,:);
    end
end

boundaries(reject) = [];

%% Concatenate all boundaries
boundariesCat.pixelIdx = [];
boundariesCat.pixelCrd = [];
boundariesCat.orient = [];
boundariesCat.bndId = [];

for i=1:numel(boundaries)
    boundariesCat.pixelIdx = [boundariesCat.pixelIdx; boundaries{i}.pixelIdx];
    boundariesCat.pixelCrd = [boundariesCat.pixelCrd; boundaries{i}.pixelCrd];
    boundariesCat.orient = [boundariesCat.orient; boundaries{i}.orient];
    boundariesCat.bndId = [boundariesCat.bndId; i*ones(size(boundaries{i}.pixelIdx))];
end

%% Get the mask at different orientations
% NOTE: we want the ratio comparisons to be consistent between each other.
% That is why we create our masks such that first halfdisc is always on the
% side of the gradient and the other halfdisc is on the side of the
% negative gradient
[orientedMasks, orientations]= getOrientedMasks(nOrient, radius);

%% Get boundary pixel orientation angles
gX = boundariesCat.orient(:,1);
gY = boundariesCat.orient(:,2);

gTheta = atan2(gY, gX);

%% Find closest mask orientation for each boundary pixel
angleDiff = min(...
    mod(repmat(gTheta+pi/2, 1, nOrient) - repmat(orientations(:)', numel(boundariesCat.bndId), 1), 2*pi), ...
    mod(repmat(gTheta+pi/2, 1, nOrient) - repmat(orientations(:)'+pi, numel(boundariesCat.bndId), 1), 2*pi));
[~, boundariesCat.orientIdx] = min(angleDiff, [], 2);

% figure
% visualizeMasks(im, orientedMasks, boundariesCat.pixelIdx, boundariesCat.pixelCrd, boundariesCat.orientIdx, gX, gY);
% pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract pixel features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intensity ratio
[tmp.intRatio, brightDarkSide] = getIntensityRatioFeatures(im, scales, orientedMasks, boundariesCat.pixelCrd, boundariesCat.orientIdx);

%% Texture
for i = 1:numel(texMode)
    curTexMode = texMode{i};
        texture.(curTexMode) = getTextureFeatures(im, scales, orientedMasks, boundariesCat.pixelCrd, boundariesCat.orientIdx, curTexMode, fbParams);
end

%% RGB ratio
tmp.RGBratio = getRGBratioFeatures(im, scales, orientedMasks, boundariesCat.pixelCrd, boundariesCat.orientIdx, brightDarkSide);

%% Chromatic alignment
tmp.chromAlign = getChromaticAlignmentFeatures(im, scales, orientedMasks, boundariesCat.pixelCrd, boundariesCat.orientIdx, brightDarkSide);

%% LAB ratio
tmp.labRatio = getLABratioFeatures(im, scales, orientedMasks, boundariesCat.pixelCrd, boundariesCat.orientIdx, brightDarkSide);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge boundary features by takin a mean over the boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get scale
featureNames = fieldnames(tmp);

for i=1:numel(boundaries)
    for s=1:numel(scales)
    
        % Non texture features
        for curFeatName = featureNames'
            curFeat = tmp.(curFeatName{1}).(scaleNames{s})(boundariesCat.bndId == i,:);
            feat.(curFeatName{1}).(scaleNames{s})(i, :) = mean(curFeat,1);
        end

        % Texture features
        for newTexMode = fieldnames(texture)'
            curTexDist = texture.(newTexMode{1}).(scaleNames{s})(boundariesCat.bndId == i);
            feat.texture.(newTexMode{1}).(scaleNames{s})(i, 1) = mean(curTexDist);
        end          
    end
    
    feat.label(i, 1) = boundaries{i}.label;    
end