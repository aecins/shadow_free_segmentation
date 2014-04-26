function [intensityRatioFeat, brightDarkSide] = getIntensityRatioFeatures(im, scales, orientedMasks, edgePixelCrd, edgePixelOrientIdx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract intensity ratio features.
% 
% Input:
%   im,                 input image (colour or grayscale)
%   sacles,             halfdisc scales
%   orientedMasks,      halfdisk masks
%   edgePixelCrd,       coordinates of boundary pixels
%   edgePixelOrientIdx, orientation indeces of boundary pixels
%
% Output:
%   intensityRatioFeat, intensity ratio feature datastructure
%   brightDarkSide,     datastructure indicating which of the two halfdisks
%                       is dark and which is bright
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numScales = numel(scales);
numEdgePixels = size(edgePixelCrd, 1);
scaleNames = cell(numScales, 1);

for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));
    intensityRatioFeat.(scaleNames{i}) = zeros(numEdgePixels, 1);
    brightDarkSide.(scaleNames{i}) = zeros(numEdgePixels, 2);
end

%% Preprocess image
if size(im, 3) == 3         % Convert to grayscale
    im = rgb2gray(im);
end
im = im2double(im);         % Convert to double

%% Convolve image with all filters
filteredIm = cellfun(@(x) filter2(x, im, 'same'), orientedMasks, 'UniformOutput', 0);

% Loop over all scales
for scaleIdx = 1:numScales

    % Loop over all edge pixels
    for pixelIdx=1:numEdgePixels
        leftIntMean = filteredIm{edgePixelOrientIdx(pixelIdx), 1, scaleIdx}(edgePixelCrd(pixelIdx, 1), edgePixelCrd(pixelIdx, 2));
        rightIntMean = filteredIm{edgePixelOrientIdx(pixelIdx), 2, scaleIdx}(edgePixelCrd(pixelIdx, 1), edgePixelCrd(pixelIdx, 2));
        
        if leftIntMean > rightIntMean
            intensityRatioFeat.(scaleNames{scaleIdx})(pixelIdx) = leftIntMean / rightIntMean;
            brightDarkSide.(scaleNames{scaleIdx})(pixelIdx, :) = [1 2];
        else
            intensityRatioFeat.(scaleNames{scaleIdx})(pixelIdx) = rightIntMean / leftIntMean;
            brightDarkSide.(scaleNames{scaleIdx})(pixelIdx, :) = [2 1];
        end
    end
end