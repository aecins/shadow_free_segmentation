function rgbRatioFeat = getRGBratioFeatures(im, scales, orientedMasks, edgePixelCrd, edgePixelOrientIdx, brightDarkSide)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract RGB colour channel ratio features.
% 
% Input:
%   im,                 input image (colour or grayscale)
%   sacles,             halfdisc scales
%   orientedMasks,      halfdisk masks
%   edgePixelCrd,       coordinates of boundary pixels
%   edgePixelOrientIdx, orientation indeces of boundary pixels
%   brightDarkSide,     datastructure indicating which of the two halfdisks
%                       is dark and which is bright
%
% Output:
%   rgbRatioFeat,       RGB ratio feature datastructure
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numScales = numel(scales);
numEdgePixels = size(edgePixelCrd, 1);
scaleNames = cell(numScales, 1);

for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));
    rgbRatioFeat.(scaleNames{i}) = zeros(numEdgePixels, 3);
end

%% Preprocess image
im = im2double(im);         % Convert to double

%% Convolve image with all filters
filteredIm{1} = cellfun(@(x) filter2(x, im(:,:,1), 'same'), orientedMasks, 'UniformOutput', 0);
filteredIm{2} = cellfun(@(x) filter2(x, im(:,:,2), 'same'), orientedMasks, 'UniformOutput', 0);
filteredIm{3} = cellfun(@(x) filter2(x, im(:,:,3), 'same'), orientedMasks, 'UniformOutput', 0);

% Loop over scales
for scaleIdx = 1:numScales
        
    % Loop over colour channels
    for clrChannel = 1:3

        % Loop over all edge pixels
        for pixelIdx=1:numEdgePixels
            
            y = edgePixelCrd(pixelIdx, 1);
            x = edgePixelCrd(pixelIdx, 2);
            orient = edgePixelOrientIdx(pixelIdx);
            brightSide  = brightDarkSide.(scaleNames{scaleIdx})(pixelIdx, 1);
            darkSide    = brightDarkSide.(scaleNames{scaleIdx})(pixelIdx, 2);
            
            brightClrMean = filteredIm{clrChannel}{orient, brightSide, scaleIdx}(y, x);
            darkClrMean = filteredIm{clrChannel}{orient, darkSide, scaleIdx}(y, x);
            
            rgbRatioFeat.(scaleNames{scaleIdx})(pixelIdx, clrChannel) = brightClrMean / darkClrMean;
        end
    end
end