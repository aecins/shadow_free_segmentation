function [boundaries, internalBoundaries] = segmentBoundary(segMask, im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function analyzes the segmentation area and extracts:
%   - segmentation boundary segments
%   - internal boundaries of the segmentation area
% Perpendiculars to segmentation boundary segments are computed so that
% they can be used later for extracting shadow detection features.
%
% Input:
%   segMask,                binary mask of the segmentation area
%   im,                     input image
%
% Output:
%   boundaries,             segmentation boundary segments
%   internalBoundarries,    boundaries from oversegmenting the segmentation
%                           area
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get internal boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Oversegment the segmentation area
tmp = imRegionHighlight(im, ~segMask, [0 0 0]);
minArea = floor(sum(segMask(:))/20);
[~, meanshiftSeg] = edison_wrapper(tmp, @RGB2Luv, ...
'SpatialBandWidth', 9, 'RangeBandWidth', 15, ...
 'MinimumRegionArea', minArea);

% Find internal boundaries of segments belonging to object
internalBoundaries = false(size(segMask));
for j = 1:max(max(meanshiftSeg(:)))
    maskTmp = meanshiftSeg == j;
    maskTmp = bwmorph(maskTmp, 'remove');
    internalBoundaries = internalBoundaries | maskTmp;
end
internalBoundaries = internalBoundaries - bwmorph(segMask, 'remove');
internalBoundaries = logical(internalBoundaries);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get segmentation boundary segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
numAdjacentPixels = 5;                          % Number of neighbouring points used to estimate contour perpendiculars
minBndLen = 20;                                 % Minimum boundary segment length
curvatureSpacing = 10;                          % Distance between neighbours used to estimate curvature

%% Find perpendiculars to boundary

segBoundary = bwmorph(segMask, 'remove');                                   % Segmentation boundary
sortedIdx = traceEdge(segBoundary);                                         % Order pixels
orient = zeros(numel(sortedIdx), 2);
[contourY, contourX] = ind2sub(size(segMask), sortedIdx);
contourX = [contourX(end-numAdjacentPixels:end); contourX; contourX(1:numAdjacentPixels)];
contourY = [contourY(end-numAdjacentPixels:end); contourY; contourY(1:numAdjacentPixels)];

% Loop through segmentation boundary pixels
for i=1:numel(sortedIdx)

    % Find edge pixel orientation
    neighbours = i:i+2*numAdjacentPixels;
    X = [contourX(neighbours) contourY(neighbours)];
    coeff = princomp(X);
    contourPerp = coeff(:,2)';

    % Make sure it is consistent with previous measurement
    if i > 1
        dist            = sum((contourPerp - orient(i-1,:)).^2);
        reversedDist    = sum((-contourPerp - orient(i-1,:)).^2);
        if reversedDist < dist
            contourPerp = -contourPerp;
        end
    end

    % Save to variable
    orient(i,:) = contourPerp;    
end

%% Find curvature
% NOTE: the 2D curvature of a point in a curve is estimated by finding the
% radius of the circle circumscribed in the triangle defined by the point
% itself and it's left and right neighbours that are at fixed distance from
% it.

curvature = zeros(numel(sortedIdx), 1);
[contourY, contourX] = ind2sub(size(segMask), sortedIdx);
contourX = [contourX(end-curvatureSpacing:end); contourX; contourX(1:curvatureSpacing)];
contourY = [contourY(end-curvatureSpacing:end); contourY; contourY(1:curvatureSpacing)];

% Loop through the point of the segmentation boundary
for i=1:numel(sortedIdx)
    
    % Get 3 triangle points
    leftPt      = [contourX(i); contourY(i)];
    centerPt    = [contourX(i+curvatureSpacing); contourY(i+curvatureSpacing)];    
    rightPt     = [contourX(i+2*curvatureSpacing); contourY(i+2*curvatureSpacing)];
        
    % Compute triangle side lengths
    a = sqrt(sum((leftPt    - centerPt).^2));
    b = sqrt(sum((centerPt  - rightPt).^2));
    c = sqrt(sum((rightPt   - leftPt).^2));
    
    % Compute traingle area using Heron's formula
    s = (a+b+c) / 2;
    S = sqrt( s * max(s-a,0) * max(s-b,0) * max(s-c,0));
    
    % Compute curvature
    curvature(i) = 4*S / (a*b*c);

end

%% Find curvature maxima
[~, peakIdx] = findpeaks(curvature, 'MINPEAKDISTANCE', minBndLen, 'MINPEAKHEIGHT', 0.03);

%% Find segments
boundaries = cell(numel(peakIdx), 1);
for i=1:numel(peakIdx)-1
    boundaries{i}.pixelIdx  = sortedIdx(peakIdx(i):peakIdx(i+1)-1);
    boundaries{i}.orient    = orient(peakIdx(i):peakIdx(i+1)-1,:);
    boundaries{i}.label     = NaN;
end
boundaries{i+1}.pixelIdx    = [sortedIdx(peakIdx(i+1):end); sortedIdx(1:peakIdx(1)-1)];
boundaries{i+1}.label       = NaN;
boundaries{i+1}.orient      = [orient(peakIdx(i+1):end,:); orient(1:peakIdx(1)-1,:)];