function visualizeMasks(im, orientedMasks, edgePixelIdx, edgePixelCrd, edgePixelOrientIdx, gX, gY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize boundaries boundary perpendiculars and corresponding halfdiscs.
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curMaskIdx = 1;
radius = (size(orientedMasks{1,1,curMaskIdx}, 1) - 1) / 2;
[yRes, xRes, ~] = size(im);

imVis = im;
idxList = [];

% Go through all pixels
for i=1:25:numel(edgePixelIdx)
    
    idxList = [idxList, i];
    
    % Get coordinates of window defined by mask
    leftX   = edgePixelCrd(i,2)-radius;
    rightX  = edgePixelCrd(i,2)+radius;
    topY    = edgePixelCrd(i,1)-radius;
    bottomY = edgePixelCrd(i,1)+radius;
    
    % Check that window is within the bounds of the image
    if leftX < 1 || rightX > xRes || topY < 1 || bottomY > yRes
        continue
    end
    
    % Get left and right masks
    maskL = zeros(yRes, xRes);
    maskL(topY:bottomY, leftX:rightX) = orientedMasks{edgePixelOrientIdx(i), 1, curMaskIdx} > 0;
    maskR = zeros(yRes, xRes);
    maskR(topY:bottomY, leftX:rightX) = orientedMasks{edgePixelOrientIdx(i), 2, curMaskIdx} > 0;
    
    tmpMask = zeros(yRes, xRes);
    tmpMask(edgePixelIdx) = 1;
    imVis = imRegionHighlight(imVis, logical(tmpMask), 'w');
    imVis = imRegionHighlight(imVis, maskL*0.5, 'r'); 
    imVis = imRegionHighlight(imVis, maskR*0.5, 'b'); 
    

end

imshow(imVis);
hold on;
quiver(edgePixelCrd(idxList,2), edgePixelCrd(idxList,1), gX(idxList), gY(idxList), 1, 'g');