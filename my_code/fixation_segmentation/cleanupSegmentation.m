function labelsCleaned = cleanupSegmentation(labels, fix_x, fix_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function cleans up segmentation results based on the following
% criteria:
% - There should be a single segmentation region
% - Segmentation region should contain the fixation point
% - All other regions that are not connected to it should be disregarded
% - Holes inside the segmented region should be filled in
%
% Input:
%   labels,             segmentation result
%   fix_x, fix_y,       fixation point coordinates
%
% Output:
%   labelsCleaned,      cleaned up segmentation
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelsCC = bwconncomp(labels, 8);

compSize = cellfun(@(x) numel(x), labelsCC.PixelIdxList);
idx = find(compSize == max(compSize));
tmp = zeros(size(labels));
tmp(labelsCC.PixelIdxList{idx}) = 1;
tmp = imfill(tmp, 'holes');

if tmp(fix_y, fix_x) == 1
    labels = tmp;
end

labelsCC = bwconncomp(labels, 8);
componentIdx = labelsCC.PixelIdxList;                               % Cell array with indices of each component
fixationIdx = sub2ind(size(labels), fix_y, fix_x);                      % Index of fixation point
for i=1:numel(componentIdx)
    if any(componentIdx{i} == fixationIdx)
        idx = i;
        break
    end
end

labelsCleaned = false(size(labels));
labelsCleaned(labelsCC.PixelIdxList{idx}) = true;

% Propagate labels near image boundaries
tmp = (labelsCleaned(:,2) == 1);
labelsCleaned(tmp, 1) = 1;
tmp = (labelsCleaned(:,end-1) == 1);
labelsCleaned(tmp, end) = 1;
tmp = (labelsCleaned(2,:) == 1);
labelsCleaned(1, tmp) = 1;
tmp = (labelsCleaned(end-1,:) == 1);
labelsCleaned(end, tmp) = 1;

% Further cleaning
if sum(labelsCleaned(:)) > 5
    labelsCleaned = bwmorph(labelsCleaned, 'majority', Inf);
end

% Fill in holes
% NOTE: this is a hack to fill in holes in the segmentaion. We just find
% all regions that are smaller than some fraction of the segmented region
% and set them as foreground.
labelsCC = bwconncomp(~labelsCleaned);
areas = regionprops(labelsCC, 'Area');
idx = [areas.Area] < floor(sum(labelsCleaned(:)) / 20);
labelsCleaned(cat(1,labelsCC.PixelIdxList{idx})) = true;