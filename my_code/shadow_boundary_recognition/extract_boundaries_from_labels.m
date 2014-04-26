function extract_boundaries_from_labels(startId, endId)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a labelled edge map extract a set of corresponding boundaries. Here
% a boundary is a data structure that corresponds to a single continuous
% contour of the same label in the edge map and holds information about:
%   pixelIdx:   contour pixel coordinates
%   orient:     contour perpendicular at each pixel
%   label:      label of the contour (1 = shadow, 2 = non-shadow)
% Boundaries are saved to files in bndDirPath directory.
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearAllPreserveBreakpoints;

%% Parameters
numAdjacentPixels = 10;                                 % Number of neighbours used to estimate contour orientation
minBndLength = numAdjacentPixels*2 + 15;                % Minimum accepted length of the boundary
maxBndLength = 100;                                     % Maximum accepted length of the boundary

%% Directories
labelDirPath = './dataset/labels/';                    % Directory containing labelled edge maps
bndDirPath = './dataset/boundaries/';                  % Output directory
labelDir = dir(fullfile(labelDirPath, '*.png'));

% If no range is provided process all edge maps in the input directory
if nargin < 2
    startId = 1;
    endId = numel(labelDir);
end
endId = min(endId, numel(labelDir));

% Count number of shadow and non shadow boundaries
numShadowBnd = 0;
numNonShadowBnd = 0;

%% Create bnd directory if it doesn't exist
if ~exist(bndDirPath, 'dir')
    mkdir(bndDirPath);
end

%% Loop through all images
for imId = 1:endId-startId+1
    
    imName = labelDir(imId+startId-1).name(1:end-8);
    fprintf('%s - %d/%d\n', imName, imId, endId-startId+1);
    labelIm = imread(fullfile(labelDirPath, [imName '_lbl.png']));
    
    % Get all edge fragments
    bndCC = bwconncomp(labelIm == 1 | labelIm == 2);
    
    % Remove short pixels
    shortBnd = cellfun(@(x) numel(x), bndCC.PixelIdxList) < minBndLength;
    bndCC.PixelIdxList(shortBnd) = [];
    bndCC.NumObjects = numel(bndCC.PixelIdxList);
    
    %% Loop through the edge fragments, find orientation
    for i=1:bndCC.NumObjects
        edgeMask = zeros(size(labelIm));
        edgeMask(bndCC.PixelIdxList{i}) = 1;
        sortedIdx = traceEdge(edgeMask);
        [contourY, contourX] = ind2sub(size(edgeMask), sortedIdx);
         
        % Find edge pixel orientation based on neighbours
        curBnd.pixelIdx = [];
        curBnd.orient   = [];
        for j=numAdjacentPixels+1 : numel(sortedIdx)-numAdjacentPixels
            
            % Find edge pixel orientation
            neighbours = j-numAdjacentPixels:j+numAdjacentPixels;
            X = [contourX(neighbours) contourY(neighbours)];
            coeff = princomp(X);
            contourPerp = coeff(:,2)';
           
            % Make sure it is consistent with previous measurement
            if ~isempty(curBnd.orient)
               
                dist            = sum((contourPerp - curBnd.orient(end,:)).^2);
                reversedDist    = sum((-contourPerp - curBnd.orient(end,:)).^2);
                if reversedDist < dist
                    contourPerp = -contourPerp;
                end
            end
                        
            % Save to variable
            curBnd.pixelIdx = [curBnd.pixelIdx; sortedIdx(j)];
            curBnd.orient   = [curBnd.orient; contourPerp];
        end        
        
        bndCC.PixelIdxList{i} = curBnd.pixelIdx;
        bndCC.OrientList{i} = curBnd.orient;
        
    end
    
    % Break up boundaries
    boundaries = {};
    orient = {};
    
    numBnd = 1;
    
    %% Loop through the edge fragments, break up those that are too long
    for i=1:bndCC.NumObjects
        if numel(bndCC.PixelIdxList{i}) <= maxBndLength
            boundaries{numBnd}.pixelIdx = bndCC.PixelIdxList{i};
            boundaries{numBnd}.orient = bndCC.OrientList{i};
            numBnd = numBnd + 1;
        else
            numNewBnd = ceil(numel(bndCC.PixelIdxList{i}) / maxBndLength);
            newBndLength = floor(numel(bndCC.PixelIdxList{i}) / numNewBnd);
            
            for j=1:numNewBnd
                boundaries{numBnd}.pixelIdx = bndCC.PixelIdxList{i}((j-1)*newBndLength+1 : j*newBndLength);
                boundaries{numBnd}.orient   = bndCC.OrientList{i}((j-1)*newBndLength+1 : j*newBndLength, :);
                numBnd = numBnd + 1;
            end
        end
    end
    
    %% Assign labels to boundaries
    for i=1:numel(boundaries)
        boundaries{i}.label = labelIm(boundaries{i}.pixelIdx(1));
    end

    %% Count number of shadow and non-shadow boundaries
    numShadowBnd = numShadowBnd + sum(cellfun(@(x) x.label == 1, boundaries));
    numNonShadowBnd = numNonShadowBnd + sum(cellfun(@(x) x.label == 2, boundaries));
        
    %% Save edge information to file
    save(fullfile(bndDirPath, [imName '.mat']), 'boundaries');
end

fprintf('%d shadow boundaries, %d non shadow boundaries\n', numShadowBnd, numNonShadowBnd);