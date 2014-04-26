function textureFeat = getTextureFeatures(im, scales, orientedMasks, edgePixelCrd, edgePixelOrientIdx, texMode, fbParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract texture features
% 
% Input:
%   im,                 input image (colour or grayscale)
%   sacles,             halfdisc scales
%   orientedMasks,      halfdisk masks
%   edgePixelCrd,       coordinates of boundary pixels
%   edgePixelOrientIdx, orientation indeces of boundary pixels
%   texMode,            texture description mode used
%   fbParams,           parameters of Martin's filterbank
%
% Output:
%   textureFeat,        texture feature datastructure
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numScales = numel(scales);
numEdgePixels = size(edgePixelCrd, 1);
scaleNames = cell(numScales, 1);

for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));
    textureFeat.(scaleNames{i}) = zeros(numEdgePixels, 1);
end

%% Preprocess image
if size(im, 3) == 3         % Convert to grayscale
    im = rgb2gray(im);
end
im = im2double(im);         % Convert to double

%% Get the required texture descriptor map

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEXTONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(texMode, 'texton')

    % Convert image if needed
    if strfind(texMode, 'intensity')               % Intensity
        % do nothing
    elseif strfind(texMode, 'density')             % Density
        r = findParameter(texMode, 'density');
        if isempty(r)
            r = [3 9 17];
            imDensity = zeros([size(im(:,:,1)) numel(r)]);
            for i=1:numel(r)
                imDensity(:,:,i) = density(im, 1:r(i));
            end
            [imDensity, maxIdx] = max(imDensity, [], 3);
        else
            im = density(im, 1:r);
        end
    end

    % Create filters
    [scale, elongate] = getFilterParameters(texMode);
    fb = fbCreate(fbParams.no, scale, fbParams.ns, fbParams.sc, elongate);

    % Filter image
    fim = fbRun(fb, im);

    %% Universal
    if ~isempty(strfind(texMode, 'universal'))

        % Assign universal textons
        textonFile = fullfile(folderPath, [texMode, '.mat']);
        textonData = load(textonFile);                                              % Load universal texton dictionary
        textons = textonData.tex;
        nTexels = size(textons, 2);
        texelMap = assignTextons(fim, textons);                                    % Assign textons
    end

    %% Imagespec
    if ~isempty(strfind(texMode, 'imagespec'))

        % Cluster textons
        nTexels = findParameter(texMode, 'texton');                 % Number of textons
        [texelMap, ~] = computeTextons(fim,nTexels);                % Compute textons
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LBP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strfind(texMode, 'lbp')

    % Convert image if needed
    if strfind(texMode, 'intensity')               % Intensity
        % do nothing
    elseif strfind(texMode, 'density')             % Density
        r = findParameter(texMode, 'density');
        im = density(im, 1:r);
    end    

    % Directory containing LBP mappings
    mappingDir = fullfile(homeroot, '/Code/lbp/mappings');

    % Get LBP parameters
    R = findParameter(texMode, 'rad');
    P = findParameter(texMode, 'pts');

    if ~isempty(strfind(texMode, 'riu2'))
        mappingMode = 'riu2';
    elseif ~isempty(strfind(texMode, 'ri'))
        mappingMode = 'ri';
    elseif ~isempty(strfind(texMode, 'u2'))
        mappingMode = 'u2';
    end

    % Loop over all radii
    texelMap = zeros([size(im), numel(R)]);
    for i=1:numel(R)

        % Load lbp mapping data
        load(sprintf('%s/%s_p%d.mat', mappingDir, mappingMode, P(i)));             

        % Get LBP for all image pixels
        im_ = reflectIm(im, R(i));
        texelMap(:,:,i) = lbp(im_,R(i),P(i),mapping, '');
        nTexels(i) = mapping.num;
    end    
end
    
%% Find textons histograms in halfdiscs and find distances between them

% Loop over all edge pixels
for scaleIdx = 1:numScales
    radius = (size(orientedMasks{1,1,scaleIdx}, 1) - 1) / 2;
%     fprintf('\tRadius %d...\n', radius);
    
    % Loop over pixels
    for pixelIdx=1:numEdgePixels
        
%         fprintf('%d/%d\n', pixelIdx, numEdgePixels);
        
        % Get coordinates of window defined by mask
        leftX   = edgePixelCrd(pixelIdx,2)-radius;
        rightX  = edgePixelCrd(pixelIdx,2)+radius;
        topY    = edgePixelCrd(pixelIdx,1)-radius;
        bottomY = edgePixelCrd(pixelIdx,1)+radius;

        % Get left and right masks
        maskL = zeros(size(im));
        maskL(topY:bottomY, leftX:rightX) = orientedMasks{edgePixelOrientIdx(pixelIdx), 1, scaleIdx} > 0;
        maskR = zeros(size(im));
        maskR(topY:bottomY, leftX:rightX) = orientedMasks{edgePixelOrientIdx(pixelIdx), 2, scaleIdx} > 0;

        % Collect histograms
        histLCat = [];
        histRCat = [];
        for j=1:size(texelMap,3)
            curTexelMap = texelMap(:,:,j);
            histL = histc(curTexelMap(logical(maskL)), 1:nTexels(j))';
            histL = histL / sum(histL);
            histR = histc(curTexelMap(logical(maskR)), 1:nTexels(j))';
            histR = histR / sum(histR);
            
            histLCat = [histLCat histL];
            histRCat = [histRCat histR];
        end
        
        histLCat = histLCat / numel(nTexels);
        histRCat = histRCat / numel(nTexels);
        
        % Calculate chi2 distance
        textureFeat.(scaleNames{scaleIdx})(pixelIdx) = chi2dist(histLCat', histRCat');
        
    end
end