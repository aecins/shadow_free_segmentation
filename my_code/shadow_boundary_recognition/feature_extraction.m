%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract features used for shadow boundary classification. Multiple
% texture and colour features are extracted and saved to files in
% featDirPath directory. During classifier testing a subset of these
% features can be used for training.
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearAllPreserveBreakpoints;

%% Directories
imDirPath           = './dataset/images/';                 % Image directory
bndDirPath          = './dataset/boundaries/';             % Boundary directory
featDirPath         = './dataset/features/';               % Feature directory

%% Create bnd directory if it doesn't exist
if ~exist(featDirPath, 'dir')
    mkdir(featDirPath);
end

%% Feature extraction parameters

% Size of the halfdisc filters
scales = [8];

% Martin's filterbank parameters
filterScale     = 0.5; 
filterElongate  = 3;

%% Strings for feature parameters used for naming
scaleStr = '';
for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));    
    scaleStr = [scaleStr '_' num2str(scales(i))];
end

filterScaleStr = num2str(filterScale);
filterScaleStr = strrep(filterScaleStr, '.', 'x');
filterElongateStr = num2str(filterElongate);
filterElongateStr = strrep(filterElongateStr, '.', 'x');

%% Texture descriptor used
texMode = [];
texMode{1} = sprintf('tex_texton128_malik_%s_%s_imagespec_density5', filterScaleStr, filterElongateStr);
texMode{2} = sprintf('tex_texton128_malik_%s_%s_imagespec_intensity', filterScaleStr, filterElongateStr);
texMode{3} = sprintf('tex_lbp_pts8_rad1_u2_intensity');

%% Loop through all images and extract features
% You can run the for loop below in parallel. To do it just replace 'for'
% with parfor' (assuming that you have a multicore processor and Matlab
% Parallel Computing Toolbox)
bndDir = dir(fullfile(bndDirPath, '*.mat'));
featuresConcat = cell(numel(bndDir), 1);
parfor imId = 1:numel(bndDir)
    
    imName = bndDir(imId).name(1:end-4);
    
    %% Read image
    im = imread(fullfile(imDirPath, [imName '.png']));
    im = im2double(im);
    
    %% Load boundaries
    tmp = load(fullfile(bndDirPath, [imName '.mat']));
    boundaries = tmp.boundaries;
    
    tic
    featuresConcat{imId} = getFeatures(im, boundaries, scales, texMode);
    fprintf('Image %s - ', imName);
    toc
end

%% Save features
for imId=1:numel(bndDir)
    
    imName = bndDir(imId).name(1:end-4);
    featureFileName = fullfile(featDirPath, [imName '.mat']);

    features = featuresConcat{imId};
    save(featureFileName, 'features');
end