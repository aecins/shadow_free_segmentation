%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train a shadow boundary classifier from the precomputed features. You can
% select the features to be used in the Classifier parameters section of
% this code.
% NOTE: the names training and prediction commands of libsvm were changed
%   svmtrain    -> libsvmtrain
%   svmpredict  -> libsvmpredict
% This is to avoid conflict with similar commands with the same name
% available in MATLAB Statistics toolbox. You can either do the same to
% your installation of LIBSVM or rename commands in this file to their
% originals.
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearAllPreserveBreakpoints;

%% Directories
featDirPath     = './dataset/features/';       % Directory with precomputed features
modelDirPath    = './models';                  % Directory containing trained models

%% Classifier parameters

% SVM parameters
svmParam = '-s 0 -t 2';        % RBF kernel

% Features used (set flag to 1 to use a feature, 0 to not use it)
featureUse.intRatio     = 1;
featureUse.RGBratio     = 1;
featureUse.chromAlign   = 0;
featureUse.labRatio     = 0;
featureUse.texture      = 1;

% Size of the halfdisc filters
scales = [8];
scaleStr = '';
for i=1:numel(scales)
    scaleNames{i} = genvarname(num2str(scales(i)));    
    scaleStr = [scaleStr '_' num2str(scales(i))];
end

% Martin's filterbank parameters
filterScale     = 0.5; 
filterElongate  = 3;

filterScaleStr = num2str(filterScale);
filterScaleStr = strrep(filterScaleStr, '.', 'x');
filterElongateStr = num2str(filterElongate);
filterElongateStr = strrep(filterElongateStr, '.', 'x');

% Texture mode used
texMode = [];
texMode{1} = sprintf('tex_texton128_malik_%s_%s_imagespec_density5', filterScaleStr, filterElongateStr);
% texMode{1} = sprintf('tex_texton128_malik_%s_%s_imagespec_intensity', filterScaleStr, filterElongateStr);
% texMode{1} = 'tex_lbp_pts8_rad1_u2_intensity';

% Model file name
modelFileName = sprintf('%d_%d_%d_%d_%d_%s%s',...
        featureUse.intRatio,...
        featureUse.RGBratio,...
        featureUse.chromAlign,...
        featureUse.labRatio,...
        featureUse.texture,...
        texMode{1},...
        scaleStr);
    
%% Prepare features
featDir = dir(fullfile(featDirPath, '*.mat'));

featuresAll = [];
labelsAll = [];

for imId = 1:numel(featDir)
    
    % Load features
    load(fullfile(featDirPath, featDir(imId).name));
    
    % Get current feature vector
    curFeatures = getFeatureVector(features, scaleNames, featureUse, texMode);

    % Concatenate current features and labels
    featuresAll = [featuresAll; curFeatures];
    labelsAll   = [labelsAll; features.label];
end

labelsAll = double(labelsAll);
labelsAll(labelsAll == 2) = -1;

%% Find SVM parameters using grid search (C and gamma)
bestcv  = 0;
tc = 2.^(1:5);      bestc = tc(1);
tg = 2.^(-4:3);     bestg = tg(1);
for c = tc,
    for g = tg,
        param = [svmParam ' -v 5 -c ', num2str(c), ' -g ', num2str(g), ' -q'];
        cv = libsvmtrain(labelsAll, featuresAll, param);
        if (cv >= bestcv),
            bestcv = cv; bestc = c; bestg = g;
        end
    end
end

%% Build model
param = [svmParam ' -b 1' ' -c ' num2str(bestc) ' -g ' num2str(bestg)];
curModel = libsvmtrain(labelsAll, featuresAll, param);

save(fullfile(modelDirPath, modelFileName), 'curModel');

%% Test on training data
[~, accuracy, probValues] = libsvmpredict(labelsAll, featuresAll, curModel, '-b 1');