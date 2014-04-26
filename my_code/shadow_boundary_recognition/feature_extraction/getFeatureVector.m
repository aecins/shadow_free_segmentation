function featureVector = getFeatureVector(features, scaleNames, featureUse, texMode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a feature vector given the precomputed feature strcuture and features
% used.
%
% Input:
%   features,           precomputed feature datastructure
%   scaleNames,         
%   featureUse,         vector indicating which features are used
%   texMode,            texture mode used
%
% Output:
%   featureVector,      feature vector
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Scale features
for s = scaleNames'
    features.intRatio.(s{1})      = scaleRatioFeature(features.intRatio.(s{1}));
    features.RGBratio.(s{1})      = scaleRatioFeature(features.RGBratio.(s{1}));
    features.chromAlign.(s{1})    = scaleRatioFeature(features.chromAlign.(s{1}));
    features.labRatio.(s{1})      = scaleRatioFeature(features.labRatio.(s{1}));
end

%% Construct feature vector

% Combine all non-texture features
featuresNonTexture = [];
for curFeatureType = fieldnames(featureUse)'
    if ~strcmp(curFeatureType{1}, 'texture') && featureUse.(curFeatureType{1})
        for s=scaleNames'
            featuresNonTexture = [featuresNonTexture features.(curFeatureType{1}).(s{1})];
        end
    end
end

% Add texture features
featuresTexture = [];
if featureUse.texture        
    for s=scaleNames'            
        featuresTexture = [featuresTexture features.texture.(texMode{1}).(s{1})];
    end
end

% Concatenate features
featureVector = [featuresNonTexture featuresTexture];