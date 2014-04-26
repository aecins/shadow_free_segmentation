function [FGweights, BGweights] = computeUnaryWeights(im, labels, mode, nSample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the unary terms used in min-cut minimization. Two
% colour models are constructed using GMM: one describing the background
% the other the foreground. Then for each pixel the probability of
% observing that pixel given foreground/background model is computed. The
% weights are set to be the negative of natural log of this probability.
%
% Input:
%   im,         image
%   labels,     labels describing which parts of the image belong to the
%               foreground and which to the background
%   mode,       colour feaure mode
%   nSample,    number of sample used by k-means to get initial colour
%               model
%   tolerance,  determinant tolerance for covariance matrix regularization
%               in GMM
% Output:
%   FGweights,  foreground weights
%   BGweights,  background weights
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check input arguments
if nargin < 4
    nSample = 0;
end

if nargin < 3
    mode = 'hs';
end

nPix = size(im,1) * size(im,2);

switch mode

    % RGB
    case 'rgb'
        X = im(:);
        X = [X(1:nPix), X(nPix+1:2*nPix), X(2*nPix+1:end)];
    
    % Hue and saturation
    case 'hs'
        imHSV = rgb2hsv(im);        
        X = imHSV(:);
        X = [X(1:nPix), X(nPix+1:2*nPix)];
        
    % L*a*b colorspace
    case 'lab'
        C = makecform('srgb2lab');
        X = applycform(im2uint8(im),C);
        X = im2double(X(:));
        X = [X(1:nPix), X(nPix+1:2*nPix), X(2*nPix+1:end)];

    % Chromaticity channels of L*a*b
    case 'ab'
        C = makecform('srgb2lab');
        X = applycform(im2uint8(im),C);
        X = im2double(X(:));
        X = [X(nPix+1:2*nPix), X(2*nPix+1:end)];        
        
    otherwise
        
        error('Unsupported colour mode %s\n', mode);
        return;
end

%% Get pixels for foreground and background
FG = X(labels == 1, :);
BG = X(labels == 0, :);

%% Calculate colour models for foreground and background
if nSample == 0
    modelFG = colourGMM(FG, size(FG,1));
    modelBG = colourGMM(BG, size(BG,1));
else
    modelFG = colourGMM(FG, nSample);
    modelBG = colourGMM(BG, nSample);
end

%% Calculate unary weights for all pixels (probability of belonging to foreground or background)
if ~isempty(modelFG)
    BGweights =  -calcLogProb(X, modelFG);
%     BGweights = replaceInfs(BGweights);
else
    BGweights = NaN;
end

if ~isempty(modelBG)
    FGweights =  -calcLogProb(X, modelBG);
%     FGweights = replaceInfs(FGweights);    
else
    FGweights = NaN;
end