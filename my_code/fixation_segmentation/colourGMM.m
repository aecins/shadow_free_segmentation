function model = colourGMM(X, sampleSize, nClusters, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the GMM colourmodel of the input pixel array.
%
% Input:
%   X,              input pixel array
%   sampleSize,     number of samples used to estimate the model
%   nClusters,      number of clusters used
%   opts,           display flag
%
% Output:
%   model,          GMM model
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check number of input argumetns
if nargin < 4
    displayFlag = false;
else
    switch opts
        case 'display'
            displayFlag = true;
            
        otherwise
            warning ('Unrecognized option\n');
    end
end

if nargin < 3
    nClusters = 5;
end

if (nargin < 2) || (nargin > 1 && isempty(sampleSize))
    sampleSize = size(X,1);
end

%% Check if there are enough pixels to build a model
if size(X,1) < nClusters
    warning('Too few points to compute colour model');
    model = [];
    return;
end

%% Downsample data
if (sampleSize < size(X,1))
    
    % Permute pixels randomly
    select = randperm(size(X,1));
    select = select(1:sampleSize);
    X = X(select,:);
end

%% Parameters

% Model parameters
model.nDim = size(X,2);                                                    % Number of dimensions
model.K = nClusters;                                                       % Numer of clusters
model.mu = zeros(model.nDim,model.K);                                      % Cluster means
model.sigma = zeros(model.nDim, model.nDim, model.K);                      % Covariance 
model.pClust = zeros(1,model.K);                                           % Cluster probabilities

% Other parameters
N = size(X,1);                                                             % Number of pixels
gamma = zeros(N, model.K);                                                 % Resposabilties

%% Estimate GMM parametrs using k-means

% Run k-means
[IDX, model.mu] = kmeans(X, model.K, 'replicates', 5);

% Calculate covariance
for k=1:model.K
    model.sigma(:,:,k) = cov(X(IDX==k,:));
end

% Cluster probabilities
for k=1:model.K
    model.pClust(k) = sum(IDX==k) / numel(IDX);
end

% Display
if displayFlag
    if ~sum(gamma(:)) == 0
        [~, IDX] = max(gamma, [], 2);
    end
    subplot(1,2,2);
    displayGMM(X, IDX, model);
    title('GMM after convergence');
end

% Display
if displayFlag
    figure;
    subplot(1,2,1);
    displayGMM(X, IDX, model);
    title('GMM at initialization');
end