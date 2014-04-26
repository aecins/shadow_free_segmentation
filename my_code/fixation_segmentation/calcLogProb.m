function logP = calcLogProb(X, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute log probability of the data given the model
%
% Input:
%       X,          data points
%       model,      GMM colour model
%
% Output:
%       logP,       log probability of each pixel according to the colour
%                   model
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    N = size(X,1);
    
    logP = zeros(N, model.K);
    
    % Calculate the probability of belonging to each cluster
    % NOTE: should replace mvnpdf with direct computation of multivariate
    % Gaussian distribution since mvnpdf requires Statistics Toolbox
    for k=1:model.K
        logP(:,k) = model.pClust(k) * mvnpdf(X, model.mu(k,:), model.sigma(:,:,k));
    end
    
    % Sum probabilities over clusters and take log
    logP = log(sum(logP, 2));
   
    % Replace -Inf with minimum value
    if any(abs(logP) == Inf)
        tmp = logP;
        tmp(tmp == -Inf) = Inf;
        logP(logP == -Inf) = min(tmp);
        
    end
end