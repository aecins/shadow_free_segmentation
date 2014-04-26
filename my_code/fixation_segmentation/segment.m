function [labels, unaryColourFG, unaryColourBG] = segment(im, edgeGrad, fix_x, fix_y, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of the fixation based figure ground
% segmentation algorithm presented in the paper:
%   A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in
%   Still Images using Local Density Measure", ICCP 2014
%
% The original fixation based segmentation algorithm is by Ajay Misra:
%   Mishra, Ajay, Yiannis Aloimonos, and Cornelia Fermuller. "Active
%   segmentation for robotics." Intelligent Robots and Systems, 2009. IROS 
%   2009. IEEE/RSJ International Conference on. IEEE, 2009.
%
% Input:
%   im,             input image (must be double)
%   edgeGrad,       image probabilistic edge map
%   fix_x,          fixation point x coordinate
%   fix_y,          fixation point y coordinate
%   opt,            options (verbose)
%
% Output:
%   labels,         binary figure ground labels
%   unaryColourFG,  foreground colourmodel probability map
%   unaryColourBG,  background colourmodel probability map
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



verbose = false;

if nargin > 4
    if strcmp(opt, 'verbose')
        verbose = true;
    end
end

%% Energy fucntion parameters

% Binary weights
nju = 5;                    % Binary weight exponent constant
k = 20;                     % Binary weight for two pixels with zero edge probability
lambda = 1000;              % Importance of binary weights

% Unary weights
D  = 1e100;                 % Fixation point unary weight

% This is a small hack that allows us to add additional point that we think
% are likely to belong to the object (i.e. high probability non-shadow
% edges).
D_ = 10;

unaryColourWeight = 1;      % Importance of colourmodel unary weights

foreground = 2;             % Object
background = 1;             % Background

%% Image parameters
[yRes, xRes] = size(edgeGrad);
nPix = xRes*yRes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST ITERATION. EDGE INFORMATION ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose 
    fprintf('First iteration. Using edge information only\n'); 
end

%% Compute binary weights
if verbose
    fprintf('Computing binary weight matrix\n');
end
E = edges4connected(yRes,xRes);                                             % Indices of adjacent pixels (Potts model)
I_pq = edgeGrad(E(:,1))+edgeGrad(E(:,2)) / 2;                               % Average edge probability at adjacent edges
V = zeros(size(E,1),1);                                                     % Weights
V(I_pq ~= 0) = exp( -nju*( I_pq(I_pq ~= 0) ));                              % Edge where at least one of the pixels belongs to the edge map
V(I_pq == 0) = k;                                                           % Edges where none of the pixels belong to the edge map

% Calculate the distance of each edge from the fixation point
[y_p, x_p] = ind2sub([yRes, xRes], E(:,1));                               % Calculate midpoint of each edge
[y_q, x_q] = ind2sub([yRes, xRes], E(:,2));
x_mid = (x_p + x_q) / 2 - fix_x(1);
y_mid = (y_p + y_q) / 2 - fix_y(1);
r = sqrt(x_mid.^2 + y_mid.^2);

% Weights are the inverse of the distance from the fixation point
W = 1./r;
W = W/max(W);                                                               % Normalize the weights to have maximum of 1
V = V.*W;

A = sparse(E(:,1),E(:,2),V,nPix,nPix,4*nPix);

%% Construct unary weights for image boundary and fixation point 
T_ = zeros(numel(edgeGrad),2);
T_(1:yRes,background) = D;                                      % Left column
T_(end-yRes+1:end,background) = D;                              % Right column
T_((0:xRes-1)*yRes+1,background) =D;                                        % Top row
T_((1:xRes)*yRes,background) =D;                                            % Bottom row
T_(sub2ind([yRes, xRes], fix_y(1), fix_x(1)), foreground) = D;            % Fixation point
T_(sub2ind([yRes, xRes], fix_y(2:end), fix_x(2:end)), foreground) = D_;   % Fixation point
T = sparse(T_);

%% Perform min-cut
if verbose
    fprintf('Performing min-cut...\n');
end

[~, labels] = maxflow(A,T);
labels = reshape(labels, [yRes, xRes]);
labels = cleanupSegmentation(labels, fix_x(1), fix_y(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND ITERATION. USING COLOUR MODEL DERIVED FROM SEGMENTATION IN FIRST ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('Second iteration. Using colour model\n');
end

% If first segmentation is too small we can not compute background model
% hence we can not do second segmentation iteration
if sum(labels(:) == 1) < 5
    fprintf('Area of segmentation obtained in first iteration is too small. Returning results from first iteration\n');
    labels = reshape(labels, [yRes, xRes]);
    segBoundary = bwmorph(labels==0, 'remove');
    segIm = imRegionHighlight(im, segBoundary == 1, 'g');
    unaryColourFG = [];
    unaryColourBG = [];
    return;
end

%% Set the importance of binary terms
A = A * lambda;

% Loop until convergence (similar to Grab Cut approach). We iterate until
% the overlap between the previous and current segmentations is above some
% threshold.

itNum = 1;
overlap = 0;

while ~(overlap>0.9)
    
    if verbose
        fprintf('Iteration %d:', itNum);
    end
    itNum = itNum+1;
    labels_prev = labels;
   
    % Compute unary terms from color model
    if verbose
        fprintf(' Computing unary terms...');
    end    
    
    [unaryColourFG, unaryColourBG] = computeUnaryWeights(im, labels(:), 'lab', 100000);  

    % Set unary colourmodel weights
    T_ = zeros(numel(edgeGrad),2);
    T_(:,foreground) = unaryColourWeight * unaryColourFG;
    T_(:,background) = unaryColourWeight * unaryColourBG;

    % Set background and foreground weights for image boundary and fixation point
    T_(1:yRes,background) = D;                                              % Left column
    T_(end-yRes+1:end,background) = D;                                      % Right column
    T_((0:xRes-1)*yRes+1,background) =D;                                    % Top row
    T_((1:xRes)*yRes,background) =D;                                        % Bottom row
    T_(sub2ind([yRes, xRes], fix_y(1), fix_x(1)), foreground) = D;          % Fixation point
    T_(sub2ind([yRes, xRes], fix_y(2:end), fix_x(2:end)), foreground) = D_; % Fixation point
    T = sparse(T_);

    %% Perform min-cut
    if verbose
        fprintf('Performing min-cut...\n');
    end
    [~, labels] = maxflow(A,T);

    labels = reshape(labels, [yRes, xRes]);
    labels = cleanupSegmentation(labels, fix_x(1), fix_y(1));
    
    % Find overlap between segmentations
    intersection    = labels & labels_prev;
    union           = labels | labels_prev; 
    overlap = sum(intersection(:)) / sum(union(:));

end

% Reshape unary weights
unaryColourFG = reshape(unaryColourFG, [yRes, xRes]);
unaryColourBG = reshape(unaryColourBG, [yRes, xRes]);

if verbose
    fprintf('Done!\n');
end