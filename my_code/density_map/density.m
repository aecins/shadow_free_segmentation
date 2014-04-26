function [slope, intercept, E] = density (im, r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the density map of an image.
%
% See:
%   A. Ecins, C. Fermüller, and Y.Aloimonos, "Shadow Free Segmentation in
%   Still Images using Local Density Measure", ICCP 2014
% 
% Input:
%   im,             input image 
%   r,              set of radii to be used for estimating mu
%
% Output:
%   slope,          density map
%                   i.e. the slope of the line in the log mu vs log r plot
%   intercept,      interecept of the line in the log mu vs log r plot
%   E,              error of fitting the stright lines to the datapoints in
%                   the log mu vs log r plot
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default radii
if nargin < 2
    r = 1:7;
end

%% Preprocessing

% Check grayscale
if (size(im,3) == 3),
    display('Warning: Only works on grayscale image, converted to grayscale image');
    im = rgb2gray(im);
end;

im = double(im);                % Convert to double
im = im*255 + 1;                % Scale to [1..256]. This is required to avoid log10(0) = -Inf
width = ceil(r(end)+2);
im = reflectIm(im, width);      % Reflect image boundaries    

%% Compute log(mu)
nScales = numel(r);                     % Number of scales
log_r = log10(r);                       % log of radius
log_mu = zeros([size(im), nScales]);    % log of the measurement function

for i=1:nScales
    
    % NOTE: in the paper the summation is done over a square of radius
    % r(i). Here instead we do a Gaussian weighted summation over a square
    % of the same radius.
    mask = fspecial('gaussian',r(i),r(i)/2) * r(i)^2;
    mu = conv2(im, mask, 'same');
    log_mu(:,:,i) = log10(mu);                                             % Take log
end

%% Find slope using least squares

% Convert to zero mean
x = log_r - sum(log_r(:))/nScales;
y = log_mu - repmat(sum(log_mu,3), [1 1 nScales])/nScales;

DD=0;
for m=1:nScales
    DD = DD + y(:,:,m).*x(m);
end

if (nScales ==1)
    slope = 10.^log_mu(:,:,1);
else
    slope = DD ./ (x*x');
end;

%% Find fitting error
E = zeros(size(im));
for m=1:nScales
    E = E + abs(y(:,:,m) - slope*x(m));
end

E = E / nScales;

%% Find y-intercept
tmp = zeros(size(log_mu));
for m=1:nScales
    tmp(:,:,m) = slope * log_r(m);
end

intercept = mean(log_mu - tmp, 3);

%% Cut off extra boundary
slope = slope(width+1:end-width, width+1:end-width);
intercept = intercept(width+1:end-width, width+1:end-width);
E = E(width+1:end-width, width+1:end-width);