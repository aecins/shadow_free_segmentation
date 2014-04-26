function [orientedMasks, orientations] = getOrientedMasks(nOrient, radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute binary masks for the halfdisk filters.
%
% Input:
%   nOrient,            number of halfdisk filter orientations
%   radius,             radius of the halfdisk filters
%
% Output:
%   orientedMasks,      halfdisk filter masks
%   orientations,       halfdisk filter orientation angles in radians
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get filter orientations
orientations = linspace(0, pi, nOrient+1);
orientations = orientations(1:end-1);

% Variable holding the masks
orientedMasks = cell(nOrient, 2, numel(radius));

% For each radius
for i = 1:numel(radius)
    curRadius = radius(i);
    
    % Create a disk meshgrid
    [u,v] = meshgrid(-curRadius:curRadius,-curRadius:curRadius);
    gamma = atan2(v,u);
    mask = (u.^2 + v.^2 <= curRadius^2);
    mask(curRadius+1,curRadius+1) = 0; % mask out center pixel to remove bias

    % Create masks for each orientation
    for j = 1:numel(orientations)
        curOrient = orientations(j);
        side = 1 + (mod(gamma-curOrient,2*pi) < pi);
        side = side .* mask;
        if sum(sum(side==1)) ~= sum(sum(side==2)), error(0, 'bug:inbalance'); end
        orientedMasks{j,1,i} = (side==1);
        orientedMasks{j,1,i} = orientedMasks{j,1,i} / sum(orientedMasks{j,1,i}(:));
        orientedMasks{j,2,i} = (side==2);
        orientedMasks{j,2,i} = orientedMasks{j,2,i} / sum(orientedMasks{j,2,i}(:));
    end
end