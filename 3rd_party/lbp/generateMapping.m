close all
clear all

%% Open matlabpool
if matlabpool('size') == 0
    matlabpool open
end

%% Calculate mapping
P = 24;
mappingMode = 'ri';
mapping = getmapping_par(P, mappingMode);            % We use rotationally invariant 

%% Save to file
fileName = sprintf('mappings/%s_p%d.mat', mappingMode, P);
save(fileName, 'mapping');