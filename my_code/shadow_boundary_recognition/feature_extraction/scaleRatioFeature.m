function data = scaleRatioFeature(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the a feature that is a ratio of two values (i.e. intensity ratio)
% to the [0, 1] range.
% 
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

magLim = 4;

dataNegMask = data < 1;
dataPosMask = data > 1;

data(dataPosMask) = data(dataPosMask) - 1;
data(dataNegMask) = -1 ./ data(dataNegMask);

data(data>magLim) = magLim;
data(data<-magLim) = -magLim;

data = data / magLim;