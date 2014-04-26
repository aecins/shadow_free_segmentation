function paramVal = findParameter(mode, paramName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given input string that containts a substring followed by a number find
% the number i.e. given string '*_weight30_*' return 30.
% 
% Input:
%   mode,           input string
%   paramName,      parameter substring
%
% Output:
%   paramVal,       value of the parameter
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = regexp(mode,'_','split');
matchString = [paramName, '\d+'];
for i=1:numel(tmp)
    matchStr = regexp(tmp{i}, matchString, 'match');
    if ~isempty(matchStr)
        paramVal = str2double(matchStr{1}(length(paramName)+1:end));
        matchStr = regexp(tmp{i}, 'X\d+', 'match');
        for j=1:numel(matchStr)
            paramVal = [paramVal, str2double(matchStr{j}(2:end))];
        end
        return
    end
end

paramVal = [];