function readBrukerStudy( studyDir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import Bruker.*
subFolders = dir(studyDir);
scans = zeros(numel(subFolders),1);
for i =1:numel(subFolders)
    scans(i) = str2double(subFolders(i).('name'));
end
scans(isnan(scans))=[];
scans = sort(scans);
for i = 1:numel(scans)
    if exist(fullfile(studyDir,num2str(scans(i)),'method'), 'file') == 2
        tmpHeader = readBrukerHeader(fullfile(studyDir,num2str(scans(i)),'method'));
    	fprintf('Scan No: %d is a %s Scan.\n',scans(i),tmpHeader.Method);
    else
        fprintf('Scan No: %d has no header.\n',scans(i));
    end
end
end

