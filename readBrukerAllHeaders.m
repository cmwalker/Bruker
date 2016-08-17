function [ headers,method,acqp ] = readBrukerAllHeaders( studyDirectory, scanNo )
%READBRUKERALLHEADERS Takes the study directory (string) and the scan number (integer)
%of any Bruker data set and returns a structure containing the header
%information stored in both acqp and method files. both in a large struct
% and seperate structs

%   reads in raw bruker data, and adjust it for real and complex valuse as well
%   as reshapes for the number of read ot points, accounting for inherent
%   zerobpadding by the digitizer. currently digitizer minimum value is
%   hardcoded probalby want to make this a optional input value
import Bruker.*
% Load in headers
headers = struct();
method = readBrukerHeader(fullfile(studyDirectory,num2str(scanNo),'method'));
acqp = readBrukerHeader(fullfile(studyDirectory,num2str(scanNo),'acqp'));
subject = readBrukerHeader(fullfile(studyDirectory,'subject'));
% Fill with subject info
tmpNames = fieldnames(subject);
for i = 1:numel(tmpNames)
    headers.(tmpNames{i}) = subject.(tmpNames{i});   
end
% Fill with method info
tmpNames = fieldnames(method);
duplicateNames = {'TITLE','JCAMPDX','DATATYPE','ORIGIN','OWNER'};
for i = 1:numel(duplicateNames)
    if(isfield(acqp,duplicateNames{i}))
        tmpNames(strcmp(tmpNames,duplicateNames{i})) = [];
    end
end
for i = 1:numel(tmpNames)
    if(isfield(headers,tmpNames{i}))
        warning('*WARNING* %s is a field name in both the method and subject file. Using the value from the method file/n',tmpNames{i})
    end
    headers.(tmpNames{i}) = method.(tmpNames{i});   
end
% Fill with acqp info
tmpNames = fieldnames(acqp);
for i = 1:numel(duplicateNames)
    if(isfield(acqp,duplicateNames{i}))
        tmpNames(strcmp(tmpNames,duplicateNames{i})) = [];
    end
end
for i = 1:numel(tmpNames)
    if(isfield(headers,tmpNames{i}))
        warning('*WARNING* %s is a field name in both the method/subjec and acqp file. Using the value from the acqp file/n',tmpNames{i})
    end
    headers.(tmpNames{i}) = acqp.(tmpNames{i});   
end
end

