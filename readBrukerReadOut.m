function [ FIDs, header ] = readBrukerReadOut( studyDirectory, scanNo )
%READBRUKERREADOUT Takes the study directory (string) and the scan number (integer)
%of any Bruker data set and returns the raw FID reshaped [nPoints, all other
%loops] and the header. This is a general function to be used by most other
%functions that read Bruker data

%   reads in raw bruker data, and adjust it for real and complex valuse as well
%   as reshapes for the number of read ot points, accounting for inherent
%   zerobpadding by the digitizer. currently digitizer minimum value is
%   hardcoded probalby want to make this a optional input value
import Bruker.*
% Load in header
header = readBrukerAllHeaders(studyDirectory,scanNo);
nPoints = header.PVM_DigNp; % number of readout points
minPoints = 128; % minimum number of digitizer points,
%% account for digitizer automatic zero padding
if nPoints < minPoints
    tmpNPoints = minPoints;
else
    tmpNPoints = nPoints;
end
if mod(tmpNPoints,minPoints)~=0
    tmpNPoints = 2^nextpow2(tmpNPoints);
end
%% load in raw data
f = fopen(fullfile(studyDirectory,num2str(scanNo),'fid')); % open the FID file
rawFID = fread(f,inf,'int32'); % read in the raw data
fclose(f); % close the file
complexFID = rawFID(1:2:end)+1i*rawFID(2:2:end); % stich real and imaginary points
remander = numel(complexFID)/tmpNPoints;
%% Reshape the FIDS
tmp = reshape(complexFID,tmpNPoints,remander);
FIDs = tmp(1:nPoints,:);
