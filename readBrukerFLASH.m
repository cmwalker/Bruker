function [ FIDs, Images, xAxis, yAxis, header ] = readBrukerFLASH( studyDirectory, scanNo )
%READBRUKERRARE Takes the study directory (string) and the scan number (integer) of a RARE scan and
%returns the reshaped raw fids [ReadOut,PhaseEncode,Slices] and the images
%[X,Y,Slices]
%   Detailed explanation goes here
import Bruker.*
% Load in header
header = Bruker.readBrukerHeader(fullfile(studyDirectory,num2str(scanNo),'method'));
% Check header method
if( ~all(header.Method == '<Bruker:FLASH>'))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
minPoints = 128; % minimum number of digitizer points,
nPhaseEncodes = header.PVM_EncMatrix(2); % number of phase encodes
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
FOV = header.PVM_Fov;
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
%% Reshape the FIDS
tmp = reshape(complexFID,tmpNPoints,nSlices,nPhaseEncodes);
tmp = tmp(1:nPoints,:,:);
Images = zeros(nPoints,nPhaseEncodes,nSlices);
FIDs = zeros(nPoints,nPhaseEncodes,nSlices);
% rearange to have slices be the outer loop
for i = 1:nSlices
    FIDs(:,:,i) = squeeze(tmp(:,i,:));
    Images(:,:,i) = abs(fftshift(fftshift(fft2(squeeze(FIDs(:,:,i))),1),2));
end
%% Get X and Y Axis assuming read out is X-direction (probably should not be hard coded
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPoints);
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nPhaseEncodes);
end

