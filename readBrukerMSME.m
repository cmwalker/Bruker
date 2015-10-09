function [ FIDs, Images, xAxis, yAxis, header ] = readBrukerMSME( studyDirectory, scanNo )
%READBRUKERRARE Takes the study directory (string) and the scan number (integer) of a RARE scan and
%returns the reshaped raw fids [ReadOut,PhaseEncode,Slices] and the images
%[X,Y,Slices]
%   Detailed explanation goes here
import Bruker.*

% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~all(header.Method == '<Bruker:MSME>'))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
nPhaseEncodes = header.PVM_EncMatrix(2); % number of phase encodes
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
FOV = header.PVM_Fov;
tmp = reshape(inFIDS,nPoints,nSlices,nPhaseEncodes);
Images = zeros(nPhaseEncodes,nPoints,nSlices);
FIDs = zeros(nPoints,nPhaseEncodes,nSlices);
% rearange to have slices be the outer loop
for i = 1:nSlices
    FIDs(:,:,i) = squeeze(tmp(:,i,:));
    Images(:,:,i) = imrotate(abs(fftshift(fftshift(fft2(squeeze(FIDs(:,:,i))),1),2)),-90);
end
%% Get X and Y Axis assuming read out is X-direction (probably should not be hard coded
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPoints);
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nPhaseEncodes);
end

