function [FIDs, spects, xAxis, yAxis, ppmAxis, header] = readBrukerCSI(studyDirectory, scanNo)
%READBRUKERCSI Summary of this function goes here
%   Detailed explanation goes here
import Bruker.*
% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method,{'<Bruker:CSI>'})))
    error('Header not for Bruker CSI method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
filterPoints = header.PVM_DigShift; %Digitizer filer points that must be removed
nPhaseEncodes = header.PVM_EncMatrix; % number of phase encodes
FOV = header.PVM_Fov;
ppmAxis = linspace(0,header.PVM_SpecSW,nPoints);
ppmAxis = ppmAxis-mean(ppmAxis);
ppmAxis = ppmAxis+header.PVM_FrqWorkOffsetPpm(1);
FIDs = circshift(inFIDS,-filterPoints,1); % Shift fid for digial filtering
FIDs(end-filterPoints:end,:) = 0; % set filter points to zero
FIDs = reshape(FIDs,[nPoints,nPhaseEncodes(1),nPhaseEncodes(2)]);
spects = fftshift(fft(FIDs,[],1),1);
%% Get X and Y Axis assuming read out is X-direction (probably should not be hard coded
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPhaseEncodes(1));
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nPhaseEncodes(2));
end

