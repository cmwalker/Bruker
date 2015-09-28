function [ FIDs, Images ] = readBrukerRARE( studyDirectory, scanNo )
%READBRUKERRARE Takes the study directory (string) and the scan number (integer) of a RARE scan and
%returns the reshaped raw fids [ReadOut,PhaseEncode,Slices] and the images
%[X,Y,Slices]
%   Detailed explanation goes here
import Bruker.*
% Load in header
header = Bruker.readBrukerHeader(fullfile(studyDirectory,num2str(scanNo),'method'));
% Check header method
if( ~all(header.Method == '<Bruker:RARE>'))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
minPoints = 128; % minimum number of digitizer points,
nPhaseEnchodes = header.PVM_Matrix(2); % number of phase encodes
nSlices = header.PVM_SPackArrNSlices; % number of slices
rareFactor = header.PVM_RareFactor; % number of phase encode lines per excitation (RARE Factor)
%% account for digitizer automatic zero padding
if nPoints < minPoints
    tmpNPoints = minPoints;
else
    tmpNPoints = nPoints;
end
if mod(tmpNPoints,minPoints)~=0
    tmpNPoints = 2^nextpow2(tmpNPoints);
end
%% Reformat Phase order to match matlab indexing
phaseEnchodeOrder = header.PVM_EncSteps1+nPhaseEnchodes/2+1; % position of each phase encode lines
tmp = zeros(size(phaseEnchodeOrder));
for i = 1:length(phaseEnchodeOrder)
    tmp(i) = find(i==phaseEnchodeOrder);
end
phaseEnchodeOrder = tmp;
%% load in raw data
f = fopen(fullfile(studyDirectory,num2str(scanNo),'fid')); % open the FID file
rawFID = fread(f,inf,'int32'); % read in the raw data
fclose(f); % close the file
complexFID = rawFID(1:2:end)+1i*rawFID(2:2:end); % stich real and imaginary points
%% Reshape the FIDS
tmp = reshape(complexFID,tmpNPoints,rareFactor,nSlices,nPhaseEnchodes/rareFactor);
tmp = tmp(1:nPoints,:,:,:);
tmp2 = zeros(nPoints,nPhaseEnchodes,nSlices);
for i = 1:nSlices
    tmp2(:,:,i) = reshape(squeeze(tmp(:,:,i,:)),nPoints,nPhaseEnchodes);
end
FIDs = tmp2(:,phaseEnchodeOrder,:); % reorder the phase encode data
Images = fftshift(fftshift(fft2(FIDs),1),2);

figure
colormap gray
for i = 1:size(FIDs,3)
    subplot(ceil(size(FIDs,3)/2),2,i),imagesc(abs(Images(:,:,i)))
end
end

