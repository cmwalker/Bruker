function [ FIDs, Images, xAxis, yAxis, header ] = readBrukerRARE( studyDirectory, scanNo )
%READBRUKERRARE Takes the study directory (string) and the scan number (integer) of a RARE scan and
%returns the reshaped raw fids [ReadOut,PhaseEncode,Slices] and the images
%[X,Y,Slices]
%   Detailed explanation goes here
import Bruker.*
% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method,{'<Bruker:RARE>','RARE'})))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
nPhaseEnchodes = header.PVM_EncMatrix(2); % number of phase encodes
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
rareFactor = header.PVM_RareFactor; % number of phase encode lines per excitation (RARE Factor)
FOV = header.PVM_Fov;

%% Reformat Phase order to match matlab indexing KAM
[~,phaseEnchodeOrder] = sort(header.PVM_EncSteps1);
%% Reformat Slice order to match matlab indexing KAM
[~,sliceReorder] = sort(header.PVM_ObjOrderList);

%% Reshape the FIDS
tmp = reshape(inFIDS,nPoints,rareFactor,nSlices,nPhaseEnchodes/rareFactor);
tmp2 = zeros(nPoints,nPhaseEnchodes,nSlices);
for i = 1:nSlices
    tmp2(:,:,i) = reshape(squeeze(tmp(:,:,i,:)),nPoints,nPhaseEnchodes);
end
% %KAM swap out
% FIDs = tmp2(:,phaseEnchodeOrder,:); % reorder the phase encode data
% %swap in
FIDs = tmp2(:,phaseEnchodeOrder,sliceReorder); % reorder the phase and slice
% Fourier transform the Data to get images that also need a -90 degree rotation
Images = imrotate(fftshift(fftshift(fft2(FIDs),1),2),-90);
%% Get X and Y Axis assuming read out is X-direction (probably should not be hard coded
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPoints);
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nPhaseEnchodes);

end

