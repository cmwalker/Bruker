function [ FIDs, Images, xAxis, yAxis, header, Sinogram ] =...
    readRadFLASH( studyDirectory, scanNo)
%READBRUKERRARE Takes the study directory (string) and the scan number (integer) of a RARE scan and
%returns the reshaped raw fids [ReadOut,PhaseEncode,Slices] and the images
%[X,Y,Slices]
%   Detailed explanation goes here
import Bruker.*
% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method,{'<User:radFLASH>'})))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
nProjections = header.PVM_EncMatrix(2); % number of phase encodes
nDummyScans = header.PVM_DummyScans;
projAngle = header.MDA_Pra1;
projAngles = (nDummyScans:(nProjections-1+nDummyScans))*projAngle;
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
% projAngles = header.PVM_EncValues1.*2*pi; % Projection Angles
FOV = header.PVM_Fov;
tmp = reshape(inFIDS,nPoints,nSlices,nProjections);
Sinogram = zeros(nPoints,nProjections,nSlices);
Images = zeros(nPoints,nPoints,nSlices);
FIDs = zeros(nPoints,nProjections,nSlices);
% rearange to have slices be the outer loop
for i = 1:nSlices
    FIDs(:,:,i) = squeeze(tmp(:,i,:));
    Sinogram(:,:,i) = fftshift(fft(squeeze(FIDs(:,:,i)),[],1),1);
    Images(:,:,i) = iradon(abs(squeeze(Sinogram(:,:,i))),projAngles,...
        'linear','Ram-Lak',1,nPoints);
%     if verbose
%         figure
%         subplot(1,2,1)
%         imagesc(abs(squeeze(Sinogram(:,:,i))))
%         subplot(1,2,2)
%         imagesc(abs(squeeze(Images(:,:,i))))
%     end
end
%% Get X and Y Axis assuming read out is X-direction (probably should not be hard coded
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPoints);
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nProjections);
end

