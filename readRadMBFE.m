function [FIDs, Images, xAxisPerBand, yAxisPerBand, header,projections] = readRadMBFE(studyDirectory, scanNo,varargin)
%READRADMBFE Summary of this function goes here
%   Detailed explanation goes here
import Bruker.*
%% Parse Input
p = inputParser;
addParameter(p,'ppmCenters',0);
addParameter(p,'projectionOffset',0)
addParameter(p,'pointsPerBand',16)
addParameter(p,'verbose',false)
addParameter(p,'params',{})
parse(p,varargin{:})
verbose = p.Results.verbose;
ppmCenters = p.Results.ppmCenters;
pointsPerBand = p.Results.pointsPerBand;
%% Read Bruker Data
% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method,{'<User:radFLASH>','<User:radFLASH2>'})))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
% Read header
nPoints = header.PVM_Matrix(1);
nBands = nPoints/pointsPerBand;
nProjections = header.PVM_Matrix(2);
nSlices = sum(header.PVM_SPackArrNSlices);
projAngle = header.MDA_Pra1;
iZeroProjection = header.MDA_ProZero+1;
readBandwidth = header.PVM_DigSw;
readBandwidthPPM = readBandwidth/75;
PVM_RepetitionTime = header.PVM_RepetitionTime;
%FOV = header.PVM_Fov;
ppmOffset = header.PVM_FrqWorkOffsetPpm(1);
%readOutTime = (0:(nPoints-1))*(1/readBandwidth);
projAngles = (0:(nProjections-1))*projAngle;
reshapedFIDs = reshape(inFIDS,nPoints,nSlices,nProjections);
ppmBandwidth = header.PVM_EffSWh/75;
%pixelShift = 2*pi/nPoints;
ppmShift = 2*pi/ppmBandwidth;
xAxisPerBand= linspace(-header.PVM_Fov(1)/nBands/2,header.PVM_Fov(1)/nBands/2,pointsPerBand);
yAxisPerBand = xAxisPerBand;
%% Calculate Offsets
offsets = zeros(size(ppmCenters));
BWperBand = ppmBandwidth/nBands;
minPpm = ppmOffset-ppmBandwidth/2;
for i = 1:length(ppmCenters)
    offsets(i) = ppmCenters(i)-minPpm-BWperBand/2;
    bandPpmRange(i,:) = linspace(ppmCenters(i)-BWperBand/2,...
        ppmCenters(i)+BWperBand/2,pointsPerBand);
end
for j = 1:nSlices
FIDs = squeeze(reshapedFIDs(:,j,:));
projections = fftshift(fft(FIDs,[],1),1);
if iZeroProjection > 0
    spectrum = projections(:,iZeroProjection);
    projections(:,iZeroProjection) = [];
end
if verbose
    figure('Position',[634 102 632 824],'Name','Raw Projections')
    xAxis = linspace(-1,1,nPoints);
    ppmAxis = linspace(-1,1,nPoints)*readBandwidthPPM/2+ppmOffset;
    timeAXis = (1:nProjections)*PVM_RepetitionTime/1000;
    subplot(2,1,1),imagesc(timeAXis,ppmAxis,abs(projections));
    ylabel('ppm shifts')
    xlabel('Time (sec)')
    subplot(2,1,2),imagesc(timeAXis,xAxis,abs(projections));
    ylabel('position (pixels)')
    xlabel('Time (sec)')
end
for i = 1:length(offsets)
    shiftedFID = FIDs.*...
        (repmat(exp(-1i*ppmShift*offsets(i)*(1:nPoints)),[nProjections,1]).');
    tmpProj = fftshift(fft(shiftedFID,[],1),1);
    tmpProj = tmpProj(1:pointsPerBand,:);
    Images(:,:,i,j) = iradon(abs(tmpProj),projAngles,'linear','Ram-Lak',1,pointsPerBand);
    if verbose
    figure('Name',sprintf('Band # %d',i),'position',[634 422 1071 503])
    subplot(1,2,1)
    imagesc(timeAXis,bandPpmRange(i,:),abs(tmpProj))
    subplot(1,2,2)
    imagesc(abs(Images(:,:,i,j)));
    end
end
end
end

