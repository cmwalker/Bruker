function [ FIDs, Spectrums, tAxis, ppmAxis, trAxis, header ] = readBrukerSPFLASH( studyDirectory, scanNo, varargin )
%READBRUKESPFLASH Takes the study directory (string) and the scan number (integer) of a SpFLASH scan and
%returns the reshaped raw fids [ReadOut,TRs,Slices] and the Spectrums
%[frequency,TRs,Slices]
%   Detailed explanation goes here
import Bruker.*
%% Parse optional inputs
p = inputParser;
addParameter(p,'lb',[]);
parse(p,varargin{:});
lineBroadening = p.Results.lb*pi;
% use subfunction to read Raw Bruker data
[inFIDS, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method, {'sPFLASH','<User:spFLASH>','<User:spFLASH2>'})))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
nPhaseEncodes = header.PVM_EncMatrix(2); % number of phase encodes
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
nTrs = nPhaseEncodes; % get number of repetitions
Tr = header.PVM_RepetitionTime/1000; % get repetition time in seconds
trAxis = 0:Tr:Tr*(nTrs-1); % calculate repetition axis
tmp = reshape(inFIDS,nPoints,nSlices,nPhaseEncodes);
Spectrums = zeros(nPoints,nPhaseEncodes,nSlices);
FIDs = zeros(nPoints,nPhaseEncodes,nSlices);
%% Get Frequency Axis
if isfield(header,'PVM_FrqRef')
    ppmAxis = linspace(-header.PVM_EffSWh/2,header.PVM_EffSWh/2,nPoints)+header.PVM_FrqWorkOffset(1);
else
    ppmAxis = linspace(-header.PVM_EffSWh/2,header.PVM_EffSWh/2,nPoints);
end
ppmAxis = ppmAxis./75;
sampleingTimes = linspace(0,1/header.PVM_DigSw*header.PVM_DigNp,header.PVM_DigNp);
tAxis = linspace(0,header.PVM_DigDur-header.PVM_DigDw,nPoints); % calculate the time axis for each FID
% rearange to have slices be the outer loop
for i = 1:nSlices
    lineBroadeingWindow = zeros(size(squeeze(tmp(:,i,:))))+1;
    % line bradening
    if(~isempty(lineBroadening))
        lineBroadeingWindow = exp(-(sampleingTimes.*sampleingTimes)...
            *lineBroadening^2/(2*4*log(2))).';
        lineBroadeingWindow = repmat(lineBroadeingWindow,1,...
            length(trAxis));
    end
    FIDs(:,:,i) = squeeze(tmp(:,i,:)).*lineBroadeingWindow;
    Spectrums(:,:,i) = fftshift(fft(squeeze(FIDs(:,:,i)),[],1),1);
end
end

