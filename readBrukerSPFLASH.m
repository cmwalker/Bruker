function [ FIDs, Spectrums, freqAxis,tAxis, header ] = readBrukerSPFLASH( studyDirectory, scanNo, varargin )
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
if( ~any(strcmp(header.Method, {'sPFLASH','<User:spFLASH>'})))
    error('Header not for Bruker Rare method. method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
nPhaseEncodes = header.PVM_EncMatrix(2); % number of phase encodes
nSlices = sum(header.PVM_SPackArrNSlices); % number of slices
tmp = reshape(inFIDS,nPoints,nSlices,nPhaseEncodes);
Spectrums = zeros(nPoints,nPhaseEncodes,nSlices);
FIDs = zeros(nPoints,nPhaseEncodes,nSlices);
%% Get Frequency Axis
if isfield(header,'PVM_FrqRef')
    freqAxis = linspace(-header.PVM_EffSWh/2,header.PVM_EffSWh/2,nPoints)+header.PVM_FrqRef(1);
else
    freqAxis = linspace(-header.PVM_EffSWh/2,header.PVM_EffSWh/2,nPoints);
end
sampleingTimes = linspace(0,1/header.PVM_DigSw*header.PVM_DigNp,header.PVM_DigNp);
tAxis = 0:header.PVM_RepetitionTime:header.PVM_RepetitionTime*(nPhaseEncodes-1);
% rearange to have slices be the outer loop
for i = 1:nSlices
    lineBroadeingWindow = zeros(size(squeeze(tmp(:,i,:))))+1;
    % line bradening
    if(~isempty(lineBroadening))
        lineBroadeingWindow = exp(-(sampleingTimes.*sampleingTimes)...
            *lineBroadening^2/(2*4*log(2))).';
        lineBroadeingWindow = repmat(lineBroadeingWindow,1,...
            length(tAxis));
    end
    FIDs(:,:,i) = squeeze(tmp(:,i,:)).*lineBroadeingWindow;
    Spectrums(:,:,i) = fftshift(fft(squeeze(FIDs(:,:,i)),[],1),1);
end
end

