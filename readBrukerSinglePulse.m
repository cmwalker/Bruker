function [FIDs, Spectrum, tAxis, ppmAxis, trAxis, header] = readBrukerSinglePulse(studyDirectory, scanNo)
%READBRUKERSINGLEPULSE Summary of this function goes here
%   Detailed explanation goes here
import Bruker.*

% use subfunction to read Raw Bruker data
[inFIDs, header] = readBrukerReadOut(studyDirectory, scanNo);
% Check header method
if( ~any(strcmp(header.Method,{'<Bruker:NSPECT>'})))
    error('Header not for Bruker NSpect method. Method is %s\n',header.Method);
end
nPoints = header.PVM_DigNp; % number of readout points
filterPoints = header.PVM_DigShift; %Digitizer filer points that must be removed
nTrs = header.PVM_NRepetitions; % get number of repetitions
Tr = header.PVM_RepetitionTime/1000; % get repetition time in seconds
trAxis = 0:Tr:Tr*(nTrs-1); % calculate repetition axis
FIDs = reshape(inFIDs,nPoints,nTrs); % reshape fid for repetitions
FIDs = circshift(FIDs,-filterPoints,1); % Shift fid for digial filtering
FIDs(end-filterPoints:end,:) = 0; % set filter points to zero
Spectrum = fftshift(fft(squeeze(FIDs),[],1),1); % Calculate the spectrum
tAxis = linspace(0,header.PVM_DigDur-header.PVM_DigDw,nPoints); % calculate the time axis for each FID
% Calculate the ppm axis
ppmAxis = linspace(0,header.PVM_SpecSW,nPoints);
ppmAxis = ppmAxis-mean(ppmAxis);
ppmAxis = ppmAxis+header.PVM_FrqWorkOffsetPpm(1);
end



