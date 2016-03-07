function [ power ] = BrukerRefPower( dataDir,start,varargin )
%BRUKERREFPOWER Processes single pulse scans to get a reference power
%   Detailed explanation goes here
import Bruker.*
p = inputParser;
p.addParameter('verbose',false)
% Load in folder Names
listing = dir(dataDir);
% Pull out directories that matter
j = 1;
for i = 1:length(listing)
    tmp = str2double(listing(i).name);
    if ~isnan(tmp) && tmp>=start
        dataDirs(j) = tmp;
        j = j+1;
    end
end
signals = zeros(size(dataDirs));
powers = zeros(size(dataDirs));
spectrum = zeros(length(dataDirs),2048);
% Load Data
for i = 1:length(dataDirs)
    % If there is no method file break out
    if ~exist(fullfile(dataDir,num2str(dataDirs(i)),'method'), 'file') == 2
        continue
    end
    % Try to load data, if failed move on
    try
        [FIDs, header] = readBrukerReadOut(dataDir, dataDirs(i));
    catch
        continue
    end
    % Make sure the moethod is a single pulse
    if ~any(strcmp(header.Method,{'<Bruker:NSPECT>'}))
        continue
    end
    % get Power
    tmp = strsplit(header.ExcPulse1,',');
    powers(i) = str2double(tmp(11));
    iDelay = header.PVM_DigShift;
    spectrum(i,:) = fftshift(fft(FIDs(iDelay:end),2048));
end
% Remove any scans not loaded
[~,I] = find(powers ~=0);
powers = powers(I);
signals = signals(I);
spectrum = spectrum(I,:);
% Sort  based on powers
[powers,I] = sort(powers);
signals = signals(I);
% Find the center and width of largest peak
spectrum = spectrum(I,:);
tmpLoc = zeros(size(powers));
tmpw = zeros(size(powers));
tmpPks = zeros(size(powers));
    for i = 1:length(powers)
        [pks,locs,w,~] = findpeaks(abs(spectrum(i,:)),'WidthReference','halfheight');
        [~,I] = max(pks);
        tmpPks = pks(I);
        tmpLoc(i) = locs(I);
        tmpw(i) = w(I);
    end
    [~,I] = max(tmpPks);
    peakCenter = tmpLoc(I);
    peakWidth = tmpw(I);
    integralRange = [peakCenter-floor(peakWidth/2):peakCenter+ceil(peakWidth/2)];
% Phase Correct and sum each peak
for i = 1:length(powers)
    % phase correct spectrum, this has no effect as abs is still used
    phase = angle(spectrum(i,peakCenter));
    spectrumCorected = spectrum(i,:).*exp(-1i*phase);
    signals(i) = trapz(integralRange,abs(spectrumCorected(integralRange)));
end
% Assume 50 ohm coil
voltages=sqrt(50*powers);   
%% Run the fit:  
ydata=signals;
xdata=voltages;
[maxv, maxi] = max(ydata);
%assume TR>>T1
%define vector x(1) = amplitude, x(2) = power level for 90 pulse.
fitFunction = @(x,xdata) x(1)*abs(sin(90/x(2)*(pi/180.*xdata)));
opts = optimset('Display','off');
[x,~]=lsqcurvefit(fitFunction,[maxv xdata(maxi)],xdata,ydata,[],[],opts);
refPower = x(2)*x(2)/50;
plot(powers,signals,'kx',powers,fitFunction(x,xdata))
xlabel('Powers'), ylabel('Signal')
fprintf('Calculated reference power: %2.3f Watts \n',refPower);
end

