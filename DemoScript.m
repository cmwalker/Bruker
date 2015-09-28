%% Demonstration script for the functions to read and procces Bruker MRI Data
% Author: Chris Walker
% Date created: 9/28/2015

%% Prep Workspace
clear all
close all
clc
import Bruker.*

%% Read RARE data
studyDirectory = 'DemoData\BrukerRAREData';
scanNo = 1;
[FIDS, Images, xAxis, yAxis] = readBrukerRARE( studyDirectory, scanNo );
figure('Name','RARE Images')
colormap gray
for i = 1:size(Images,3)
    subplot(ceil(size(Images,3)/2),2,i),imagesc(xAxis, yAxis, abs(Images(:,:,i)));
    title(sprintf('Slice %d',i))
end
