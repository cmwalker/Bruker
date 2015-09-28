function [ims, out2] = rbdo2(ifil,imsh,flg,cflg,slrf,f3d)

echo on
% Function [ims,out2] = readbruker(ifil,imsh,flg,dflg,slrf,f3d);
%
% Where:
%
% ifil		is the input filename
% imsh      header
% flg    data mode flag: 1 = 4-byte integer MAGNITUDE (default)
%                        2 = complex data
%                        3 = 2-byte integer MAGNITUDE
%                        4 = raw data
%                        5 = 1-byte integer magnitude
% dflg  data type flag: 'b'=big endian (SGI), 'l'=little endian (Linux)
% slrf   slice reorder flag: set to 1 to reorder 2-d multislice raw data.
% f3d    3-d flag.  Set to 1 for 3d.
%
% Output:
% [ims]	are the [NPE x NOP x NPE2 x npts] images
% [out2] are raw images, if applicable.
%
%Mods:
% 10-31-2002: Corrected bugs in reading of complex data, 
%             added processing for raw... JAB
%  6-12-2003: Added f3d and ability to read 3D raw data.
% 12-31-2003: Added dflg to select big/little endian (SGI or Linux)
%  3-31-2004: Added check to see if raw data is used before it's read,
%             added limit to how much data is read at once before handled
%  8-19-2004: Fixed bug that set raw data bayond first slice to zeros.
%  1-25-2005: Renamed file to separate reading of header and data 
%
echo off

debugf=0;   % 1=verbose, 0=silent;
max_ntp=16; % Max # of timepoints to read at a time.

if ~exist('cflg','var')
    cflg='l';
end

NPE=imsh.npe;
NOP=imsh.nop;
NPE2=imsh.nsl;
NTP=imsh.nr;
if imsh.shuffled==1
    NCH=sum(bitget(imsh.chs,1:8));
else
    NCH=1;
end

if strcmp(imsh.method,'MGE')
    if debugf
        fprintf('Swapping parms for MGE...\n')
    end
    NPE=imsh.nop;
    NOP=imsh.npe;
    NTP=imsh.nsl;
    NPE2=imsh.ne;
end

if strcmp(imsh.method,'fSPGR')
    if debugf
        fprintf('Swapping parms for fSPGR (temporary fix)...\n')
    end
    NPE=imsh.nop;
    NOP=imsh.npe;
end

d=dir(ifil);
if ~exist('flg','var')
   flg=1;
end
if ((d(1).bytes)/NOP/NPE/NPE2/NTP/NCH==2)&(flg==1)
   flg=3;
end
if ~exist('slrf','var')
   slrf=0;
end
if ~exist('f3d','var')
    f3d=0;
end

% Set the slice order, if flag is set:             
if slrf==1
    if mod(NPE2,2)==0
        slord=(mod([1:NPE2],2)==1).*(1+[1:NPE2])/2 + ...
            (mod([1:NPE2],2)==0).*(NPE2+[1:NPE2])/2;
    else
        slord=(mod([1:NPE2],2)==1).*(1+[1:NPE2])/2 + ...
            (mod([1:NPE2],2)==0).*(1+NPE2+[1:NPE2])/2;
    end 
else
    slord=[1:NPE2];
end %slrf 

if debugf
    fprintf('ifil=%s\n',ifil)
    fprintf('[%d %d %d %d]: flg = %d; cflg = %d ; slrf = %d ; f3d = %d\n',...
        NPE, NOP, NPE2, NTP, flg, cflg, slrf, f3d);
end

% Read the data:
FID=fopen(ifil,'r',cflg);
switch flg
    case 1
        dsizestr='int32';
    case 2
        dsizestr='int16';
    case 3
        dsizestr='int16';
    case 4
        dsizestr='int32';
    case 5
        dsizestr='uint8';
end
flg
if (flg==1)|(flg==3)|(flg==5)       % Image-domain magnitude integer data
    d=zeros(NPE,NOP,NPE2,NTP,NCH,dsizestr);
    size(d);
    for mm=1:NCH
        ntpread=0;
        while ntpread<NTP
            ntp2read=min([max_ntp NTP-ntpread]);       
            [big sz]=fread(FID,NOP*NPE*NPE2*ntp2read,dsizestr);
            for kk=1:ntp2read
                for ii=1:NPE2
                    %Not sure why these can't be consistently used:
                    for jj=1:NOP
                        ind1=1+(jj-1)*NPE+(ii-1)*NOP*NPE+(kk-1)*NOP*NPE*NPE2;
                        ind2=ind1+NPE-1;
                        d(:,jj,ii,kk+ntpread,mm) = big(ind1:ind2);
                    end
%                    for jj=1:NPE
%                        ind1=1+(jj-1)*NOP+(ii-1)*NOP*NPE+(kk-1)*NOP*NPE*NPE2;
%                        ind2=ind1+NOP-1;
%                        d(jj,:,ii,kk+ntpread,mm) = big(ind1:ind2);
%                    end
                end
            end
            ntpread=ntpread+ntp2read;
        end
    end
    out2=[];
end
%imagesc(abs(d(:,:,ceil(NPE2/2),ceil(NTP/2),ceil(NCH/2))))
if flg==2                                    % complex reconstructed data;
    d=zeros(NPE,NOP,NPE2,NTP,'int16');
    ntpread=0;
    while ntpread<NTP
        ntp2read=min([max_ntp NTP-ntpread]);
        [big sz]=fread(FID,NOP*NPE*NPE2*2*ntp2read,'int16');
        for kk=1:ntp2read
            for ii=1:NPE2
                for jj=1:NPE
                    ind1=1+(jj-1)*NOP+(ii-1)*NOP*NPE+(kk-1)*NOP*NPE*NPE2;
                    ind2=ind1+(NOP-1);
                    d(jj,:,ii,kk+ntpread) = big(ind1:ind2)+ ...
                        j*big(ind1+NOP*NPE*NPE2*NTP:ind2+NOP*NPE*NPE2*NTP);
                end
            end
        end
        ntpread=ntpread+ntp2read;
    end
    out2=[];
end

if flg==4                                    % Raw data.
    NOPr = 2^ceil(log(NOP)/log(2));          % Data stored with zeros padded...I think!
%    NOP=NOPr;
    d=zeros(NPE,NOPr,NPE2,NTP);
    if debugf
        fprintf('Raw Data: ')
    end
    if ~f3d
        if debugf
            fprintf('2-D\n')
        end
        ntpread=0;
        slrf=0;
        if nargout==2
            out2=d;
        else
            out=0;
        end
        while ntpread<NTP
            ntp2read=min([max_ntp NTP-ntpread]);
            %[big]=fread(FID,2*NOP*NPE*NPE2*NTP,'int32');
            [big]=fread(FID,2*NOPr*NPE*NPE2*ntp2read,'int32');
            outt=zeros(NPE,NOPr,NPE2,ntp2read);
            for kk=1:ntp2read
                for ii=1:NPE2
                    for jj=1:NPE
                        ind1=1+(ii-1)*2*NOPr+(jj-1)*2*NOPr*NPE2+(kk-1)*2*NOPr*NPE*NPE2;
                        ind2=ind1+2*NOPr-1;
%                        [ii jj ind1 ind2 length(big)]
                        outt(jj,:,slord(ii),kk) = big(ind1:2:ind2)+j*big(ind1+1:2:ind2+1);
                    end %for jj
                end %%ii
            end % kk
            for kk=1:ntp2read
                for ii=1:NPE2
                    tmp=squeeze(outt(:,:,slord(ii),kk));
                    if nargout==2
                        out2(:,:,slord(ii),kk+ntpread)=tmp;
                    end
                    d(:,:,slord(ii),kk+ntpread)=fftshift(fft2(fftshift(tmp)));
                end %NPE2
            end %NTP
            ntpread=ntpread+ntp2read;
        end
    else % if ~f3d
        if debugf
            fprintf('3-D Selected\n')
        end
        ntpread=0;
        slrf=0;
        while ntpread<NTP
            outt=zeros(NOP,NPE,NPE2,ntp2read);
            ntp2read=min([max_ntp NTP-ntpread]);
            %[big]=fread(FID,2*NOP*NPE*NPE2*NTP,'int32');
            [big]=fread(FID,2*NOP*NPE*NPE2*ntp2read,'int32');
            if nargout==2
                out2=d;
            else
                out2=0;
            end
            for kk=1:ntp2read
                for ii=1:NPE2
                    for jj=1:NPE
                        ind1=1+(jj-1)*2*NOP+(ii-1)*2*NOP*NPE+(kk-1)*2*NOP*NPE*NPE2;
                        ind2=ind1+2*NOP-1;
                        outt(jj,:,ii,kk) = big(ind1:2:ind2)+j*big(ind1+1:2:ind2+1);
                    end %for jj
                end %%ii
            end % kk
            for kk=1:ntpread
                if nargout==2
                    out2(:,:,:,kk+ntpread)=outt(:,:,:,kk);
                end
                tmp12=squeeze(fftshift(fft(fftshift(outt(:,:,:,kk),3),[],3),3));
                tmp12=fftshift(fft(fftshift(tmp12,2),[],2),2);
                d(:,:,:,kk)=fftshift(fft(fftshift(tmp12,1),[],1),1);
            end %NTP
            ntpread=ntpread+ntp2read;
        end % while ntpread
    end % if ~f3d
end
%ims = squeeze(d);

if strcmp(imsh.method,'MGE')
    if debugf
        fprintf('Reshaping matrix for MGE...\n')
    end
    tmp=d;
    [dim1 dim2 dim3 dim4] = size(d);
    d=zeros(dim1,dim2,dim4,dim3);
    for ii=1:dim3
        for jj=1:dim4
            d(:,:,jj,ii)=tmp(:,:,ii,jj);
        end
    end
end
if strcmp(imsh.method,'fSPGR')
    if debugf
        fprintf('Transposing matrix for fSPGR...\n')
    end
    tmp=d;
    [dim1 dim2 dim3 dim4 dim5] = size(d);
    d=zeros(dim2,dim1,dim3,dim4,dim5);
    for ii=1:dim3
        for jj=1:dim4
            for ijk=1:dim5
                d(:,:,ii,jj,ijk)=squeeze(tmp(:,:,ii,jj,ijk)).';
            end
        end
    end
end
ims=d;
fclose(FID);
