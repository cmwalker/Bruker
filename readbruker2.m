function [ims, imsh, imr] = readbruker2(filn,flg,specflag)

echo on
%Usage:
% function [ims imsh imr] = readbruker(filn,flg)
%
% Input variables:
%  filn     =   input file or directory name
%  flg      =   0 for raw, 1 for image data.
%
% Output variables:
%  ims     = image domain data [npe nop nsl nr]
%  imsh    = image information (header)
%  imr     = raw data, if needed.
%
% JAB May 9,2003
%
% Modifications:
% 8-19-2004:    Fixed bug in object reordering code.
% 1-25-2005:    Separated Header, Data input modules
%
echo off
debugf=0; % 1 for verbose, 0 for silent

switch nargin
    case 0
        imsh=rbho;
    case 1
        imsh=rbho(filn);
    case 2
        imsh=rbho(filn,flg);
end

if ~exist('specflag','var')
    specflag=0;
end

npe=imsh.npe;
nop=imsh.nop;
nsl=imsh.nsl;
nr=imsh.nr;

flg=(bitand(imsh.typflg,4))>0;
magflag=(bitand(imsh.typflg,2))>0;
%specflag=bitand(imsh.typflg,1);

if ispc
    ifil=sprintf('%s\\%d\\acqp',imsh.indir,imsh.scanno);
else
    ifil=sprintf('%s/%d/acqp',imsh.indir,imsh.scanno);
end
ifilh=fopen(ifil,'rt');
if (ifilh<1)
    error('Cannot open ACQP file!')
end   
eoflag=0;eoflago=0;
while (eoflag<3),
    [instr,count]=fscanf(ifilh,'%s',1);
    if count==0
        error('End of ACQP Reached without success!')
    end
    if strncmp(instr,'##$ACQ_obj_order=',16)
        tmp2=fscanf(ifilh,'%s',1);
        tmp=fscanf(ifilh,'%s',1);
        slord=zeros(1,str2num(tmp2));
        for ii=1:str2num(tmp2)
            tmp=fscanf(ifilh,'%s',1);
            slord(ii)=str2num(tmp);
        end
        eoflag=bitor(eoflag,1);
    end
    if strncmp(instr,'##$BYTORDA=',11)
        %instr
        bof=instr(12:12);
        eoflag=bitor(eoflag,3);
    end
    if (eoflag~=eoflago)&debugf
        fprintf('Acquired %s\n',instr);
    end
    eoflago=eoflag;
end
fclose(ifilh);

%fprintf('------after or before -------')
% If we're supposed to read the 2dseq (reconstructed) file, we need to
% check the d3proc/reco headers...
if flg==1
    if ispc
        ifil=sprintf('%s\\%d\\pdata\\%d\\reco',imsh.indir,imsh.scanno,imsh.procno);
        ifil3d=sprintf('%s\\%d\\pdata\\%d\\d3proc',imsh.indir,imsh.scanno,imsh.procno);
    else
        ifil=sprintf('%s/%d/pdata/%d/reco',imsh.indir,imsh.scanno,imsh.procno);
        ifil3d=sprintf('%s/%d/pdata/%d/d3proc',imsh.indir,imsh.scanno,imsh.procno);
    end 
    
    if exist(ifil,'file')
        ifilh=fopen(ifil,'rt');
        if (ifilh<1)
            error('Cannot open RECO file!')
        end
        eoflag=0;eoflago=0;
        while (eoflag<1),
            [instr,count]=fscanf(ifilh,'%s',1);
            if count==0
                error('End of RECO Reached without success!')
            end
            if strncmp(instr,'##$RECO_wordtype',16)
                if strncmp(instr(18:end),'_32BIT_SGN_INT',14)
                    bytespixel=4;
                elseif strncmp(instr(18:end),'_16BIT_SGN_INT',14)
                    bytespixel=2;
                end
                eoflag=bitor(eoflag,1);
            end
            %fprintf('%s\n',instr)
            if (eoflag~=eoflago)&debugf
                fprintf('Acquired %s\n',instr);
            end
            eoflago=eoflag;
        end
        fclose(ifilh);
    elseif exist(ifil3d,'file')
        ifilh=fopen(ifil3d,'rt');
        if (ifilh<1)
            error('Cannot open d3proc file!')
        end
        eoflag=0;eoflago=0;
        while (eoflag<1),
            [instr,count]=fscanf(ifilh,'%s',1);
            if count==0
                error('End of d3proc Reached without success!')
            end
            if strncmp(instr,'##$DATTYPE',10)
                bytespixel=1;
                magflag=1; %Until I can figure out what it's supposed to be...
                if strncmp(instr(15:end),'u_byte',6)
                    bytespixel=1;
                    magflag=1;
                end
                eoflag=bitor(eoflag,1);
            end
            if (eoflag~=eoflago)&debugf
                fprintf('Acquired %s\n',instr);
            end
            eoflago=eoflag;
        end
        fclose(ifilh);
    else
        error('Cannot open reco or d3proc files!')
    end
end

%Now read the data
if flg==0   %Read raw data
    if debugf
        fprintf('readbrukerh: read raw data ');
    end
    if ispc
        ifilnm=sprintf('%s\\%d\\fid',imsh.indir,imsh.scanno');
    else
        ifilnm=sprintf('%s/%d/fid',imsh.indir,imsh.scanno');
    end
    if imsh.nd==2
        if debugf
            fprintf('2-D\n')
        end
        if nargout>2
            [ims imr]=rbdo2(ifilnm,imsh,4,bof,0);
        else
            [ims]=rbdo2(ifilnm,imsh,4,bof,0);
        end
        for slordc=1:imsh.nr
            ims(:,:,slord+1,slordc)=ims(:,:,[1:imsh.nsl],slordc);
            if nargout>2
                imr(:,:,slord+1,slordc)=imr(:,:,[1:imsh.nsl],slordc);
            end
        end
    elseif imsh.nd==3
        if debugf
            fprintf('3-D\n')
        end
        if nargout>2
            [ims imr] = rbdo2(ifilnm,imsh,4,bof,0,1);
        else
            [ims] = rbdo2(ifilnm,imsh,4,bof,0,1);
        end
%        error('3D Raw yet to be implemented...')
    end
end

if ispc
    if specflag
        ifilnm=sprintf('%s\\%d\\pdata\\%d\\1',imsh.indir,imsh.scanno,imsh.procno);
    else
        ifilnm=sprintf('%s\\%d\\pdata\\%d\\2dseq',imsh.indir,imsh.scanno,imsh.procno);
    end
else
    if specflag
        ifilnm=sprintf('%s/%d/pdata/%d/1',imsh.indir,imsh.scanno,imsh.procno);
    else
        ifilnm=sprintf('%s/%d/pdata/%d/2dseq',imsh.indir,imsh.scanno,imsh.procno);
    end
end

if specflag
    imsh.nop=imsh.nop*2;
    ftmp=sprintf('%sr',ifilnm)
    %[imsh.nop imsh.npe imsh.nsl imsh.nr]
    ims=rbdo2(ftmp,imsh,1,bof);
    ftmp=sprintf('%si',ifilnm)
    imstmp=rbdo2(ftmp,imsh,1,bof);
    ims=ims+j*imstmp;
end

if (magflag)&(~specflag) %Read magnitude data, not spectro
    if bytespixel==1
        tmp=5;
    elseif bytespixel==2
        tmp=3;
    elseif bytespixel==4
        tmp=1;
    end
    if debugf
        fprintf('read magnitude image data:rbdo2(fname,imsh,tmp=%d,bof=%s)\n',tmp,bof)
    end
    [ims] = rbdo2(ifilnm,imsh,tmp,bof);
    imr=[];
end

if (~magflag)&(flg) %Not magnitude, not raw...must be phase
    if debugf
        fprintf('readbrukerh: read complex image')
    end
    [ims] = rbdo2(ifilnm,imsh,2,bof);
    imr=[];
end

if debugf
    if specflag
        figure(1)
        xax=linspace(-imsh.sw/2,imsh.sw/2,imsh.nop);
        subplot('position',[0.05 0.5 0.85 0.4])
        plot(xax,real(ims))
        title(sprintf('%s: %s',imsh.subj,imsh.scan));
        axis tight
        ylabel('real')
        subplot('position',[0.05 0.05 0.85 0.4])
        plot(xax,imag(ims))
        axis tight
        ylabel('imag')
    else
        figure(1)
        imagesc([-imsh.fovro/2 imsh.fovro/2], [-imsh.fovpe/2 imsh.fovpe/2],...
            abs(squeeze(ims(:,:,ceil(imsh.nsl/2),1))));
        title(sprintf('%s: %s',imsh.subj,imsh.scan));
        [npe nop nsl nr] = size(ims);
        ylabel(sprintf('NPE = %d points',npe));
        xlabel(sprintf('NOP = %d points',nop));
        colormap(gray)
        shg
    end
end