function [imsh] = readbrukerhos(filn,flg)

echo on
%Usage:
% function [ims imsh imr] = readbrukerh(filn,flg)
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
%
echo off
debugf=0; % 1 for verbose, 0 for silent

if nargin<1
    if ispc
        filn=uigetdir('c:\scans\lucy\data','Select input directory')
    else
        filn=uigetdir('/scans/Lucy/data','Select input directory')
    end
end

if nargin<2
    flg=1;
    reconum=1;
end

magflag=0;
specflag=1;

if strcmp(filn(end:end),'/')|strcmp(filn(end:end),'\')
    filn=filn(1:end-1);
end

tmp=length(filn);
sindx=[1 1 1 1]*tmp;
sindxi=1;
notdone=1;
while (sindxi<5)&(notdone)                          %find the last 4 path separators
    %fprintf('%d --> %s\n',tmp,filn(tmp:tmp));
    if strcmp(filn(tmp:tmp),'/')|strcmp(filn(tmp:tmp),'\')
        sindx(sindxi)=tmp;
        switch sindxi
            case 1,
                if strcmp(filn(sindx(1)-4),'.')      %scan directory name was chosen
                    indir=filn(1:sindx(1)-1);
                    scanno=str2num(filn(sindx(1)+1:end));
                    procno=1;
                    notdone=0;
                    mark = 'scan directory name was chosen'
                end
            case 2,
                if strcmp(filn(sindx(2)-4),'.')      % FID file was chosen
                    indir=filn(1:sindx(2)-1);
                    scanno=str2num(filn(sindx(2)+1:sindx(1)-1));
                    procno=1;
                    flg=0;
                    notdone=0;
                    mark = 'FID file was chosen'
                end
            case 3,
                if strcmp(filn(sindx(3)-4),'.')      % procno directory was chosen
                    indir=filn(1:sindx(3)-1);
                    scanno=str2num(filn(sindx(3)+1:sindx(2)-1));
                    procno=str2num(filn(sindx(1)+1:end));
                    notdone=0;
                    mark = 'procno directory was chosen'
                end
            case 4,
                if strcmp(filn(sindx(4)-4),'.')      % 2dseq file was chosen
                    indir=filn(1:sindx(4)-1);
                    scanno=str2num(filn(sindx(4)+1:sindx(3)-1));
                    procno=str2num(filn(sindx(2)+1:sindx(1)-1));
                    notdone=0;
                    mark = '2dseq file was chosen'
                end
            otherwise,
                error('Unrecognizeable input path/file!')
        end
        sindxi=sindxi+1;
    end
    tmp=tmp-1;
end

%%First, let's read in the headers so we know what's coming:
%The header values that we want to fill in:
imsh=struct('subj',char(ones(1,64)*32),'scan',char(ones(1,64)*32),'method',char(ones(1,64)*32),...
    'indir',indir,'typflag',0,'scanno',scanno,'procno',procno,'te',1.1,'tr',40,'acqsz',[128 128],...
    'phasfact',1,'npe',128,'nop',128,'nsl',2,'BF1',123.456789,...
    'ne',1,'nr',2,'nd',2,'fovpe',3,'fovro',4,'st',1,'si',1.25,'sw',50000.0,'invdel',0.0);

%imsh=struct('te',1.1,'tr',40,'scan',char(ones(1,64)*32),'npe',96,...
%    'nop',128,'nsl',2,'nr',2,'subj',char(ones(1,64)*32),'nd',2,...
%    'fovpe',3,'fovro',4,'st',1,'si',1.25,'indir',sprintf('%s/%d/pdata/%d',indir,scanno,procno));

%First let's get the SUBJECT name:
dirlen=length(indir);
if ispc
    ifil=sprintf('%s\\subject',indir)
else
    ifil=sprintf('%s/subject',indir)
end
ifilh=fopen(ifil,'rt');
if (ifilh<1)
    error('Cannot open subject file!')
end   
eoflag=1;
while eoflag,
    instr=fscanf(ifilh,'%s',1);
    if isempty(instr)
        error('Error reading subject file!')
    end
    if strncmp(instr,'##$SUBJECT_study_name',21)
        tmp=fscanf(ifilh,'%s',3);
        fi=findstr(tmp,'<');
        imsh.subj=tmp(fi+1:end-1);
        eoflag=0;
    end
end
fclose(ifilh);

%Next, info from ACQP file:
if ispc
    ifil=sprintf('%s\\%d\\acqp',indir,scanno);
else
    ifil=sprintf('%s/%d/acqp',indir,scanno);
end
ifilh=fopen(ifil,'rt');
if (ifilh<1)
    error('Cannot open ACQP file!')
end   
eoflag=4+64+128+4096;eoflago=eoflag;
while (eoflag<131071),
    [instr,count]=fscanf(ifilh,'%s',1);
    if count==0
        error('End of ACQP Reached without success!')
    end
    if strncmp(instr,'##$BF1=',7)
        imsh.BF1=str2num(instr(8:end));
        eoflag=bitor(eoflag,65536);
    end
    if strncmp(instr,'##$ACQ_phase_factor=',20)
        imsh.phasfact=str2num(instr(21:end));
        eoflag=bitor(eoflag,32768);
    end     
    if strncmp(instr,'##$ACQ_method',13)
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.method=tmp(2:end-1);
        eoflag=bitor(eoflag,16384);
    end  
    if strncmp(instr,'##$ACQ_fov',10);
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.fovro=str2num(tmp);
        tmp=fscanf(ifilh,'%s',1);
        imsh.fovpe=str2num(tmp);
        eoflag=bitor(eoflag,1);
    end
    if strncmp(instr,'##$ACQ_size',11);
        tmp=fscanf(ifilh,'%s',1);
        imsh.nd=str2num(tmp);
        tmp=fscanf(ifilh,'%s',1);
        tmp=fscanf(ifilh,'%s',1);
        imsh.nop=str2num(tmp)/2;
        imsh.acqsz(1)=imsh.nop;
        if (imsh.nd==1)
            imsh.npe=1;
            imsh.acqsz(2)=1;
        end
        if (imsh.nd>1)
            tmp=fscanf(ifilh,'%s',1);
            imsh.npe=str2num(tmp);
            imsh.acqsz(2)=imsh.npe;
        end
        if (imsh.nd==3)
            tmp=fscanf(ifilh,'%s',1);
            imsh.nsl=str2num(tmp);
            imsh.acqsz(3)=imsh.nsl;
        end
        eoflag=bitor(eoflag,2);
    end
    if strncmp(instr,'##$ACQ_echo_time',16)
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.te=str2num(tmp);
        eoflag=bitor(eoflag,4);
    end
    if strncmp(instr,'##$NI',5)
        if imsh.nd~=3
             imsh.nsl=str2num(instr(7:end));
        end
        eoflag=bitor(eoflag,8);
    end
    if strncmp(instr,'##$NR',5)
        imsh.nr=str2num(instr(7:end));
        eoflag=bitor(eoflag,16);
    end
    if strncmp(instr,'##$ACQ_repetition_time',22)
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.tr=str2num(tmp);
        eoflag=bitor(eoflag,32);
    end
    if strncmp(instr,'##$ACQ_slice_thick',18)
        imsh.st=str2num(instr(20:end));
        eoflag=bitor(eoflag,64);
        if (imsh.st<1e-4)      % This is probably spectro
            imsh.si=0;
            eoflag=bitor(eoflag,128);
            specflag=1;
        end
    end
    if strncmp(instr,'##$ACQ_slice_sepn=',18)
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.si=str2num(tmp);
        eoflag=bitor(eoflag,128);
    end
    if strncmp(instr,'##$ACQ_obj_order=',16)
        tmp2=fscanf(ifilh,'%s',1);
        tmp=fscanf(ifilh,'%s',1);
        slord=zeros(1,str2num(tmp2));
        for ii=1:str2num(tmp2)
            tmp=fscanf(ifilh,'%s',1);
            slord(ii)=str2num(tmp);
        end
        eoflag=bitor(eoflag,256);
    end
    if strncmp(instr,'##$ACQ_scan_name=',17)
        tmp=fscanf(ifilh,'%s',2);
        tmp=fscanf(ifilh,'%s',1);
        imsh.scan=tmp(2:end-1);
        eoflag=bitor(eoflag,512);
    end
    if strncmp(instr,'##$SW_h=',8)
        %instr
        imsh.sw=str2num(instr(9:end));
        eoflag=bitor(eoflag,1024);
    end
    if strncmp(instr,'##$BYTORDA=',11)
        %instr
        bof=instr(12:12);
        eoflag=bitor(eoflag,2048);
    end
    if strncmp(instr,'##$NECHOES=',11)
        %instr
        imsh.ne=str2num(instr(12:end));
        eoflag=bitor(eoflag,4096);
    end
    if strncmp(instr,'##$ACQ_inversion_time=',22)
        %instr
        ninv=str2num(fscanf(ifilh,'%s',1));
        tmp=fscanf(ifilh,'%s',1);
        invtimes=zeros(1,ninv);
        ninvc=0;
        while ninvc<ninv,
            ninvc=ninvc+1;
            invtimes(ninvc)=fscanf(ifilh,'%f',1);
        end
        imsh.invdel=invtimes;
        eoflag=bitor(eoflag,8192);
    end %fprintf('%s\n',instr)
    if (eoflag~=eoflago)&debugf
        fprintf('Acquired %6d: %s\n',eoflag-eoflago,instr);
    end
 %   if (eoflag==eoflago)&(strncmp(instr(1:3),'##$',3))&debugf
 %       fprintf('Nothing for: %s\n',instr);
 %   end
    eoflago=eoflag;
end
fclose(ifilh);

%Now we check to see if there is any other special information
%that should be read from other headers:

%Is this a diffusion set?  If so need to read in b-values
if strcmp(imsh.method,'DtiEpi')|strcmp(imsh.method,'DtiStandard')
    if ispc
        ifil=sprintf('%s\\%d\\method',indir,scanno);
    else
        ifil=sprintf('%s/%d/method',indir,scanno);
    end
    ifilh=fopen(ifil,'rt');
    if (ifilh<1)
        error('Cannot open method file!')
    end   
    eoflag=0;eoflago=0;
    while (eoflag<31),
        [instr,count]=fscanf(ifilh,'%s',1);
        if count==0
            error('End of method file reached without success!')
        end
        if strncmp(instr,'##$PVM_DwNDiffDir',17)
            imsh.dwndir=str2num(instr(19:end));
            eoflag=bitor(eoflag,1);
        end
        if strncmp(instr,'##$PVM_DwNDiffExpEach',21)
            imsh.dwnbvals=str2num(instr(23:end));
            eoflag=bitor(eoflag,2);
        end
        if strncmp(instr,'##$PVM_DwBMat=',14)
            tmp=fscanf(ifilh,'%s',1);
            ndwims=str2num(tmp(1:end-1));
            tmp=fscanf(ifilh,'%s',3);
            bmats=zeros(3,3,ndwims);
            for ii=1:ndwims
                bmats(1,1,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(1,2,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(1,3,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(2,1,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(2,2,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(2,3,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(3,1,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(3,2,ii)=str2num(fscanf(ifilh,'%s',1));
                bmats(3,3,ii)=str2num(fscanf(ifilh,'%s',1));
            end
            imsh.dwbmats=bmats;
            eoflag=bitor(eoflag,4);
        end
        if strncmp(instr,'##$PVM_DwEffBval=',17)
            ndwims=str2num(fscanf(ifilh,'%s',1));
            tmp=fscanf(ifilh,'%s',1);
            effbvals=zeros(1,ndwims);
            for ii=1:ndwims
                effbvals(ii)=str2num(fscanf(ifilh,'%s',1));
            end
            imsh.dweffbval=effbvals;
            eoflag=bitor(eoflag,8);
            %imsh.dwnims=ndwims;
        end
        
        if strncmp(instr,'##$PVM_DwGradVec',16)
            tmp=fscanf(ifilh,'%s',1);
            ndwims=str2num(tmp(1:end-1));
            tmp=fscanf(ifilh,'%s',2);
            gvecs=zeros(3,ndwims);
            for ii=1:ndwims
                gvecs(1,ii)=str2num(fscanf(ifilh,'%s',1));
                gvecs(2,ii)=str2num(fscanf(ifilh,'%s',1));
                gvecs(3,ii)=str2num(fscanf(ifilh,'%s',1));
            end
            imsh.dwgradvecs=gvecs;
            eoflag=bitor(eoflag,16);
        end
        if (eoflag~=eoflago)&debugf
            fprintf('Acquired %s\n',instr);
        end
        eoflago=eoflag;
    end
    fclose(ifilh);
end



npe=imsh.npe;
nop=imsh.nop;
nsl=imsh.nsl;
nr=imsh.nr;
%fprintf('------after or before -------')
% If we're supposed to read the 2dseq (reconstructed) file, we need to
% check the d3proc/reco headers...
if flg==1
    if ispc
        ifil=sprintf('%s\\%d\\pdata\\%d\\reco',indir,scanno,procno);
    else
        ifil=sprintf('%s/%d/pdata/%d/reco',indir,scanno,procno);
    end 
    ifilh=fopen(ifil,'rt');
    if (ifilh<1)
        error('Cannot open RECO file!')
    end   
    eoflag=0;eoflago=0;
    while (eoflag<15),
        [instr,count]=fscanf(ifilh,'%s',1);
        if count==0
            error('End of RECO Reached without success!')
        end
        if strncmp(instr,'##$RECO_size',12);
            tmp=fscanf(ifilh,'%s',1);
            imsh.nd=str2num(tmp);
            tmp=fscanf(ifilh,'%s',1);
            tmp=fscanf(ifilh,'%s',1);
            imsh.nop=str2num(tmp);
            if imsh.nd==1
                imsh.npe=1;
            end
            if imsh.nd>1
                tmp=fscanf(ifilh,'%s',1);
                imsh.npe=str2num(tmp);
            end
            if imsh.nd==3
                tmp=fscanf(ifilh,'%s',1);                
                imsh.nsl=str2num(tmp);       
            end
            eoflag=bitor(eoflag,1);
        end
        if strncmp(instr,'##$RECO_fov',11)
            tmp=fscanf(ifilh,'%s',1);
            imsh.nd=str2num(tmp);
            tmp=fscanf(ifilh,'%s',1);
            tmp=fscanf(ifilh,'%s',1);
            imsh.fovro=str2num(tmp)*10;
            if (imsh.nd>1)
                tmp=fscanf(ifilh,'%s',1);
                imsh.fovpe=str2num(tmp)*10;
            end
            if (imsh.nd==3)
                tmp=fscanf(ifilh,'%s',1);
                imsh.st=str2num(tmp)*10;
            end
            eoflag=bitor(eoflag,2);
        end
        if strncmp(instr,'##$RECO_image_type',18)
            if strncmp(instr(20:end),'MAGNITUDE_IMAGE',15)
                magflag=1;
            else
                magflag=0;
            end
            eoflag=bitor(eoflag,4);
        end
        if strncmp(instr,'##$RECO_wordtype',16)
            if strncmp(instr(18:end),'_32BIT_SGN_INT',14)
                bytespixel=4;
            elseif strncmp(instr(18:end),'_16BIT_SGN_INT',14)
                bytespixel=2;
            end
            eoflag=bitor(eoflag,8);
        end
        %fprintf('%s\n',instr)
        if (eoflag~=eoflago)&debugf
            fprintf('Acquired %s\n',instr);
        end
        eoflago=eoflag;
    end
    fclose(ifilh);
    imsh.st=imsh.st/imsh.nsl;
end

imsh.typflg = flg*4+magflag*2+specflag;