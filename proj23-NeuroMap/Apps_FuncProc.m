% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Source Code Provided By: Sridhar Kandala
%
% fix_3_clean(fixlist,aggressive,domot,hp) - apply the FIX cleanup to filtered_func_data
% fix_3_clean(fixlist,aggressive,domot,hp,0) - CIFTI processing only (do not process volumetric data)
%
% fixlist is a vector of which ICA components to remove (starting at 1 not 0)
%
% aggressive = 0 or 1 - this controls whether cleanup is aggressive (all variance in confounds) or not (only unique variance)
%
% mot = 0 or 1 - this controls whether to regress motion parameters out of the data (24 regressors)
%
% hp determines what highpass filtering had been applied to the data (and so will get applied to the motion confound parameters)
% hp=-1 no highpass
% hp=0 linear trend removal
% hp>0 the fullwidth (2*sigma) of the highpass, in seconds (not TRs)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a modified script that will allow for processing volumetric data
% from the HCP-A and/or HCP-D data sets. The CIFTI processing portion will
% be disabled so that only volume data is processed (the default behavior
% is to generate CIFTI files, and *optionally also* NIFTI files).
%
% There is one 'major' modification, which is to to use the MATLAB nifti
% functions niftiinfo, niftiread and niftiwrite instead of the older FSL
% read_avw and save_avw (which require fslcpgeom to make a compliant NIFTI
% header
%
% the first input is the path to a .fix file which is the
% timeseries from the ICA components. This script assumes an HCP style
% setup, find the relevant data and motion files, loads them, produces the
% output and then write the output to the second input, the outputdir
% while this script is deployed (compiled) we will keep these options
% hardcoded regardless of interactive or compiled matlab
%
%if (isdeployed)
%    aggressive = str2double(aggressive);
%    domot = str2double(domot);
%    hp = str2double(hp);
%end
%%% setup the following variables for your site

hp = -1; 
DOvol=1; 
domot = 1; 
aggressive = 0; 
[~,whoami] = system('whoami'); whoami = strcat(whoami);
[~,WBC] = system('which wb_command'); WBC = strcat(WBC); FSLDIR = getenv('FSLDIR'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixlist='/ceph/intradb/archive/CinaB/CCF_HCD_STG/SUBID/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_hp0.ica/.fix'
addpath('/scratch/jirsaraie/RobertJirsaraie/toolbox/bids_apps/dependencies/designs_xcp')
outputdir='/scratch/jirsaraie/study-HCD/bids/SUBID'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get subject id name
pathparts = strsplit(fixlist,'/'); % split file path by '/' character
subject = char(pathparts(find(strcmp(pathparts,'MNINonLinear'))-1));
mkdir(outputdir); 

% fileparts will go up to .ica->ConcatDir->ResultsDir directories
[ica,~,~] = fileparts(fixlist); % ica is the directory where fixlist resides
[ConcatDir,fMRI_CONCATvol,~] = fileparts(ica); % concat dir is where the data is, and concatvol is the volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a mask that keeps resting state runs and cuts the task fMRI out

ConcatCsv = textscan(fopen([ConcatDir '/' strrep(fMRI_CONCATvol,'hp0','Runs.csv')]),'%s%f%f','delimiter',',');
RunNames = ConcatCsv{1}(contains(ConcatCsv{1},'rfMRI'));
RunBoundaries = [ConcatCsv{2}(contains(ConcatCsv{1},'rfMRI')) ConcatCsv{3}(contains(ConcatCsv{1},'rfMRI'))];

rmask = zeros(ConcatCsv{3}(end),1);
for i = 1:size(RunBoundaries,1)
    rmask(RunBoundaries(i,1):RunBoundaries(i,2)) = 1;
end
rmask = logical(rmask);

RunLengths = diff(RunBoundaries,[],2)+1; %#ok<*NASGU>

%%%%%%%%%%%%%%%% the following code is from fix_3_clean %%%%%%%%%%%%%%%%%%
%%%%  read set of bad components
DDremove=load(fixlist);

%%%%  unzip data in user /scratch to prevent default /tmp r+w space issue
tmpdir = [outputdir,'/temp']; mkdir(tmpdir);
OG_FILE=[ConcatDir,'/',fMRI_CONCATvol,'.nii.gz']
TMP_FILE=[tmpdir,'/',fMRI_CONCATvol,'.nii.gz']
copyfile(OG_FILE,TMP_FILE);
unzip_command = char(strcat("gunzip -f ",TMP_FILE));
system(unzip_command);
UNZIP_FILE=[tmpdir,'/',fMRI_CONCATvol,'.nii']

%%%%  find TR of data
cts_info = niftiinfo(UNZIP_FILE); 
TR = cts_info.PixelDimensions(4);
sprintf('%s TR=%0.2f',fMRI_CONCATvol,TR); DObrainord=0;

%%%%  read NIFTI version of the data
if DOvol
    if(~exist('cts_info','var'))        
        unzip(strcat(ConcatDir,'/',fMRI_CONCATvol,'.nii.gz'),tmpdir)
        cts_info=niftiinfo(strcat(tmpdir,'/',fMRI_CONCATvol,'.nii'));
    end
  cts = niftiread(cts_info);
  ctsX=size(cts,1); ctsY=size(cts,2); ctsZ=size(cts,3); ctsT=size(cts,4); 
  cts=reshape(cts,ctsX*ctsY*ctsZ,ctsT)';
end

%%%%  read and prepare motion confounds
confounds=[];
if domot == 1
  confounds = niftiread([ica '/mc/prefiltered_func_data_mcf_conf_hp.nii.gz']);
  size(confounds)
  confounds = functionnormalise(squeeze(confounds)');
  size(confounds)
end

%%%%  read ICA component timeseries
ICA=functionnormalise(load(sprintf([ ica '/filtered_func_data.ica/melodic_mix'])));

%%%%  do the cleanup
if aggressive == 1
	sprintf('aggressive cleanup')
	confounds=[confounds ICA(:,DDremove)];
	if DOvol
        cts = cts - (confounds * (pinv(confounds,1e-6) * cts));
	end
	if DObrainord
        BO.cdata = BO.cdata - (confounds * (pinv(confounds,1e-6) * BO.cdata'))'; %#ok<*UNRCH>
	end
else
	sprintf('unaggressive cleanup')
    if domot == 1
        % aggressively regress out motion parameters from ICA and from data
        % HERE - Matrix Multiplation Not Permitted on Line Below!!!
        ICA = ICA - (confounds * (pinv(confounds,1e-6) * ICA));
        if DOvol
            cts = cts - (confounds * (pinv(confounds,1e-6) * cts));
        end
        if DObrainord
            BO.cdata = BO.cdata - (confounds * (pinv(confounds,1e-6) * BO.cdata'))';
        end
    end
    if DOvol
        betaICA = pinv(ICA,1e-6) * cts;                         % beta for ICA (good *and* bad)
        cts = cts - (ICA(:,DDremove) * betaICA(DDremove,:));    % cleanup
    end
    if DObrainord
        betaICA = pinv(ICA,1e-6) * BO.cdata';                              % beta for ICA (good *and* bad)
        BO.cdata = BO.cdata - (ICA(:,DDremove) * betaICA(DDremove,:))';    % cleanup
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load a brain mask to find non-zero voxels for computing global signal
final_mask = niftiread([ConcatDir '/fMRI_CONCAT_ALL_brain_mask.nii.gz']);
final_mask = logical(final_mask(:));
MGTR = zeros(ctsT,2); 

% compute global signal, and regress it out
cts=cts';
MGTR(:,1) = mean(cts(final_mask,:));
MGTR(:,2) = functionnormalise([0;diff(MGTR(:,1))]);
cts = cts - (MGTR * (pinv(MGTR,1e-6) * cts'))';

%%%%%%%%%%%%%%%%%%%%%%%%% highpassfilter the data %%%%%%%%%%%%%%%%%%%%%%%%%
[butter2,butter1] = butter(1,0.009*2*TR,'high');
cts = transpose(single(filtfilt(butter2,butter1,transpose(double(cts)))));
  
%%%% save cleaned data to file
if DOvol
    cts = reshape(cts(:,rmask),[ctsX ctsY ctsZ sum(rmask)]);
    cts_info.ImageSize(4) = sum(rmask);
    outname = [outputdir '/func/' subject '_rfMRI_REST_hp0_clean_mgtr_hpss'];
	niftiwrite(cts,outname,cts_info);
    system(['gzip ' outname '.nii']);
	%call_fsl('fslcpgeom filtered_func_data filtered_func_data_clean');
end
if DObrainord
    B1 = BO;
    BN = BO;
    BO.cdata = BO.cdata(:,rmask);
	ciftisavereset(BO,[outputdir '/func/' subject '_rfMRI_REST_Atlas_hp0_clean_mgtr_hpss.dtseries.nii'],WBC);
    for i = 1:length(RunNames)
        BN.cdata = B1.cdata(:,RunBoundaries(i,1):RunBoundaries(i,2));
        ciftisavereset(BN,[outputdir '/func/' subject '_' RunNames{i} '_Atlas_hp0_clean_mgtr_hpss.dtseries.nii'],WBC);
    end
end

%%%%%%%%%%%%%%%% Compute the Filtered FD and write it out %%%%%%%%%%%%%%%%ƒ
% to compute FD, read in the movement_regressors file
% convert the Pitch (P), Yaw (Ya), and Roll (R) from degrees to arclength 
% assuming r=50mm
% compute the 1st moment difference for the 6 parameters (XYZPYaR), take
% the absolute value, and sum across the parameters to compute FD
% Filtered FD is compute by filtering movement regressors, and computing
% the difference between (t) and (t-4), which most approximates FD in
% single-band data w/ TR ~= 2.5s-3s

MVM = dlmread([ConcatDir '/Movement_Regressors_demean.txt']);
MVM = [MVM(:,1:3) (50*pi/180)*MVM(:,4:6)];

% generate filter parameters and make filtered movement matrix
nyq = (1/TR)/2; 
stopband = [0.2 0.5]; 
Wn = stopband/nyq; 
filtN = 10; 
[B,A] = butter(filtN,Wn,'stop'); 
MVM_filt = filtfilt(B,A,MVM);

% create standard FD
dMVM = diff(MVM,1);
FD = [0;sum(abs(dMVM),2)];

% create the 1-back filtered FD
dMVM_filt = diff(MVM_filt,1);
FD_filt1 = [0;sum(abs(dMVM_filt),2)];

% create the 4-back filtered FD
FD_filt4 = [zeros(4,1);sum(abs(MVM_filt(1+4:end,:) - MVM_filt(1:end-4,:)),2)];

% and write them to disk
dlmwrite([outputdir '/func/FD.txt'],FD(rmask),'precision','%.4f');
dlmwrite([outputdir '/func/FDfilt1.txt'],FD_filt1(rmask),'precision','%.4f');
dlmwrite([outputdir '/func/FDfilt4.txt'],FD_filt4(rmask),'precision','%.4f');

% remove tmpdir
rmdir(tmpdir,'s');

%#######⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
%###           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
%#######⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######