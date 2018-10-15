
% Requirements:
% 1) DICOM MR images reconstructed by the MR console are required to obtain MR
% image coordinates
% 
% 2) NIFTI and SPM(12) packages
            
io.spm_path = 'c:\spm12\';
io.nifit_path = 'c:\nifti\';
io.asl_lib_path = 'c:\ASL_Recon\';
io.rawdata_path = 'c:\ASL_Recon\example_dataset\';
io.rawdata_flname = [io.rawdata_path 'meas_MID00165_FID10100_Head_jj_3Dpcasl_1X2pe_kb_1500_1800_BW3126_etd511.dat'];

% add path
addpath(io.spm_path);
addpath(io.nifit_path)
addpath(genpath(io.asl_lib_path))

io.t1mprage.dicom_path = [io.rawdata_path 'dicom_t1mprage' os_bar];
io.t1mprage.nii_path = [io.t1mprage.dicom_path 'nii' os_bar];
io.t1mprage.nii_flname =[];
io.results_path = [io.rawdata_path 'results' os_bar];
io.timeseries.dicom_path = [io.rawdata_path 'dicom_time_series' os_bar];


%% convert mprage dicom into 3D nifit and do segementgation
io.t1mprage.nii_flname = spm_dicom2nifti(io.t1mprage.dicom_path);

% segmentation into Gray matter, white matter and CSF
job = spm_create_seg_job(io);
spm('defaults', 'FMRI');
spm_jobman('run', job, {});

[pathname, fname,ext] = fileparts(io.t1mprage.nii_flname);
io.t1mprage.gm_flname = [pathname os_bar 'c1' fname ext];
io.t1mprage.wm_flname = [pathname os_bar 'c2' fname ext];
io.t1mprage.csf_flname = [pathname os_bar 'c3' fname ext];

gm = load_untouch_nii(io.t1mprage.gm_flname);
wm = load_untouch_nii(io.t1mprage.wm_flname);
csf = load_untouch_nii(io.t1mprage.csf_flname);
brainExtractMask = (gm.img + wm.img + csf.img)>0;
brainExtractMask = logical(imdilate(single(brainExtractMask),rStrel()));
%% Convert time series into 4D nifit
spm_dicom2nifti(io.timeseries.dicom_path);
out_dir = [io.timeseries.dicom_path 'nii' os_bar];
%% registration of M0 to anatomical image
if 0 % in some cases the registeration is not satisfactory, hence better to rely on the imperfect registerations of M0 and T1
    imgs=spm_select('FPList',out_dir, '\w*\.nii$');
    VG = spm_vol(io.t1mprage.nii_flname);
    VF = spm_vol(imgs(1,:)); % first image in the series is asummed to be M0, otherwise it needs to be specfified by the user
    x  = spm_coreg(VG, VF);
    M  = inv(spm_matrix(x));
    spm_get_space(VF.fname,M*VF.mat)
end
%% make 4D nifti
imgs=spm_select('FPList',out_dir, '\w*\.nii$');

io.timeseries.nii_path = [out_dir '4D' os_bar];
io.timeseries.nii_flname =  [io.timeseries.nii_path '4DTimeSeries.nii'];
mkdir(io.timeseries.nii_path)
spm_file_merge(imgs,io.timeseries.nii_flname);
%% estimate affine transfromations
% estimate affine transfromations to referece motion state
% by registeratig all time frames to first one (m0)
reaFlags = struct('quality', 0.9, 'fwhm', 4,'rtm', 1,'PW','');
Pt = spm_select('ExtFPList',io.timeseries.nii_path,'4DTimeSeries.nii');
P = spm_realign_asl_2(Pt, reaFlags);
save([io.timeseries.nii_path 'transformix.mat'],'P')
%% build aslRecon object 
[MrInfo, MrNifti] = getNiftiDataInfo(io.t1mprage.nii_flname);
[tsInfo,tsNifit] = getNiftiDataInfo(io.timeseries.nii_flname);
twix = mapVBVD(io.rawdata_flname);

removeOS = 1;
[Data,hdr,imgOut] = getwixASLdata(twix{1,2},removeOS);
Data.centerOfkSpaceMask = 1;
% %
sensOpt.null = 0;
aslObj = AslReconClass(sensOpt,Data,hdr);
aslObj.getCoilSensitivityMap() % estimated from the mean of all Reps, you may need to use the individual coils for each Rep measurement.

aslObj.setSuperResolutionImgRef(io.t1mprage.nii_flname,io.timeseries.nii_flname) % the first dataset (m0) is taken as refernce motion state
% set motion states
aslObj.SR.doMotionCorrection = 1;
aslObj.setMotionAffineTransforms([io.timeseries.nii_path 'transformix.mat'])
mkdir(io.results_path)
save([io.results_path 'aslObj_R1.mat'],'aslObj', 'io','-v7.3')
save([io.results_path 'masks.mat'],'brainExtractMask')




