
% set parameters for CBF quantification
opt.TR = 4000;
opt.T1t = 1300;
opt.lambda=0.9;
opt.alpha=0.85;
opt.LabDur=1800;
opt.PLD=1800;
%% root sum of squres of image channels
aslObj.SR.interp = 0;
[~,imgRSOS]=aslObj.reconCoilsImages();
imgRSOS_nifit = abs(aslObj.mapNativeSpaceToNiftiRefSpace(imgRSOS));
 
 imgRSOS_nifit_mc = 0*imgRSOS_nifit;

for i = 1:aslObj.nReps
    imgRSOS_nifit_mc(:,:,:,i) = mapSpaceAToSpaceBspm(abs(imgRSOS_nifit(:,:,:,i)),P(i),P(1),1); % motion correction with estimated motion vectors
end

brainExtractMask_asl = mapSpaceAToSpaceBspm(brainExtractMask,MrInfo,tsInfo)>0;
brainExtractMask_asl = imdilate(brainExtractMask_asl,rStrel())>0;
[x,y,z] = meshgrid(-4:1:4); r = (x.^2+y.^2+z.^2)<10;
brainExtractMask_asl_extended = imdilate(brainExtractMask_asl,r)>0;

M0 = imgRSOS_nifit_mc(:,:,:,1);
tmp = gauss3DFilter(M0,tsInfo.pixdim,20);
M0_extended = M0.*brainExtractMask_asl + tmp.*~brainExtractMask_asl.*brainExtractMask_asl_extended;  
M0_extended_smoothed = gauss3DFilter(M0_extended,tsInfo.pixdim,5);

M0_extended_smoothed_scaled = M0_extended_smoothed./(1-exp(-opt.TR/opt.T1t));

temp = mean(imgRSOS_nifit_mc(:,:,:,4:2:aslObj.nReps)-imgRSOS_nifit_mc(:,:,:,3:2:aslObj.nReps-1),4);
temp = max(0,temp);
temp = temp./abs(M0_extended_smoothed_scaled) .* brainExtractMask_asl;
temp(isnan(temp)) = 0;
temp(isinf(temp)) = 0;

cbf_imgRSOS_mc = CBF(temp,opt);
cbf_imgRSOS_mc_mrs = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc, tsInfo,MrInfo,0); 
% cbf_imgRSOS_mc_mrs_tl = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc, tsInfo,MrInfo,1); 
% %
save_nifti('cbf_imgRSOS_mc',cbf_imgRSOS_mc,[io.results_path 'nii_files' os_bar],tsNifit,tsInfo)
save_nifti('cbf_imgRSOS_mc_mrs',cbf_imgRSOS_mc_mrs,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)
% save_nifti('cbf_imgRSOS_mc_mrs_tl',cbf_imgRSOS_mc_mrs_tl,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)
%% SENSE reconstruction and motion correction
aslObj.SR.is = 0;
opt.RepKspaceData = inf;
opt.niter = 10;
imgSENSE_R1 = aslObj.SENSE_CG(opt);
imgSENSE_R1_nifti = abs(aslObj.mapNativeSpaceToNiftiRefSpace(imgSENSE_R1));
imgSENSE_R1_nifti_mc = 0*imgSENSE_R1_nifti;
%%
load([io.timeseries.nii_path 'transformix.mat'])
for i = 1:aslObj.nReps
    imgSENSE_R1_nifti_mc(:,:,:,i) = mapSpaceAToSpaceBspm(abs(imgSENSE_R1_nifti(:,:,:,i)),P(i),P(1),1); % motion correction with estimated motion vectors
end

M0 = imgSENSE_R1_nifti_mc(:,:,:,1);
tmp = gauss3DFilter(M0,tsInfo.pixdim,20);
M0_extended = M0.*brainExtractMask_asl + tmp.*~brainExtractMask_asl.*brainExtractMask_asl_extended;  
M0_extended_smoothed = gauss3DFilter(M0_extended,tsInfo.pixdim,5);

M0_extended_smoothed_scaled = M0_extended_smoothed./(1-exp(-opt.TR/opt.T1t));

temp = mean(imgSENSE_R1_nifti_mc(:,:,:,4:2:aslObj.nReps)-imgSENSE_R1_nifti_mc(:,:,:,3:2:aslObj.nReps-1),4);
temp = max(0,temp);
temp = temp./abs(M0_extended_smoothed_scaled) .* brainExtractMask_asl;
temp(isnan(temp)) = 0;
temp(isinf(temp)) = 0;

cbf_imgSENSE_mc = CBF(temp,opt);
cbf_imgSENSE_mc_mrs = mapSpaceAToSpaceBspm(cbf_imgSENSE_mc, tsInfo,MrInfo,0); 

save_nifti('cbf_imgSENSE_mc',cbf_imgSENSE_mc,[io.results_path 'nii_files' os_bar],tsNifit,tsInfo)
save_nifti('cbf_imgSENSE_mc_mrs',cbf_imgSENSE_mc_mrs,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)
%% 3D Linear regression PVC using FSL FAST PV estimates

pvgm = load_nii([io.rawdata_path 'fast_pv_estimates' os_bar 'pvgminasl.nii.gz']);
pvgm = flip(pvgm.img,1);
pvwm = load_nii([io.rawdata_path 'fast_pv_estimates' os_bar 'pvwminasl.nii.gz']);
pvwm = flip(pvwm.img,1);
%
opt.imCropFactor = [8,8,0];
opt.kernelSize = 5;
[lr_gm,lr_wm] = aslObj.LR3D_2(temp, pvgm,pvwm,opt);
gm_mask = pvgm >0.1;
wm_mask = pvwm>0.1;
cbf_imgRSOS_mc_LR_gm = CBF(lr_gm.*gm_mask,opt);
cbf_imgRSOS_mc_LR_gm_mrs = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc_LR_gm, tsInfo,MrInfo,0); 
% cbf_imgRSOS_mc_LR_gm_mrs_tl = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc_LR_gm, tsInfo,MrInfo,1); 
cbf_imgRSOS_mc_LR_wm = CBF(lr_wm.*wm_mask,opt);
cbf_imgRSOS_mc_LR_wm_mrs = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc_LR_wm, tsInfo,MrInfo,0); 
% cbf_imgRSOS_mc_LR_wm_mrs_tl = mapSpaceAToSpaceBspm(cbf_imgRSOS_mc_LR_wm, tsInfo,MrInfo,1); 

pvgm_mrs = mapSpaceAToSpaceBspm(pvgm,tsInfo,MrInfo,0);
pvwm_mrs = mapSpaceAToSpaceBspm(pvwm,tsInfo,MrInfo,0);

save_nifti('cbf_imgRSOS_mc_LR_gm',cbf_imgRSOS_mc_LR_gm,[io.results_path 'nii_files' os_bar],tsNifit,tsInfo)
save_nifti('cbf_imgRSOS_mc_LR_gm_mrs',cbf_imgRSOS_mc_LR_gm_mrs,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)

save_nifti('cbf_imgRSOS_mc_LR_wm',cbf_imgRSOS_mc_LR_wm,[io.results_path 'nii_files' os_bar],tsNifit,tsInfo)
save_nifti('cbf_imgRSOS_mc_LR_wm_mrs',cbf_imgRSOS_mc_LR_wm_mrs,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)
save_nifti('pvgm_mrs',uint16(pvgm_mrs*10^4),[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)
save_nifti('pvwm_mrs',uint16(pvwm_mrs*10^4),[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)

%% high-resolution guided reconstruction
% build a crop handel function to reduce computational load
if exist([io.results_path '\imCropHandelFunction.mat'],'file')
    load([io.results_path '\imCropHandelFunction.mat'])
else
    imCropHandelFunction = imcrop3d(abs(cbf_imgRSOS_mc_mrs));
    save([io.results_path '\imCropHandelFunction.mat'],'imCropHandelFunction')
end

%}
aslObj.SR.is = 1;
aslObj.SR.doMotionCorrection = 1;
optPrior.imCropHandle = imCropHandelFunction;
optPrior.sWindowSize = 3;
% build prior
aslObj.BuildSuperResolutionPrior(optPrior)
% calculate weights from Mr
optHrGr.MrSigma = 0.15;
% optHrGr.fSigma = 0.4;
MrImgPrior = gauss3DFilter(MrInfo.img,MrInfo.pixdim,2);
W = aslObj.getMrGaussianWeights(MrImgPrior,optHrGr.MrSigma);
%set PSF
PSF = [3,3,8];
aslObj.setPointSpreadFunction(PSF)
aslObj.SR.interp = 1;
% call SENSE_CG
optHrGr.niter = 100;
optHrGr.RepKspaceData = [3,4];%3:2:nReps-1
optHrGr.MrRegularizationParameter = 20;
optHrGr.MrPriorType = 'Quadratic'; % 'joint'
optHrGr.MrPreCompWeights = W;
optHrGr.display = 100;
optHrGr.stepSizeOptimization = 1;

optHrGr.psf = PSF;
optHrGr.report = 1; 
[Xc.img, Xc.M0, Xc.report] = aslObj.gradientDescent_4D_diff(optHrGr);
Xc.opt = optHrGr;
Xc.opt = rmfield(Xc.opt,'MrPreCompWeights');

imgHrGr_R1_4D_beta20_psf338 = Xc;
% %
if ~exist('trim_mask','var')
    trim_mask = sum(cbf_imgRSOS_mc_mrs)>0;
    save([io.results_path 'masks.mat'],'trim_mask','-append')
else
    load([io.results_path 'masks.mat'],'trim_mask')
end

% convert to cbf and save in nifti
brainExtractMask_extended = imdilate(brainExtractMask,r)>0;

tmp = gauss3DFilter(Xc.M0,tsInfo.pixdim,20);
M0_extended = Xc.M0.*brainExtractMask + tmp.*~brainExtractMask.*brainExtractMask_extended;  
M0_extended_smoothed = gauss3DFilter(M0_extended,tsInfo.pixdim,5);

M0_extended_smoothed_scaled = M0_extended_smoothed./(1-exp(-opt.TR/opt.T1t));

temp = Xc.img./abs(M0_extended_smoothed_scaled) .* brainExtractMask;
temp = max(0,temp);
temp(isnan(temp)) = 0;
temp(isinf(temp)) = 0;

cbf_imgHrGr_R1_4D_mc = CBF(temp,opt).*trim_mask;
save_nifti('cbf_imgHrGr_R1_4D_mc',cbf_imgHrGr_R1_4D_mc,[io.results_path 'nii_files' os_bar],MrNifti,MrInfo)

%% display all results
figure,
i = 100;
window = [0,80];
subplot(241),imshow(fliplr(rot90(MrInfo.img(:,:,i))),[]),title('T1-MPRAGE');
subplot(242),imshow(fliplr(rot90(abs(cbf_imgRSOS_mc_mrs(:,:,i)))),window),title('Standard');
subplot(243),imshow(fliplr(rot90(abs(cbf_imgRSOS_mc_LR_gm_mrs(:,:,i)))),window),title('3DLR');
subplot(244),imshow(fliplr(rot90(abs(cbf_imgHrGr_R1_4D_mc(:,:,i)))),window),title('Proposed');
i = 175;
subplot(245),imshow(squeeze(MrInfo.img(:,i,:)),[])
subplot(246),imshow(abs(squeeze(cbf_imgRSOS_mc_mrs(:,i,:))),window)
subplot(247),imshow(abs(squeeze(cbf_imgRSOS_mc_LR_gm_mrs(:,i,:))),window)
subplot(248),imshow(abs(squeeze(cbf_imgHrGr_R1_4D_mc(:,i,:))),window)

% %
% % analysis of CBF
gm = load_untouch_nii(io.t1mprage.gm_flname);
wm = load_untouch_nii(io.t1mprage.wm_flname);
gm = (gm.img./max(gm.img(:)))>0.7;
wm = (wm.img./max(wm.img(:)))>0.7;
Mean = @(x,mask) mean(min(abs(x(mask)),150));
gMean = [Mean(cbf_imgRSOS_mc_mrs,gm),Mean(cbf_imgRSOS_mc_LR_gm_mrs,gm),Mean(cbf_imgHrGr_R1_4D_mc,gm)];
wMean = [Mean(cbf_imgRSOS_mc_mrs,wm),Mean(cbf_imgRSOS_mc_LR_wm_mrs,wm),Mean(cbf_imgHrGr_R1_4D_mc,wm)];

figure,bar([gMean;wMean])
legend('Standard','3DLR','Proposed')
xticklabels({'Grey matter','White matter'})
setplot('','CBF (mL/100g/min')
grid on
colormap jet
