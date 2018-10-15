

function [Data,hdr,imgPcSos]=getwixASLdata(twixh,RemoveOS)

% only tested for data with 2 PE segments and 3 phasecor reference lines,
% need different data to generalize it


twixh.('image').flagRampSampRegrid = 1;
twixh.('phasecor').flagRampSampRegrid = 1;

if RemoveOS
    twixh.('image').flagRemoveOS = RemoveOS;
    twixh.('phasecor').flagRemoveOS = RemoveOS;
end
if ~isempty(twixh.('image').RawDataCorrectionFactors)
    twixh.('image').flagDoRawDataCorrect = 1;
end
if ~isempty(twixh.('phasecor').RawDataCorrectionFactors)
    twixh.('phasecor').flagDoRawDataCorrect = 1;
end

hdr.data.NCol = twixh.('image').NCol;
hdr.data.NLin = twixh.('image').NLin;
hdr.data.NCha = twixh.('image').NCha;
hdr.data.NPar = twixh.('image').NPar;
hdr.data.NRep = twixh.('image').NRep;
hdr.data.NSeg = twixh.('image').NSeg;
hdr.data.NSli = twixh.('image').NSli;
hdr.data.NAve = twixh.('image').NAve;
hdr.data.NAcq = twixh.('image').NAcq;
hdr.data.is3d = twixh.('image').NPar>1;
hdr.data.RemoveOS = RemoveOS;

imageData = twixh.image{''};
phasecorData = twixh.phasecor{''};

if RemoveOS
    OS = 2;
else
    OS = 1;
end
data = zeros(hdr.data.NCol/OS, hdr.data.NLin, hdr.data.NPar, hdr.data.NCha, hdr.data.NRep);
% phase correction
for r = 1:hdr.data.NRep
    for c = 1:hdr.data.NCha
        
        imgDataS1 = squeeze(imageData(:,c,:,:,r,1));
        imgDataS2 = squeeze(imageData(:,c,:,:,r,2));
        
        ref = [phasecorData(:,c,1,r,1),...% seg 1, ave 1
            phasecorData(:,c,2,r,1),...% seg 1, ave 2
            phasecorData(:,c,1,r,2)]';% seg 2, ave 1
        
        fftseg1 = fft((ref(1,:)+ref(2,:))/2);
        fftseg2 = fft(ref(3,:));
        
        reflect = 0;
        if reflect %the 1st and 3rd are reverse scan
            mfilter = fftseg2./fftseg1;
        else
            mfilter = fftseg1./fftseg2;
        end
        mfilter = mfilter./abs(mfilter);
        
        imgDataS1C = ifft(fft(imgDataS1,[],1).*repmat(mfilter',[1,hdr.data.NLin,hdr.data.NPar]),[],1);
        
        data(:,:,:,c,r) = imgDataS1C + imgDataS2;
        
    end
end

hdr.img.matrixSize = [twixh.hdr.Meas.NImageCols,twixh.hdr.Meas.NImageLins,twixh.hdr.Meas.NImagePar];
hdr.img.FOV = [twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV,twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV,twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness];
hdr.img.voxelSize = [hdr.img.FOV./hdr.img.matrixSize];
sPosition = twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
position = [nan, nan, nan];
if isfield(sPosition,'dSag'), position(1) = sPosition.dSag; end
if isfield(sPosition,'dCor'), position(2) = sPosition.dCor; end
if isfield(sPosition,'dTra'), position(3) = sPosition.dTra; end
hdr.img.position = position;

hdr.data.ftSize = [twixh.hdr.Config.NColMeas/OS, twixh.hdr.Config.NPeftLen,twixh.hdr.Config.NPaftLen];
hdr.data.sqzSize = twixh.('image').sqzSize([1,3,4,2,5:ndims(data)]);
hdr.data.sqzDim = twixh.('image').sqzDims([1,3,4,2,5:ndims(data)]);
if isfield(twixh.hdr.MeasYaps.sKSpace,'dSliceOversamplingForDialog')
    hdr.data.sliceOverSampling = twixh.hdr.MeasYaps.sKSpace.dSliceOversamplingForDialog;
end
hdr.data.centerOfkSpace = [twixh.('image').centerCol(1),twixh.('image').centerLin(1),twixh.('image').centerPar(1)];

hdr.data.centerOfkSpace(1) = (hdr.data.centerOfkSpace(1)-1)/OS+1;

% centrelize k-sapce data and create a corresponding mask (used for
% itertaive reconstruction)
[data, centerOfkSpaceMask] = centrelizeKspace(data,hdr.data.sqzSize,hdr.data.ftSize,hdr.data.centerOfkSpace);
Data.centerOfkSpaceMask = centerOfkSpaceMask;
Data.ftKspaceData= data;

% single channel reconstruction
imgPC = zeros(size(Data.ftKspaceData));
for r = 1:hdr.data.NRep
    for c = 1:hdr.data.NCha
        imgPC(:,:,:,c,r) = fftshift(ifftn(ifftshift( Data.ftKspaceData(:,:,:,c,r))));
    end
end
imgPcSos = sqrt(squeeze(sum(abs(imgPC).^2,4)));

% imgPcSos = cropASLimg(imgPcSos,hdr);

end 


function [kSpaceDataCentre, centerOfkSpaceMask] = centrelizeKspace(kSpaceData,sqzSzie,FtMatrixSize,centerOfkSpace)

% sqzSzie: size of k-space data returned by twix
% FtMatrixSize: size of k-space for FT trasnfrom
% centerOfkSpace: center of k-space data retured by twix

centerOfkSpaceMask = true(FtMatrixSize(1:3));
delta = floor(FtMatrixSize(1:3)/2+1- centerOfkSpace) ;
if all(delta ==0)
    kSpaceDataCentre = kSpaceData;
    return;
end

kSpaceDataCentre = zeros([FtMatrixSize,sqzSzie(4:end)],'single');

if delta(1)>0
    centerOfkSpaceMask(1:delta(1),:,:) = 0;
elseif delta(1)<0
    centerOfkSpaceMask(FtMatrixSize(1)+delta(1)+1:end,:,:) = 0;
end
if delta(2)>0
    centerOfkSpaceMask(:,1:delta(2),:) = 0;
elseif delta(2)<0
    centerOfkSpaceMask(:,FtMatrixSize(2)+delta(2)+1:end,:) = 0;
end
if delta(3)>0
    centerOfkSpaceMask(:,:,1:delta(3)) = 0;
elseif delta(3)<0
    centerOfkSpaceMask(:,:,FtMatrixSize(3)+delta(3)+1:end) = 0;
end
% %

kSpaceDataCentre(1:sqzSzie(1),1:sqzSzie(2),1:sqzSzie(3),:,:,:,:) = kSpaceData;

for i = 1:length(delta)
    if delta(i)~=0
        kSpaceDataCentre = circshift(kSpaceDataCentre,delta(i),i);
    end
end

end


% function imagOut = cropASLimg(imagTimeSeries, hdr)
% 
% 
% 
% if hdr.data.RemoveOS
%     imagOut = imagTimeSeries;
% else
%     idx = hdr.data.NCol/4;
%     imagOut = imagTimeSeries(idx+1:end-idx,:,:,:);
% end
% imagOut = flip(imagOut,3);
% imagOut = imagOut(:,:,1+(1:hdr.img.matrixSize(3)),:);
% temp = zeros([hdr.img.matrixSize,size(imagOut,4)]);
% for i = 1:size(imagOut,4)
%     temp(:,:,:,i) = resize(imagOut(:,:,:,i),hdr.img.matrixSize);
% end
% imagOut = temp;
% 
% end




