
classdef AslReconClass < matlab.mixin.Copyable
    
    
    properties (SetAccess = public)
        
        isSimulation
        is3D
        hdr
        kScaleFactor
        kUnderSampling
        CentralizingMask
        PhaseFoV
        ReadFoV
        dThickness
        Coil
        Prior
        nReps
        SR
        ftKspaceData
        ftKspaceDataSize
        ftImgSize
        reconImgSize
        reconVoxelSize
        
    end
    
    methods
        function ObjMRI = AslReconClass(varargin)
            % call the constructor
            ObjMRI.isSimulation     = 0;
            ObjMRI.is3D             = [];
            ObjMRI.kScaleFactor     = 1;
            ObjMRI.CentralizingMask             = 1;
            ObjMRI.kUnderSampling.is            = 0;
            ObjMRI.kUnderSampling.trajectory    = 1;
            ObjMRI.kUnderSampling.factor        =0;
            ObjMRI.kUnderSampling.method        = 'radial';
            ObjMRI.Coil.N                       = 1;
            ObjMRI.Coil.sensitivityMap           = [];
            ObjMRI.Coil.estimationMethod         = 'rsos';
            ObjMRI.Coil.supportMask              = [];
            
            ObjMRI.SR.is = 0;
            ObjMRI.SR.imgRefH =[];
            ObjMRI.SR.imgRefL =[];
            ObjMRI.SR.method = 'matlab-spm'; %matlab
            ObjMRI.SR.interp = 1;
            ObjMRI.SR.PSF = [0 0 0]; %mm anisotric Gaussian
            ObjMRI.SR.affineTransforms = [];
            ObjMRI.SR.doMotionCorrection = 0;
            
            
            if ~isempty(varargin) && isstruct(varargin{1})
                % update object's properties from user's options
                ObjMRI = getFiledsFromUsersOpt(ObjMRI,varargin{1});
            else
                error('input is required')
            end
            if ObjMRI.isSimulation
                
            else
                data  = varargin{2};
                hdr  = varargin{3};
                getwixMRdata(ObjMRI,data,hdr);
            end
        end
    end
    
    
    
    methods (Access = public)
        % ///////////////////////////////////////////////////////////////////////////////
        %                           DATA PROCESSING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        function ObjMRI = getwixMRdata(ObjMRI,data,hdr)
            
            ObjMRI.hdr = hdr;
            ObjMRI.Coil.N = hdr.data.NCha;
            ObjMRI.is3D    = hdr.data.is3d;
            ObjMRI.nReps = hdr.data.NRep;
            ObjMRI.ftKspaceData  = data.ftKspaceData;%ObjMRI.kScale(data.ftKspaceData);
            ObjMRI.ftKspaceDataSize = hdr.data.ftSize;
            ObjMRI.ftImgSize = ObjMRI.ftKspaceDataSize(1:3);
            ObjMRI.reconImgSize = hdr.img.matrixSize;
            ObjMRI.reconVoxelSize  = hdr.img.voxelSize;
            ObjMRI.CentralizingMask = 1;%logical(data.centerOfkSpaceMask);
            ObjMRI.kUnderSampling.trajectory = data.centerOfkSpaceMask;
        end
        
        function scaledData = kScale(ObjMRI,data)
            ObjMRI.kScaleFactor  = ObjMRI.kScaleFactor ./(max(abs(data(:))));
            scaledData = ObjMRI.kScaleFactor * data;
        end
        
        function doRestropetiveKspaceUndersampling(ObjMRI,arg)
            
            arg.null = 0;
            opt.undersamplingMethod = 'radial'; % 'cartesian', 'spiral', 'vd-random'
            opt.radial_nspokes = 40;
            opt.cartesian_R_pe = 2; % phase encoding
            opt.cartesian_R_se = 1; % slice encoding
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if strcmpi(opt.undersamplingMethod, 'cartesian')%
                ObjMRI.kUnderSampling.cartesian_R_pe = opt.cartesian_R_pe;
                ObjMRI.kUnderSampling.cartesian_R_se = opt.cartesian_R_se;
                Xi = false([ObjMRI.ftKspaceDataSize]);
                
                Xi(:,1:opt.cartesian_R_pe:end,1:opt.cartesian_R_se:end) = 1;
                
            elseif strcmpi(opt.undersamplingMethod, 'radial')% stack-of-radials
                ObjMRI.kUnderSampling.radial_nspokes = opt.radial_nspokes;
                % stack of 2D radial trajectories
                Col = ObjMRI.ftKspaceDataSize(1);
                Lin = ObjMRI.ftKspaceDataSize(2);
                n = max(Col,Lin);
                Theta = linspace(0,pi,opt.radial_nspokes+1);Theta(end) = [];
                Xi = zeros(n,n);
                for theta = Theta
                    t = linspace(-1,1,3*n)*n;
                    X = round(t.*cos(theta)) + n/2+1;
                    Y = round(t.*sin(theta)) + n/2+1;
                    I = find(X>0 & X<=n & Y>0 & Y<=n);
                    X = X(I);
                    Y = Y(I);
                    Xi(X+(Y-1)*n) = 1;
                end
                dx = Col - Lin;
                if dx>0
                    subset = (abs(dx)/2+1):(Col-abs(dx)/2);
                    Xi = Xi(:,subset);
                elseif dx<0
                    subset = (abs(dx)/2+1):(Lin-abs(dx)/2);
                    Xi = Xi(subset,:);
                end
                Xi = repmat(Xi,[1,1,ObjMRI.ftKspaceDataSize(3)]);
            elseif strcmpi(opt.undersamplingMethod, 'spiral')%
                error('not impemented yet')
            elseif strcmpi(opt.undersamplingMethod, 'vd-random')%
                error('not impemented yet')
            else
                error('unknown')
            end
            ObjMRI.kUnderSampling.is = 1;
            ObjMRI.kUnderSampling.method = opt.undersamplingMethod;
            ObjMRI.kUnderSampling.trajectory = Xi;
            
            for rep = 1:ObjMRI.nReps
                for coil = 1:ObjMRI.Coil.N
                    ObjMRI.ftKspaceData(:,:,:,coil,rep) = ObjMRI.ftKspaceData(:,:,:,coil,rep).*ObjMRI.kUnderSampling.trajectory;
                end
            end
            ObjMRI.kUnderSampling.factor = prod(ObjMRI.ftKspaceDataSize)/numel(find(ObjMRI.kUnderSampling.trajectory));
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           FOURIER ENCODING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function m = F0(ObjMRI,x,mask)
            if ObjMRI.isSimulation
                m = fftn(x).*mask;
            else
                m = ifftshift(fftn(fftshift(x))).*mask;%
            end
        end
        
        function x = FH0(ObjMRI,m,mask)
            if ObjMRI.isSimulation
                x = ifftn(m.*mask);
            else
                x = fftshift(ifftn(ifftshift(m.*mask)));
            end
        end
        
        function m = F(ObjMRI,x,jm)
            % perfroms MRI forward operator
            if nargin==2
                jm = 1; % by default 1th motion state is considered as reference state to which all other motion states are mapped
            end
            if ObjMRI.kUnderSampling.is
                mask = ObjMRI.kUnderSampling.trajectory;
            else
                mask = 1;%ObjMRI.CentralizingMask;
            end
            if ObjMRI.SR.is
                x = ObjMRI.DownSampling(x,ObjMRI.SR.interp,jm);
            end
            m = zeros([ObjMRI.ftKspaceDataSize, ObjMRI.Coil.N]);
            
            for i = 1:ObjMRI.Coil.N
                if isempty(ObjMRI.Coil.sensitivityMap)
                    Sen = 1;
                else
                    if ObjMRI.is3D
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,:,i);
                    else
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,i);
                    end
                end
                K = ObjMRI.F0(x.*Sen,mask);
                if ObjMRI.is3D
                    m(:,:,:,i) = K;
                else
                    m(:,:,i) = K;
                end
                
            end
        end
        
        function x = FH(ObjMRI,m,jm)
            % perfroms MRI adjoint operator
            if nargin==2
                jm = 1; % by default 1th motion state is considered as reference state to which all other motion states are mapped
            end
            if ObjMRI.kUnderSampling.is
                mask = ObjMRI.kUnderSampling.trajectory;
            else
                mask = 1;%ObjMRI.CentralizingMask;
            end
            x = zeros(ObjMRI.ftKspaceDataSize);
            
            for i = 1: ObjMRI.Coil.N
                if isempty(ObjMRI.Coil.sensitivityMap)
                    Sen = 1;
                else
                    if ObjMRI.is3D
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,:,i);
                        K = m(:,:,:,i);
                    else
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,i);
                        K = m(:,:,i);
                    end
                end
                x = x + conj(Sen).*ObjMRI.FH0(K,mask);
            end
            if ObjMRI.SR.is
                x = ObjMRI.UpSampling(x,ObjMRI.SR.interp,jm);
            end
        end
        
        function x = FHF(ObjMRI,m,jm)
            if nargin==2
                jm=1;
            end
            % perfroms MRI Gram Matrix opertor
            x = ObjMRI.FH(ObjMRI.F(m,jm),jm);
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           RECONSTRUCTION FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function [ImgCoils, ImgCoilsSOS] = reconCoilsImages(ObjMRI,RepKspaceData)
            if nargin==1
                %for 3D data, if there are 5th dimension,
                %RepKspaceData=1,2,3,.... else all
                RepKspaceData = 1:ObjMRI.nReps;
                
            end
            % ImgCoils from fully-sampled k-space, so use CentralizingMask
            ImgCoils = zeros([ObjMRI.ftKspaceDataSize, ObjMRI.Coil.N,length(RepKspaceData)]);
            for ii = 1:length(RepKspaceData)
                for i = 1:ObjMRI.Coil.N
                    if ObjMRI.is3D
                        ImgCoils(:,:,:,i,ii) = ObjMRI.FH0(ObjMRI.ftKspaceData(:,:,:,i,RepKspaceData(ii)),ObjMRI.CentralizingMask);
                    else
                        ImgCoils(:,:,i) = ObjMRI.FH0(ObjMRI.ftKspaceData(:,:,i),ObjMRI.CentralizingMask);
                    end
                end
            end
            if ObjMRI.is3D
                ImgCoilsSOS = zeros([ObjMRI.ftKspaceDataSize, length(RepKspaceData)]);
                for ii = 1:length(RepKspaceData)
                    ImgCoilsSOS(:,:,:,ii) = ObjMRI.RSS(ImgCoils(:,:,:,:,RepKspaceData(ii)));
                end
                if isempty(ObjMRI.Coil.supportMask)
                    ObjMRI.Coil.supportMask = ObjMRI.MakeCoilSupportMask(ImgCoilsSOS(:,:,:,1));
                end
            else
                ImgCoilsSOS = ObjMRI.RSS(ImgCoils);
                if isempty(ObjMRI.Coil.supportMask)
                    ObjMRI.Coil.supportMask = ObjMRI.MakeCoilSupportMask(ImgCoilsSOS);
                end
            end
            
        end
        
        function setPointSpreadFunction(ObjMRI, psf)
            if nargin==1
                psf = 0;
            end
            ObjMRI.SR.PSF = psf;
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           HELPER FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function Y = RSS(ObjMRI,X)
            % Calculates root of sum of squares
            if ndims(X)==2 %#ok<ISMAT> % Gradient vectors
                dims = 2;
            else % images
                if ObjMRI.is3D
                    dims = 4;
                else
                    dims = 3;
                end
            end
            Y  = sqrt(sum(abs(X).^2,dims));
        end
        
        function x = Magnitude(~,x)
            x = sqrt(sum(abs(x(:)).^2));
        end
        
        function display(ObjMRI)
            disp(ObjMRI)
            methods(ObjMRI)
        end
        
        function ObjMRI = Revise(ObjMRI,opt)
            % to revise the properties of a given object without
            % re-instantiation
            vfields = fieldnames(opt);
            prop = properties(ObjMRI);
            for i = 1:length(vfields)
                field = vfields{i};
                if sum(strcmpi(prop, field )) > 0
                    ObjMRI.(field) = opt.(field);
                end
            end
        end
        % ///////////////////////////////////////////////////////////////////////////////
        %                           COIL FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        function setCoilSensitivityMap(ObjMRI,CSM)
            % in the case of a pre-calculated coil map
            ObjMRI.Coil.sensitivityMap = CSM;
        end
        
        function getCoilSensitivityMap(ObjMRI,arg)
            arg.null = 0;
            opt.coilEstimationMethod = 'rsos';
            %opt.useCenterOfKspacePercentage = 100; to be done
            opt = getFiledsFromUsersOpt(opt,arg);
            
            [ImgCoils, ImgCoilsSOS] = reconCoilsImages(ObjMRI);
            
            if strcmpi(opt.coilEstimationMethod,'rsos')
                if ndims(ImgCoils)==5
                    CoilEst = 0*ImgCoils;
                    for i = 1:size(ImgCoils,5)
                        CoilEst(:,:,:,:,i) =  ImgCoils(:,:,:,:,i)./repmat(ImgCoilsSOS(:,:,:,i),[1,1,1,ObjMRI.Coil.N]);
                    end
                    CoilEst = mean(CoilEst,5);
                else
                    CoilEst = ImgCoils./repmat(ImgCoilsSOS,[1,1,1,ObjMRI.Coil.N]);
                end
                CoilEstFilt = 0*CoilEst;
                for i = 1:ObjMRI.Coil.N
                    CoilEstFilt(:,:,:,i) = gauss3DFilter(CoilEst(:,:,:,i),[1,1,1],2);
                end
                CoilEst = CoilEstFilt;
                % mask the CSMs
                for i = 1:ObjMRI.Coil.N
                    CoilEst(:,:,:,i) = CoilEst(:,:,:,i).*ObjMRI.Coil.supportMask;
                end
                ObjMRI.Coil.sensitivityMap = (CoilEst);
            else
                
            end
        end
        
        function mask = MakeCoilSupportMask(ObjMRI,imgCoilSOS)
            
            if ObjMRI.is3D
                mask = imgaussfilt3(imgCoilSOS./max(imgCoilSOS(:)),8);
            else
                mask = imgaussfilt(imgCoilSOS./max(imgCoilSOS(:)),8);
            end
            
            level = graythresh(mask);
            if level==0
                fprintf('MakeCoilSupportMask:: zero threshold level !!\n')
                mask = 1;
                return
            end
            mask = mask >(level*0.25);
            [x,y,z] = meshgrid(-1:1:1);
            r = (x.^2+y.^2+z.^2)<2;
            
            for i = 1:5
                mask = imdilate(single(mask),r);
            end
            if ObjMRI.is3D
                for i=1:size(mask,3)
                    mask(:,:,i) = imfill(mask(:,:,i),'holes');
                end
            else
                mask= imfill(mask,'holes');
            end
            mask = logical(mask);
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           PRIOR OBJECT FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function setSuperResolutionImgRef(ObjMRI,imgRefH,imgRefL,opt)
            
            % 1) imgRefH and imgRefH: nifti (.nii) file address of the high
            % resolution image (e.i T1-weighted MR) and the low-resolution
            % image (i.e. one of ASL time serie image), use SPM's nifti
            % function to get info, sampling method is set to 'matlab-spm'
            
            % 2) imgRefH and imgRefH: spatial structre generated by
            % getNifiDataInfo (which calls SPM's nifti function, sampling
            % method is set to 'matlab-spm'
            
            % 3) imgRefH/imgRefH: matlab's imref3d (Reference 3-D image to
            % world coordinate object), sampling method is set to 'matlab-imref3d'
            
            opt.null = 0;
            if ischar(imgRefH) % --> imgRefH and imgRefL are address files
                imgRefH = ObjMRI.getNiftiDataInfo(imgRefH);
                imgRefL = ObjMRI.getNiftiDataInfo(imgRefL);
                ObjMRI.SR.method = 'matlab-spm';
                if isfield(opt,'interp') %0: nearest neighbours, 1: Trilinear
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 1;
                end
            elseif isstruct(imgRefH)
                ObjMRI.SR.method = 'matlab-spm';
                if isfield(opt,'interp') %0: nearest neighbours, 1: Trilinear
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 1;
                end
            elseif strcmpi(class(imgRefH),'imref3d')
                ObjMRI.SR.method = 'matlab-imref3d';
                if isfield(opt,'interp')
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 'linear';
                end
            end
            ObjMRI.SR.imgRefH = imgRefH;
            ObjMRI.SR.imgRefL = imgRefL;
        end
        
        function BuildSuperResolutionPrior(ObjMRI,arg)
            arg.null = 0;
            opt.imCropFactor = [8,8,8];
            opt.sWindowSize = 3;
            opt.lWindowSize = 1;
            opt.imCropHandle = [];
            
            opt = getFiledsFromUsersOpt(opt,arg);
            
            ObjMRI.SR.is = 1;
            if isempty(ObjMRI.SR.imgRefH)
                error('First set setSuperResolutionImgRef');
            end
            if strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                opt.ImageSize = ObjMRI.SR.imgRefH.ImageSize;
                ObjMRI.reconVoxelSize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY, ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
            elseif strcmpi(ObjMRI.SR.method,'matlab-spm')
                opt.ImageSize = ObjMRI.SR.imgRefH.dim;
                ObjMRI.reconVoxelSize = ObjMRI.SR.imgRefH.pixdim;
            end
            if ~isempty(opt.imCropHandle)
                opt.imCropFactor = 0; % for Siemens mMR
            end
            ObjMRI.Prior = PriorsClass(opt);
        end
        
        function BuildNativeResolutionPrior(ObjMRI,arg)
            ObjMRI.SR.is = 0;
            arg.null = 0;
            opt.imCropFactor = [8,8,0];
            opt.sWindowSize = 3;
            opt.lWindowSize = 1;
            opt.imCropHandle = [];
            
            opt = getFiledsFromUsersOpt(opt,arg);
            
            opt.ImageSize = ObjMRI.reconImgSize;
            
            if ~isempty(opt.imCropHandle)
                opt.imCropFactor = 0; % for Siemens mMR
            end
            ObjMRI.Prior = PriorsClass(opt);
        end
        
        function W = getMrGaussianWeights(ObjMRI,MrImage,MrSigma)
            % normalize MrImage to [0,1]
            MrImage = abs(MrImage)./max(abs(MrImage(:)));
            
            W = ObjMRI.Prior.W_GaussianKernel(MrImage,MrSigma);
            W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
        end
        % ///////////////////////////////////////////////////////////////////////////////
        %                           SPACE MAPPING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        function setMotionAffineTransforms(ObjMRI,matfile)
            s = load(matfile);
            
            ObjMRI.SR.affineTransforms = s.P(3:2:end); % skip controls, assuming that there's motion
        end
        
        function imagOut = mapNativeSpaceToNiftiRefSpace(ObjMRI,nativeSpaceImgs,interp)
            % this function maps the reconstcuted images from the kspace into the
            % correposnding nitifi space of the the same images reconstructed by the
            % scanner, need validation for other dataset, currently is hard coded,
            % better to figure out the coordinate system of the images and do
            % resampling
            if nargin ==2
                interp = 'linear';
            end
            if ObjMRI.hdr.data.RemoveOS
                imagOut = nativeSpaceImgs;
            else
                idx = ObjMRI.hdr.data.NCol/4;
                imagOut = nativeSpaceImgs(idx+1:end-idx,:,:,:);
            end
            imagOut = flip(imagOut,3);
            imagOut = imagOut(:,:,1+(1:ObjMRI.hdr.img.matrixSize(3)),:);
            nRe = size(imagOut,4);
            temp = zeros([ObjMRI.hdr.img.matrixSize,nRe]);
            for i = 1:nRe
                temp(:,:,:,i) = resize(imagOut(:,:,:,i),ObjMRI.hdr.img.matrixSize,interp);
            end
            
            temp = circshift(temp,-1,1);
            temp = circshift(temp,-1,2);
            
            imagOut = temp;
        end
        
        function imagOut = mapNiftiRefSpaceToNativeSpace(ObjMRI,refSpaceImgs,interp)
            if nargin ==2
                interp = 'linear';
            end
            % inverts the mapNativeSpaceToNiftiRefSpace
            temp = refSpaceImgs;
            temp = circshift(temp,+1,2);
            temp = circshift(temp,+1,1);
            
            ftSize0 = [ObjMRI.hdr.data.ftSize(1),ObjMRI.hdr.img.matrixSize(2:3)];
            nRe = size(temp,4);
            imagOut = zeros([ObjMRI.hdr.data.ftSize,nRe]);
            ii = 1+(1:ObjMRI.hdr.img.matrixSize(3));
            for i = 1:nRe
                imagOut(:,:,ii,i) = resize(temp(:,:,:,i),ftSize0,interp);
            end
            imagOut = flip(imagOut,3);
            if ~ObjMRI.hdr.data.RemoveOS
                idx = ObjMRI.hdr.data.NCol/4;
                temp = imagOut;
                imagOut = zeros([ObjMRI.hdr.data.ftSize, nRe]);
                imagOut (idx+1:end-idx,:,:,:) = temp;
            end
        end
        
        function imgOut = mapNativeSpaceToMrSpace(ObjMRI, nativeSpaceImgs)
            
            temp = ObjMRI.mapNativeSpaceToNiftiRefSpace(nativeSpaceImgs);
            
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                imgOut = mapSpaceAToSpaceBspm(temp,ObjMRI.SR.imgRefL,ObjMRI.SR.imgRefH,ObjMRI.SR.interp);
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                error('todo');
            end
        end
        
        function imgOut = mapMrSpaceToNativeSpace(ObjMRI, mrImg)
            % inverts mapNativeSpaceToMrSpace
            
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                temp = mapSpaceAToSpaceBspm(mrImg,ObjMRI.SR.imgRefH,ObjMRI.SR.imgRefL,ObjMRI.SR.interp);
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                error('todo');
            end
            imgOut = ObjMRI.mapNiftiRefSpaceToNativeSpace(temp);
            
        end
        
        function imgOut = mrSpaceToNativeSpaceJthMotionState(ObjMRI, Img,inverse,jm)
            if nargin==2, inverse = 0;  end
            if nargin<=3, jm = 1;  end
            
            if inverse
                Img = ObjMRI.mapNativeSpaceToNiftiRefSpace(Img);
                %                 imgOut =zeros(ObjMRI.SR.imgRefH.dim(1:3));
                %             else
                %                 imgOut =zeros(ObjMRI.SR.imgRefL.dim(1:3));
            end
            P = ObjMRI.SR.affineTransforms;
            Pmr = ObjMRI.SR.imgRefH;
            
            
            if inverse
                imgOut = mapSpaceAToSpaceBspm(Img,P(jm),Pmr,ObjMRI.SR.interp);
            else % map mr space to nifit reference space jth motion state
                imgOut = mapSpaceAToSpaceBspm(Img,Pmr,P(jm),ObjMRI.SR.interp);
            end
            
            % map nifit reference space to native space
            if ~inverse
                imgOut = ObjMRI.mapNiftiRefSpaceToNativeSpace(imgOut);
            end
        end
        
        function imagOut = T1_mapNiftiRefSpaceToNativeSpace(ObjMRI,RefSpaceImg)
            % this function maps the reconstcuted images from the kspace into the
            % correposnding nitifi space of the the same images reconstructed by the
            % scanner, need validation for other dataset, currently is hard coded,
            % better to figure out the coordinate system of the images and do
            % resampling
            % if nargin ==2
            %     interp = 'linear';
            % end
            imagOut = circshift(RefSpaceImg,+1,2);
            imagOut = circshift(imagOut,-1,1);
            imagOut = flip(flip(imagOut,1),3);
            imagOut = permute(imagOut,[3,2,1]);
            
            temp = imagOut;
            
            imagOut = zeros([size(temp,1),size(temp,2),ObjMRI.ftKspaceDataSize(3)]);
            nSlice = ObjMRI.hdr.img.matrixSize(3);
            idx = (ObjMRI.ftKspaceDataSize(3) - nSlice )/2;
            
            imagOut(:,:,idx:nSlice+idx-1) = temp;
            
            if ObjMRI.hdr.data.RemoveOS==0
                idx = ObjMRI.hdr.data.NCol/4;
                temp = imagOut;
                imagOut = zeros([2*size(temp,1),size(temp,2),size(temp,3)]);
                imagOut(idx+1:end-idx,:,:,:) = temp;
                
            end
            
        end
        
        function imagOut = T1_mapNativeSpaceToNiftiRefSpace(ObjMRI,nativeSpaceImg)
            % this function maps the reconstcuted images from the kspace into the
            % correposnding nitifi space of the the same images reconstructed by the
            % scanner, need validation for other dataset, currently is hard coded,
            % better to figure out the coordinate system of the images and do
            % resampling
            % if nargin ==2
            %     interp = 'linear';
            % end
            if ObjMRI.hdr.data.RemoveOS
                imagOut = nativeSpaceImg;
            else
                idx = ObjMRI.hdr.data.NCol/4;
                imagOut = nativeSpaceImg(idx+1:end-idx,:,:,:);
            end
            
            % remove slice oversampling
            nSlice = ObjMRI.hdr.img.matrixSize(3);
            idx = (size(imagOut,3) - nSlice )/2;
            imagOut = imagOut(:,:,idx:nSlice+idx-1);
            
            imagOut = permute(imagOut,[3,2,1]);
            imagOut = flip(flip(imagOut,3),1);
            
            imagOut = circshift(imagOut,1,1);
            imagOut = circshift(imagOut,-1,2);
        end
        
        function ImgNew = DownSampling(ObjMRI,Img,interp,jm)
            
            
            if nargin==2, interp = ObjMRI.SR.interp; end
            if nargin<=3, jm = 1; end
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                
                if ObjMRI.SR.is == 0 % in case of low-res deconv
                    ImgNew = Img;
                else
                    if ObjMRI.SR.doMotionCorrection
                        ImgNew = ObjMRI.mrSpaceToNativeSpaceJthMotionState(Img,0,jm);
                    else
                        ImgNew = ObjMRI.mapMrSpaceToNativeSpace(Img);
                    end
                end
                
                if any(ObjMRI.SR.PSF)
                    voxelsize = ObjMRI.SR.imgRefL.pixdim;
                    ImgNew = gauss3DFilter(ImgNew,voxelsize,ObjMRI.SR.PSF);
                end
                
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                if any(ObjMRI.SR.PSF)
                    voxelsize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY,ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
                    Img = gauss3DFilter(Img,voxelsize,ObjMRI.SR.PSF);
                end
                ImgNew = ImageResample(Img, ObjMRI.SR.imgRefH, ObjMRI.SR.imgRefL,interp);
                
            end
            
        end
        
        function ImgNew = UpSampling(ObjMRI,Img,interp,jm)
            
            if nargin==2, interp = ObjMRI.SR.interp;  end
            if nargin<=3, jm = 1; end
            if ObjMRI.SR.is
                % trim, --< during Usampling and Downsampling there are some
                % very high-intesne planes, need to be removed.
                a = 0*Img;
                pitch = 2;
                if ObjMRI.is3D
                    a(pitch:end-pitch,pitch:end-pitch,pitch:end-pitch) = 1;
                else
                    a(pitch:end-pitch,pitch:end-pitch) = 1;
                end
                Img = a.*Img;
            end
            %
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                
                if any(ObjMRI.SR.PSF)
                    voxelsize = ObjMRI.SR.imgRefL.pixdim;
                    Img = gauss3DFilter(Img,voxelsize,ObjMRI.SR.PSF);
                end
                
                if ObjMRI.SR.is == 0 % in case of low-res deconv
                    ImgNew = Img;
                else
                    if ObjMRI.SR.doMotionCorrection
                        ImgNew = ObjMRI.mrSpaceToNativeSpaceJthMotionState(Img,1,jm);
                    else
                        ImgNew = ObjMRI.mapNativeSpaceToMrSpace(Img);
                    end
                end
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                ImgNew = ImageResample(Img, ObjMRI.SR.imgRefL, ObjMRI.SR.imgRefH,interp);
                if any(ObjMRI.SR.PSF)
                    voxelsize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY,ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
                    ImgNew = gauss3DFilter(ImgNew,voxelsize,ObjMRI.SR.PSF);
                end
            end
            
            if ObjMRI.SR.is
                a = 0*ImgNew;
                pitch = 3;
                if ObjMRI.is3D
                    a(pitch:end-pitch,pitch:end-pitch,pitch:end-pitch) = 1;
                else
                    a(pitch:end-pitch,pitch:end-pitch) = 1;
                end
                ImgNew = a.*ImgNew;
            end
            
        end
        
        function [out, tmp]= getNiftiDataInfo(ObjMRI,flname)
            % This function requires both NIFIT and SPM libraries
            nifiObj =nifti(flname);
            out.mat = nifiObj.mat;
            
            tmp = load_untouch_nii(flname);
            %             out.img = double(tmp.img);
            out.dim = tmp.hdr.dime.dim(2:tmp.hdr.dime.dim(1)+1);
            out.pixdim = tmp.hdr.dime.pixdim(2:4);
            out.flname = flname;
            % pulling out 4th dimesion,in case of 4D nifti datafile
            out.dim = out.dim(1:3);
            % for compatibility
            out.ImageSize = out.dim;
            out.dataType = tmp.hdr.dime.datatype;
            
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           OPTIMIZATION ALGORITHMS
        % ///////////////////////////////////////////////////////////////////////////////
        function X = SENSE_CG(ObjMRI,arg,initialEstimate,RHS)
            
            opt.RepKspaceData = 1; %or 'inf' for reconstruction of all Reps
            opt.niter = 10;
            % opt.ReconUnderSampledkSpace = 0;
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.display = 0;
            opt.save = 0;
            
            arg.null = 0;
            
            opt = getFiledsFromUsersOpt(opt,arg);
            
            
            if isinf(opt.RepKspaceData) % if inf, then reconstruct all Reps
                nRe = 1:ObjMRI.nReps;
            else
                nRe = opt.RepKspaceData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            if nargin>=3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize, length(nRe)],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize length(nRe)],'single');
                end
            end
            
            X = img;
            
            % some reports
            if ObjMRI.SR.is
                fprintf('High-resolution reconstruction\n');
            else
                fprintf('Native-resolution reconstruction\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('Undersampled data\n');
            else
                fprintf('Fully sampled\n');
            end
            
            if opt.MrRegularizationParameter, fprintf('%s Regualrization...\n',opt.MrPriorType); end
            
            for i=1:length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                % K-space data
                Data = ObjMRI.ftKspaceData(:,:,:,:,nRe(i)); % [nCol,nLin,nPar,nCoil,nReps]
                
                % solve Ax = b, where A = FHF + beta*DTwD, b = FH
                if nargin ==4
                    b = RHS(:,:,:,nRe(i));
                else
                    b = ObjMRI.FH(Data);
                end
                
                if opt.MrRegularizationParameter
                    if strcmpi(opt.MrPriorType,'Quadratic')
                        if size(opt.MrPreCompWeights,1) ==1
                            W = 1./ObjMRI.Prior.nS;
                        else
                            W = opt.MrPreCompWeights;
                        end
                        A = @(x, dummy)ObjMRI.FHF(x) + ...
                            opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
                    else % for joint weight calculation
                        A = @(x,W)ObjMRI.FHF(x) + ...
                            opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
                    end
                else
                    A = @(x,dummy)ObjMRI.FHF(x);
                end
                
                
                X(:,:,:,i) = ObjMRI.PCG(b,A,img(:,:,:,i),opt.niter,1,opt);
            end
        end
        
        function x = PCG(ObjMRI,a,A,x0,Nit,P,arg)
            
            % a: RHS of the normal equation, i.e. FH(mri data)
            % A: FHF() handel object
            % x0: initial image estimate
            % Nit: number of sense iterations
            % P: Precondictioner
            % opt: options
            
            opt.message =[];
            opt.MrPriorType = [];
            opt.save = 0;
            opt.display = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
            
            %{
if strcmpi(opt.MrPriorType,'joint')
    Wp = ObjMRI.Prior.W_JointEntropy(opt.PriorMrImage./max(opt.PriorMrImage(:)),opt.priorImgSigma);
elseif strcmpi(opt.MrPriorType,'self-guided-wQ')
    Wp = 1;
end
            %}
            
            W = 1;
            a = P.*a;
            u = x0./P;
            r = a - P.*A(P.*u,W);
            p = r;
            for i=1:Nit
                uOld = u;
                q = P.*A(P.*p,W);
                alpha = r(:)'*r(:)/(p(:)'*q(:));
                u = uOld + alpha*p;
                rnew = r - alpha*q;
                p = rnew + rnew(:)'*rnew(:)/(r(:)'*r(:))*p;
                r = rnew;
                x = P.*u;
                
                %{
    if strcmpi(opt.MrPriorType,'joint')
        W0 = Wp .* ObjMRI.Prior.W_GaussianKernel(abs(u)./max(abs(u(:))),opt.MrSigma);
        W = W0./repmat(sum(W0,2),[1,ObjMRI.Prior.nS]);
    end
                %}
                if opt.display
                    if ObjMRI.is3D
                        drawnow, imshow(abs(u(:,:,opt.display)),[])
                        title([opt.message ' Iteration: #' num2str(i)]),pause(0.1)
                    end
                end
                
            end
        end
        
        function [X, Report] = gradientDescent(ObjMRI,arg,initialEstimate)
            
            arg.null = 0;
            opt.stepSize = 0.2;
            opt.niter = 20;
            opt.RepKspaceData = 1; %or 'inf' for reconstruction of all Reps
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            
            if opt.display, figure; end
            if isinf(opt.RepKspaceData) % if inf, then reconstruct all Reps
                nRe = 1:ObjMRI.nReps;
            else
                nRe = opt.RepKspaceData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            if nargin==3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize, length(nRe)],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize length(nRe)],'single');
                end
            end
            
            X = img;
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,length(nRe),'single');
                Report.stepSize = zeros(opt.niter,length(nRe),'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
                Xp = X;
            else
                Report = [];
            end
            % some report ----------------------------
            if ObjMRI.SR.is
                fprintf('- High-resolution reconstruction.\n');
            else
                fprintf('- Native-resolution reconstruction.\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('- Undersampled data.\n');
            else
                fprintf('- Fully sampled.\n');
            end
            if opt.MrRegularizationParameter
                fprintf('- %s Regularization.\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            if opt.stepSizeOptimization
                fprintf('- Step size optimization\n');
                col = @(x) x(:);
                H = @(x,W) col(ObjMRI.FHF(x) + D(x,W));
            end
            
            if ~ObjMRI.SR.is, warning('SR is off'); end
            for i=1:length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                % K-space data
                Data = ObjMRI.ftKspaceData(:,:,:,:,nRe(i)); % [nCol,nLin,nPar,nCoil,nReps]
                
                for j =1:opt.niter
                    % get weights coeffients
                    if strcmpi(opt.MrPriorType , 'Quadratic')
                        W = opt.MrPreCompWeights;
                    elseif strcmpi(opt.MrPriorType , 'Joint')
                        x = X(:,:,:,i);
                        W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(x)./max(abs(x(:))),opt.fSigma);
                        W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                    else
                        W = 1;
                    end
                    % do GD update
                    dPhix =(ObjMRI.FH(ObjMRI.F(X(:,:,:,i))- Data) + D(X(:,:,:,i),W));
                    if opt.stepSizeOptimization
                        stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), H(dPhix,W));
                        Report.stepSize(j,i) = stepSize;
                    else
                        stepSize = opt.stepSize;
                    end
                    X(:,:,:,i) = X(:,:,:,i) - stepSize*dPhix;
                    
                    % display and report relative L2-norm between iterates
                    if opt.display
                        if ObjMRI.is3D
                            x = X(:,:,:,i);
                            drawnow, imshow(abs(x(:,:,opt.display)),[])
                            title([opt.message ' Iteration: #' num2str(j)]),pause(0.1)
                        end
                    end
                    if opt.report
                        Report.relativeError(j,i) = Norm(X(:,:,:,i),Xp(:,:,:,i));
                        Xp(:,:,:,i) = X(:,:,:,i);
                    end
                end
            end
        end
        
        function [X, Report] = gradientDescent_4D(ObjMRI,arg,initialEstimate)
            
            arg.null = 0;
            opt.stepSize = 0.2;
            opt.niter = 20;
            opt.RepKspaceData = 'all_labels'; %all_controls or a subset vector, e.g.[2,4,5...], of all Reps
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
            if strcmpi('all_labels',opt.RepKspaceData)
                nRe = 2:2:ObjMRI.nReps;
            elseif strcmpi('all_controls',opt.RepKspaceData)
                nRe = 3:2:ObjMRI.nReps;
            else
                nRe = opt.RepKspaceData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            if nargin==3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize],'single');
                end
            end
            
            X = img;
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,1,'single');
                Report.stepSize = zeros(opt.niter,1,'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
                Xp = X;
            else
                Report = [];
            end
            % some report ----------------------------
            if ObjMRI.SR.is
                fprintf('- High-resolution reconstruction.\n');
            else
                fprintf('- Native-resolution reconstruction.\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('- Undersampled data.\n');
            else
                fprintf('- Fully sampled.\n');
            end
            if opt.MrRegularizationParameter
                fprintf('- %s Regularization.\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            if opt.stepSizeOptimization
                fprintf('- Step size optimization\n');
                col = @(x) x(:);
                H = @(x,W) col(ObjMRI.FHF(x) + D(x,W));
            end
            
            if ~ObjMRI.SR.is, warning('SR is off'); end
            
            dPhix = 0;
            
            for j =1:opt.niter
                % get weights coeffients
                if strcmpi(opt.MrPriorType , 'Quadratic')
                    W = opt.MrPreCompWeights;
                elseif strcmpi(opt.MrPriorType , 'Joint')
                    W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(X)./max(abs(X(:))),opt.fSigma);
                    W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                else
                    W = 1;
                end
                
                for i=1:length(nRe)
                    fprintf('Rep: #%d\n',nRe(i))
                    % do GD update
                    dPhix = dPhix + ObjMRI.FH(ObjMRI.F(X)- ObjMRI.ftKspaceData(:,:,:,:,nRe(i))) ; % [nCol,nLin,nPar,nCoil,nReps]
                end
                
                dPhix = dPhix/length(nRe) + D(X,W);
                
                if opt.stepSizeOptimization
                    stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), H(dPhix,W));
                    Report.stepSize(j,i) = stepSize;
                else
                    stepSize = opt.stepSize;
                end
                X = X - stepSize*dPhix;
                
                % display and report relative L2-norm between iterates
                if opt.display
                    if ObjMRI.is3D
                        drawnow, imshow(abs(X(:,:,opt.display)),[])
                        title([opt.message ' Iteration: #' num2str(j)]),pause(0.1)
                    end
                end
                if opt.report
                    Report.relativeError(j) = Norm(X,Xp);
                    Xp = X;
                end
            end
        end
        
        function [X, M0,Report] = gradientDescent_4D_diff(ObjMRI,arg,initialEstimate)
            
            arg.null = 0;
            opt.stepSize = 0.2;
            opt.niter = 20;
            opt.RepKspaceData = [3,4]; %[n,m], first label index, and first control indenx
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
            nRe = (ObjMRI.nReps -2)/2;
            DataL  = ObjMRI.ftKspaceData(:,:,:,:,opt.RepKspaceData(1):2:ObjMRI.nReps-1);
            DataC = ObjMRI.ftKspaceData(:,:,:,:,opt.RepKspaceData(2):2:ObjMRI.nReps);
            Data = DataC - DataL;
            DataM0 = ObjMRI.ftKspaceData(:,:,:,:,1);
            clear DataL DataC
            
            if nargin==3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize],'single');
                end
            end
            
            X = img;
            M0 = img;
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,1,'single');
                Report.stepSize = zeros(opt.niter,1,'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
                Xp = X;
            else
                Report = [];
            end
            % some report ----------------------------
            if ObjMRI.SR.is
                fprintf('- High-resolution reconstruction.\n');
            else
                fprintf('- Native-resolution reconstruction.\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('- Undersampled data.\n');
            else
                fprintf('- Fully sampled.\n');
            end
            if opt.MrRegularizationParameter
                fprintf('- %s Regularization.\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            if opt.stepSizeOptimization
                fprintf('- Step size optimization\n');
                col = @(x) x(:);
                H = @(x,W) col(ObjMRI.FHF(x) + D(x,W));
            end
            
            if ~ObjMRI.SR.is, warning('SR is off'); end
            
            dPhix = 0;
            for j =1:opt.niter
                % get weights coeffients
                if strcmpi(opt.MrPriorType , 'Quadratic')
                    W = opt.MrPreCompWeights;
                elseif strcmpi(opt.MrPriorType , 'Joint')
                    W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(X)./max(abs(X(:))),opt.fSigma);
                    W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                else
                    W = 1;
                end
                
                
                for i=1:nRe
                    fprintf('Rep: #%d\n',i)
                    % do GD update
                    dPhix = dPhix + ObjMRI.FH(ObjMRI.F(X,i)- Data(:,:,:,:,i),i); % [nCol,nLin,nPar,nCoil,nReps]
                end
                dPhix = dPhix/nRe + D(X,W);
                
                % 1th motion state (i.e. M0's motion state) is default
                % reference motion state, hence jm was not argumented into FH
                % and F operators
                dPhixM0 = ObjMRI.FH(ObjMRI.F(M0)- DataM0) + D(M0,W); % [nCol,nLin,nPar,nCoil,nReps]
                
                if opt.stepSizeOptimization
                    stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), H(dPhix,W));
                    stepSizeM0 = dot(dPhixM0(:), dPhixM0(:))./ dot(dPhixM0(:), H(dPhixM0,W));
                    Report.stepSize(j) = stepSize;
                else
                    stepSize = opt.stepSize;
                    stepSizeM0 = opt.stepSize;
                end
                X = X - stepSize*dPhix;
                M0 = M0 - stepSizeM0*dPhixM0;
                % display and report relative L2-norm between iterates
                if opt.display
                    if ObjMRI.is3D
                        
                        subplot(121),imshow(abs(X(:,:,opt.display)),[])
                        drawnow
                        subplot(122),imshow(abs(M0(:,:,opt.display)),[])
                        title([opt.message ' Iteration: #' num2str(j)])
                    end
                end
                
                if opt.report
                    Report.relativeError(j) = Norm(X,Xp);
                    Xp = X;
                end
            end
        end
        
        
        function [imageOut,Report] = deconvolution(ObjMRI,Imgs,arg)
            arg.null = 0;
            opt.niter = 20;
            opt.RepImageData = 1; %or 'inf' for reconstruction of all Reps
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.optimizationMethod = 'gradientDescent';% 'LucyRichardson-OSL'; conjugateGradient
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.imgMask = 1;
            opt.stepSize = 0.2;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
            if isinf(opt.RepImageData) % if inf, then reconstruct all Reps
                nRe = 1:size(Imgs,4);
            else
                nRe = opt.RepImageData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,length(nRe),'single');
                Report.stepSize = zeros(opt.niter,length(nRe),'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
            else
                Report = [];
            end
            
            if ObjMRI.SR.is
                fprintf('High-resolution deconvolution\n');
                imageOut = zeros([ObjMRI.SR.imgRefH.ImageSize,length(nRe)],'single');
                %                 Ft = @(x) ObjMRI.UpSampling(x);
                %                 F = @(x) ObjMRI.DownSampling(x);
            else
                fprintf('Native-resolution deconvolution\n');
                imageOut = Imgs;
                %                 F = @(x) gauss3DFilter(x,ObjMRI.SR.imgRefL.pixdim,ObjMRI.SR.PSF);
                %                 Ft = @(x) F(x);
            end
            %
            if opt.MrRegularizationParameter
                fprintf('%s Regularization...\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            Ft = @(x) ObjMRI.UpSampling(x);
            F = @(x) ObjMRI.DownSampling(x);
            
            if opt.stepSizeOptimization
                col = @(x) x(:);
                H = @(x,W)Ft(F(x)) + D(x,W);
            end
            
            
            for i = 1: length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                Ft1 = Ft(ones(size(Imgs(:,:,:,nRe(i)))));
                Ftx = Ft(Imgs(:,:,:,nRe(i)));
                ImageOut = Ftx;
                if opt.report, ImageOutP = ImageOut; end
                for j = 1:opt.niter
                    
                    % get weights coeffients
                    if strcmpi(opt.MrPriorType , 'Quadratic')
                        W = opt.MrPreCompWeights;
                    elseif strcmpi(opt.MrPriorType , 'Joint')
                        W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(ImageOut)./max(abs(ImageOut(:))),opt.fSigma);
                        W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                    else
                        W = 1;
                    end
                    % optimization method
                    if strcmpi(opt.optimizationMethod, 'LucyRichardson')
                        ImageOut =  (ImageOut.*opt.imgMask./(Ft1+D(ImageOut,W)+1e-8)).*Ft(Imgs(:,:,:,nRe(i))./(F(ImageOut)+1e-8));
                    elseif strcmpi(opt.optimizationMethod, 'gradientDescent')
                        dPhix = Ft(F(ImageOut) - Imgs(:,:,:,nRe(i))) + D(ImageOut,W);
                        if opt.stepSizeOptimization
                            stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), col(H(dPhix,W)));
                            Report.stepSize(j,i) = stepSize;
                        else
                            stepSize = opt.stepSize;
                        end
                        ImageOut = ImageOut - stepSize*dPhix;
                    elseif strcmpi(opt.optimizationMethod, 'conjugateGradient')
                        [ImageOut] = ObjMRI.PCG(Ftx, H,ImageOut,opt.niter,1,opt);
                        % currently no report for PCG
                        Report =[];
                        break
                    else
                        error('unknown OptimizationMethod')
                    end
                    
                    % display and report relative L2-norm between iterates
                    if opt.display
                        if ObjMRI.is3D
                            drawnow, imshow(abs(ImageOut(:,:,opt.display)),[])
                            title([opt.message ' Iteration: #' num2str(j)]),pause(0.1)
                        end
                    end
                    if opt.report
                        Report.relativeError(j,i) = Norm(ImageOut,ImageOutP);
                        ImageOutP = ImageOut;
                    end
                    
                end
                imageOut(:,:,:,i) = ImageOut;
            end
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           PVC AND DENOISING METHODS
        % ///////////////////////////////////////////////////////////////////////////////
        
        % function mLTS3D %modified least trimmed squares
        function [perfu_gm,perfu_wm] = LR3D_2(ObjMRI, perfu, pvGm,pvWm,arg)
            
            opt.imCropFactor = [8,8,0];
            opt.kernelSize = 3;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            opt.sWindowSize = opt.kernelSize;
            ObjMRI.BuildNativeResolutionPrior(opt);
            
            I = prod(ObjMRI.Prior.CropedImageSize);
            [perfu_gm,perfu_wm] = deal(zeros([ObjMRI.Prior.CropedImageSize]));
            
            perfu = ObjMRI.Prior.imCrop(perfu);
            pvGm = ObjMRI.Prior.imCrop(pvGm);
            pvWm = ObjMRI.Prior.imCrop(pvWm);
            
            for i = 1 : I
                s = ObjMRI.Prior.SearchWindow(i,:);
                G =[pvGm(s)',pvWm(s)'];%,pvCsf(s)'
                nNonZeroNhb = numel(find(G(:)~=0));
                if nNonZeroNhb > (3*6)
                    G = pinv(G);
                    tmp = G*perfu(s)';
                    perfu_gm(i) = tmp(1);
                    perfu_wm(i) = tmp(2);
                    %         tissue_csf(i) = tmp(3);
                    
                end
            end
            perfu_gm = ObjMRI.Prior.UndoImCrop(perfu_gm);
            perfu_wm = ObjMRI.Prior.UndoImCrop(perfu_wm);
            
            perfu_gm(isinf(perfu_gm)) = 0;
            perfu_gm(isnan(perfu_gm)) = 0;
            perfu_gm = max(0,perfu_gm);
            
            perfu_wm(isinf(perfu_wm)) = 0;
            perfu_wm(isnan(perfu_wm)) = 0;
            perfu_wm = max(0,perfu_wm);
            
        end
        
        function ImageOut = nonLocalMeans(ObjMRI,Image,PriorImg,GaussianWeightSigma)
            Image = ObjMRI.Prior.imCrop(Image);
            if ~isempty(PriorImg) % if empty [], self-similarities
                % Calculate non-local similarity coeff from a prior image
                PriorImg = ObjMRI.Prior.imCrop(PriorImg);
            else
                % Calculate non-local similarity coeff from the image
                % itself
                PriorImg = Image;
            end
            w = ObjMRI.Prior.W_GaussianKernel(PriorImg,GaussianWeightSigma);
            w = w./repmat(sum(w,2),[1,ObjMRI.Prior.nS]);
            ImageOut = reshape(sum(Image(ObjMRI.Prior.SearchWindow(:,:)).*w,2),ObjMRI.Prior.CropedImageSize);
            ImageOut = ObjMRI.Prior.UndoImCrop(ImageOut);
        end
        
    end
    
end

