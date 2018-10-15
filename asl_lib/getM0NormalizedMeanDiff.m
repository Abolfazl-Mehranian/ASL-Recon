
function imgOut = getM0NormalizedMeanDiff(imgSeries,brainExtractMask, nReps,fwhm,pixdim)

imgOut = abs(mean(abs(imgSeries(:,:,:,4:2:nReps))-abs(imgSeries(:,:,:,3:2:nReps-1)),4));
temp = gauss3DFilter(imgOut,pixdim,fwhm);
imgOut = temp./abs(imgSeries(:,:,:,1)) .* brainExtractMask;
imgOut(isnan(imgOut)) = 0;
imgOut(isinf(imgOut)) = 0;