clear
clc
close all

PSF_SizeRaw = zeros(1,5);
load PSF_Final_Batch1
PSF_SizeRaw(1,1) = PSF_Final;
load PSF_Final_Batch2
PSF_SizeRaw(1,2) = PSF_Final;
load PSF_Final_Batch3
PSF_SizeRaw(1,3) = PSF_Final;
load PSF_Final_Batch4
PSF_SizeRaw(1,4) = PSF_Final;
load PSF_Final_Batch5
PSF_SizeRaw(1,5) = PSF_Final;


PSF_SizeMean = mean(PSF_SizeRaw);
PSF_SizeStd = std(PSF_SizeRaw);

save('PSF_Final.mat','PSF_SizeRaw','PSF_SizeMean','PSF_SizeStd');

