clear
clc
close all

listing = dir(pwd);
filenames = {listing.name}';

fileindex1 = ~cellfun('isempty',strfind(filenames,'AngioSurveyScan_G16L0_'));
fileindex2 = ~cellfun('isempty',strfind(filenames,'.mat'));
fileindex3 = find(fileindex1 & fileindex2,1,'last');
load(filenames{fileindex3});
FOV = double(GalvoPara.FOV);
AI_SaveData = double(AI_SaveData);
AI_SaveData(AI_SaveData<0)=0;
AI_SaveData = medfilt2(AI_SaveData,[2 2]);
I = fliplr(AI_SaveData(128-90:128+90,128-90:128+90)); 
I_norm_max = I/max(I(:));
I_mask = I_norm_max>=0.1353;
I_mask = imfill(I_mask, 'holes');
I_mask = bwareafilt(I_mask, 1, 8);
subplot(1,2,1);
imshow(I_mask);
props = regionprops(I_mask, 'Area', 'EquivDiameter', 'Centroid');
area = props.Area;
dia_pixel = props.EquivDiameter;
boundaries = bwboundaries(I_mask);
subplot(1,2,2);
imshow(I_norm_max);
for p = 1 : length(boundaries)
	boundary_coordinates = boundaries{p};
	x = boundary_coordinates(:, 2);
	y = boundary_coordinates(:, 1);
	hold on;
	plot(x, y, 'g-', 'LineWidth', 2);
end
viscircles([props.Centroid(1), props.Centroid(2)], dia_pixel/2);
BeamDia_Gaussian = dia_pixel*FOV/double(GalvoPara.PixelX);

PSF_Final = BeamDia_Gaussian/1.699;
save('PSF_Final_Batch1.mat','PSF_Final');


