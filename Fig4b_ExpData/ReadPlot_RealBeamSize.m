clear
clc
close all

listing = dir(pwd);
filenames = {listing.name}';
   
VortexCharge = [0 1 3 5 7 9 11 13 15 17 19 20 21 23 25];
BeamDia = zeros(size(VortexCharge));

for m = 1:length(VortexCharge)
    if VortexCharge(m) == 0 
       %%%%%%%%%%%%%%   Vortex charge = 0
       fileindex1 = ~cellfun('isempty',strfind(filenames,['AngioSurveyScan_G16L',num2str(VortexCharge(m)),'_']));
       fileindex2 = ~cellfun('isempty',strfind(filenames,'.mat'));
       fileindex3 = find(fileindex1 & fileindex2,1,'last');
       load(filenames{fileindex3});
       FOV = double(GalvoPara.FOV);
       AI_SaveData = double(AI_SaveData);
       AI_SaveData(AI_SaveData<0)=0;
       AI_SaveData = medfilt2(AI_SaveData,[2 2]);
       I = fliplr(AI_SaveData); 
       I_norm_max = I/max(I(:));
       I_mask = I_norm_max>=0.1353;
       I_mask = imfill(I_mask, 'holes');
       I_mask = bwareafilt(I_mask, 1, 8);
       props = regionprops(I_mask, 'Area', 'EquivDiameter', 'Centroid');
       dia_pixel = props.EquivDiameter;
       boundaries = bwboundaries(I_mask);
       BeamDiaTest = dia_pixel*FOV/256;
       BeamDia(m) = BeamDiaTest;
       
       figure;
       imshow(I_norm_max);
       for p = 1 : length(boundaries)
           boundary_coordinates = boundaries{p};
           x = boundary_coordinates(:, 2);
           y = boundary_coordinates(:, 1);
           hold on;
           plot(x, y, 'g-', 'LineWidth', 2);
       end
       viscircles([props.Centroid(1), props.Centroid(2)], dia_pixel/2);
    else
       %%%%%%%%%%%%%%   Vortex charge > 0
       fileindex1 = ~cellfun('isempty',strfind(filenames,['AngioSurveyScan_G16L',num2str(VortexCharge(m)),'_']));
       fileindex2 = ~cellfun('isempty',strfind(filenames,'.mat'));
       fileindex3 = find(fileindex1 & fileindex2,1,'last');
       load(filenames{fileindex3});
       FOV = double(GalvoPara.FOV);
       AI_SaveData = double(AI_SaveData);
       AI_SaveData(AI_SaveData<0)=0;
       AI_SaveData = medfilt2(AI_SaveData,[2 2]);
       I = fliplr(AI_SaveData); 
       I_norm_max = I/max(I(:));
       %%% Find centroid
       I_mask0 = I_norm_max>=0.01;
       I_mask0 = bwareafilt(I_mask0, 1, 8);
       I0 = I.*I_mask0;
       I_norm_sum0 = double(I_mask0)/sum(I_mask0(:));
       [row0,col0] = size(I_norm_sum0);
       [X0,Y0] = ndgrid(1:row0,1:col0);
       row_centroid = sum(sum(X0.*I_norm_sum0));
       col_centroid = sum(sum(Y0.*I_norm_sum0));
       %%% Find beam diameter
       I_mask = I_norm_max>=0.01;
       I_mask = bwareafilt(I_mask, 1, 8);
       I = I.*I_mask;
       ind = find(I_mask == 1);
       [row_beam col_beam]= ind2sub(size(I),ind);
       BeamDia(m) = 2*FOV/double(GalvoPara.PixelX)*mean(mean(sqrt((row_beam-row_centroid).^2+(col_beam-col_centroid).^2)));

       figure;
       imagesc(I);
       hold on
       plot(col_centroid,row_centroid,'g*');
       colormap jet        
    end

end

BeamDiaSmallNA = [VortexCharge' BeamDia'];
save('BeamDiaSamllNA_Batch1','BeamDiaSmallNA');
