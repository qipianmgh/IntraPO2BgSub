clear
clc
close all

%%%%%%%%%%%%%%%%   G16L0 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L0_223419_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.1353;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask = SurveyImgNormMaxTemp>=Threshold;
SurveyImgMask = bwareafilt(SurveyImgMask, 1, 8);
SurveyImgNormMaxTemp = SurveyImgNormMaxTemp.*SurveyImgMask;
SurveyImgNormSum = SurveyImgNormMaxTemp/sum(SurveyImgNormMaxTemp(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L0 Plot Transverse   %%%%%%%%%%%%%%%%
PixelNum_5um_X = 5/FOV*double(GalvoPara.PixelX);
y_PixelNum_5um = 230*ones(1,2*round(PixelNum_5um_X));
PixelNum_5um_Y = 5/FOV*double(GalvoPara.PixelY);
x_PixelNum_5um = 20*ones(1,2*round(PixelNum_5um_Y));

I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_0 = figure;
h_fig3a_0.Units ='inch';
h_fig3a_0.Position=[2 2 1.5 1.5];
h_fig3a_0.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_0, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
%%% Plot line
hold on
plot(linspace(20,20+PixelNum_5um_X,2*round(PixelNum_5um_X)),y_PixelNum_5um,'w-','Linewidth',3);
plot(x_PixelNum_5um,linspace(230-PixelNum_5um_Y,230,2*round(PixelNum_5um_Y)),'w-','Linewidth',3);
text(75,232,'X','Fontname', 'Calibri','FontSize',14,'FontWeight','bold','Color','w');
text(12,165,'Y','Fontname', 'Calibri','FontSize',14,'FontWeight','bold','Color','w');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_0,'-dtiffn','-r600','Fig3a_Tran_0.tif');

%%%%%%%%%%%%%%%%    G16L0 Plot Axial   %%%%%%%%%%%%%%%%
load('AngioStack_G16L0_223802_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

PixelNum_5um_Tran = 5/FOV*double(GalvoPara.PixelX);
y_PixelNum_5um_Tran = 230*ones(1,2*round(PixelNum_5um_Tran));
PixelNum_5um_Axial = 5/25*256;
x_PixelNum_5um_Axial = 20*ones(1,2*round(PixelNum_5um_Axial));

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_0 = figure;
h_fig3a_Axial_0.Units ='inch';
h_fig3a_Axial_0.Position=[2 2 1.5 1.5];
h_fig3a_Axial_0.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_0, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
%%% Plot line
hold on
plot(linspace(20,20+PixelNum_5um_Tran,2*round(PixelNum_5um_Tran)),y_PixelNum_5um_Tran,'w-','Linewidth',3);
plot(x_PixelNum_5um_Axial,linspace(230-PixelNum_5um_Axial,230,2*round(PixelNum_5um_Axial)),'w-','Linewidth',3)
text(75,232,'X','Fontname', 'Calibri','FontSize',14,'FontWeight','bold','Color','w');
text(12,165,'Z','Fontname', 'Calibri','FontSize',14,'FontWeight','bold','Color','w');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_0,'-dtiffn','-r600','Fig3a_Axial_0.tif');

%%%%%%%%%%%%%%%%   G16L5 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L5_230336_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.1353;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask(SurveyImgNormMaxTemp>=Threshold) = 1;
SurveyImgMask(SurveyImgNormMaxTemp<Threshold) = 0;
SurveyImgNormSum = SurveyImgMask/sum(SurveyImgMask(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L5 Plot Transverse   %%%%%%%%%%%%%%%% 
I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_5 = figure;
h_fig3a_5.Units ='inch';
h_fig3a_5.Position=[2 2 1.5 1.5];
h_fig3a_5.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_5, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_5,'-dtiffn','-r600','Fig3a_Tran_5.tif');

%%%%%%%%%%%%%%%%   G16L5 Plot Axial   %%%%%%%%%%%%%%%% 
load('AngioStack_G16L5_230650_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_5 = figure;
h_fig3a_Axial_5.Units ='inch';
h_fig3a_Axial_5.Position=[2 2 1.5 1.5];
h_fig3a_Axial_5.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_5, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_5,'-dtiffn','-r600','Fig3a_Axial_5.tif');

%%%%%%%%%%%%%%%%   G16L11 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L11_215548_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.05;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask(SurveyImgNormMaxTemp>=Threshold) = 1;
SurveyImgMask(SurveyImgNormMaxTemp<Threshold) = 0;
SurveyImgNormSum = SurveyImgMask/sum(SurveyImgMask(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L11 Plot Transverse   %%%%%%%%%%%%%%%% 
I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_11 = figure;
h_fig3a_11.Units ='inch';
h_fig3a_11.Position=[2 2 1.5 1.5];
h_fig3a_11.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_11, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_11,'-dtiffn','-r600','Fig3a_Tran_11.tif');

%%%%%%%%%%%%%%%%   G16L11 Plot Axial   %%%%%%%%%%%%%%%% 
load('AngioStack_G16L11_215921_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_11 = figure;
h_fig3a_Axial_11.Units ='inch';
h_fig3a_Axial_11.Position=[2 2 1.5 1.5];
h_fig3a_Axial_11.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_11, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_11,'-dtiffn','-r600','Fig3a_Axial_11.tif');

%%%%%%%%%%%%%%%%   G16L15 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L15_214144_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.05;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask(SurveyImgNormMaxTemp>=Threshold) = 1;
SurveyImgMask(SurveyImgNormMaxTemp<Threshold) = 0;
SurveyImgNormSum = SurveyImgMask/sum(SurveyImgMask(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L15 Plot Transverse   %%%%%%%%%%%%%%%% 
I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_15 = figure;
h_fig3a_15.Units ='inch';
h_fig3a_15.Position=[2 2 1.5 1.5];
h_fig3a_15.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_15, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_15,'-dtiffn','-r600','Fig3a_Tran_15.tif');

%%%%%%%%%%%%%%%%   G16L15 Plot Axial   %%%%%%%%%%%%%%%% 
load('AngioStack_G16L15_214504_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_15 = figure;
h_fig3a_Axial_15.Units ='inch';
h_fig3a_Axial_15.Position=[2 2 1.5 1.5];
h_fig3a_Axial_15.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_15, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_15,'-dtiffn','-r600','Fig3a_Axial_15.tif');

%%%%%%%%%%%%%%%%   G16L20 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L20_203752_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.05;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask(SurveyImgNormMaxTemp>=Threshold) = 1;
SurveyImgMask(SurveyImgNormMaxTemp<Threshold) = 0;
SurveyImgNormSum = SurveyImgMask/sum(SurveyImgMask(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L20 Plot Transverse   %%%%%%%%%%%%%%%% 
I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_20 = figure;
h_fig3a_20.Units ='inch';
h_fig3a_20.Position=[2 2 1.5 1.5];
h_fig3a_20.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_20, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_20,'-dtiffn','-r600','Fig3a_Tran_20.tif');

%%%%%%%%%%%%%%%%   G16L20 Plot Axial   %%%%%%%%%%%%%%%% 
load('AngioStack_G16L20_204156_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_20 = figure;
h_fig3a_Axial_20.Units ='inch';
h_fig3a_Axial_20.Position=[2 2 1.5 1.5];
h_fig3a_Axial_20.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_20, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_20,'-dtiffn','-r600','Fig3a_Axial_20.tif');

%%%%%%%%%%%%%%%%   G16L25 Centroid   %%%%%%%%%%%%%%%%
load('AngioSurveyScan_G16L25_201036_0um.mat');
FOV = double(GalvoPara.FOV);
SurveyImg = double(fliplr(AI_SaveData));
SurveyImg(SurveyImg<0)=0;
SurveyImgFlt = medfilt2(SurveyImg,[2 2]);
SurveyImgNormMax = SurveyImgFlt/max(SurveyImgFlt(:));
Threshold = 0.05;
SurveyImgNormMaxTemp = SurveyImgNormMax; 
SurveyImgMask(SurveyImgNormMaxTemp>=Threshold) = 1;
SurveyImgMask(SurveyImgNormMaxTemp<Threshold) = 0;
SurveyImgNormSum = SurveyImgMask/sum(SurveyImgMask(:));
[row,col] = size(SurveyImgNormSum);
[X,Y] = ndgrid(1:row,1:col);
row_centroid = sum(sum(X.*SurveyImgNormSum));
col_centroid = sum(sum(Y.*SurveyImgNormSum));
%%%%%%%%%%   Plot Centroid
h1 = figure;
imagesc(SurveyImgFlt);
hold on
plot(col_centroid,row_centroid,'g*');
colormap jet

%%%%%%%%%%%%%%%%   G16L25 Plot Transverse   %%%%%%%%%%%%%%%% 
I0 = SurveyImgFlt;
I0 = I0/max(I0(:));
h_fig3a_25 = figure;
h_fig3a_25.Units ='inch';
h_fig3a_25.Position=[2 2 1.5 1.5];
h_fig3a_25.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_25, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_25,'-dtiffn','-r600','Fig3a_Tran_25.tif');

%%%%%%%%%%%%%%%%   G16L25 Plot Axial   %%%%%%%%%%%%%%%% 
load('AngioStack_G16L25_201455_0_24_1um.mat');
FOV = double(GalvoPara.FOV);
for i = 1:size(AI_SaveData,3)
    AI_SaveData(:,:,i) = fliplr(AI_SaveData(:,:,i));
end
StackImg = double(AI_SaveData);
StackImg(StackImg<0) = 0;
StackImgFlt = medfilt3(StackImg,[3 3 3]); 
StackImgMid = squeeze(StackImgFlt(round(row_centroid),:,:));
StackImgMid = StackImgMid';
[X,Y] = meshgrid(-127.5:127.5,-12:12);
[Xq,Yq] = meshgrid(-127.5:127.5,linspace(-12,12,256)); 
StackImgMidInterp = interp2(X,Y,StackImgMid,Xq,Yq,'linear');

I0 = StackImgMidInterp;
I0 = I0/max(I0(:));
h_fig3a_Axial_25 = figure;
h_fig3a_Axial_25.Units ='inch';
h_fig3a_Axial_25.Position=[2 2 1.5 1.5];
h_fig3a_Axial_25.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_Axial_25, 'Position', [0 0 1 1]);
imagesc(I0,[0 1]);
colormap jet;
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');
axis tight
axis equal
axis off
print(h_fig3a_Axial_25,'-dtiffn','-r600','Fig3a_Axial_25.tif');

%%% Plot colorbar only
h_fig3a_00 = figure;
h_fig3a_00.Units ='inch';
h_fig3a_00.Position=[2 2 2.62 2];
h_fig3a_00.PaperPositionMode='auto';
ax_handle = axes(h_fig3a_00, 'Position', [0 0.05 1 0.82]);
imagesc(I0,[0 1]);
colormap jet;
h_colorbar = colorbar;
set(h_colorbar,'Fontname', 'Calibri','FontSize',8,'FontWeight','bold','Position',[0.85 0.05 0.03 0.82]);
axis tight
axis equal
axis off
print(h_fig3a_00,'Fig3a_Colorbar','-dmeta');  