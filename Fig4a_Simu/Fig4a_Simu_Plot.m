clear
clc
close all

FOV = 25;    %%% unit: um
depth = 500;
z_limit = 12;
spa_pts = 256;
stackdepth = linspace(depth-z_limit,depth+z_limit,spa_pts);
stacksize = length(stackdepth);
FocalPlaneInd = 128;
AxialPlaneColInd = 128;
PixelNumX = 255;
PixelNumY = 255;
%%%%%%%%%%%%%%%%   G16L0 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack0_depth500um.mat');
SurveyImg = squeeze(Fluostack0(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
PixelNum_5um_X = 5/FOV*PixelNumX;
y_PixelNum_5um = 230*ones(1,2*round(PixelNum_5um_X));
PixelNum_5um_Y = 5/FOV*PixelNumY;
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
print(h_fig3a_0,'-dtiffn','-r600','Fig3a_Tran_0_Emi.tif');

%%%%%%%%%%%%%%%%    G16L0 Plot Axial   %%%%%%%%%%%%%%%%
StackImgFlt = Fluostack0;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

PixelNum_5um_Tran = 5/FOV*PixelNumX;
y_PixelNum_5um_Tran = 230*ones(1,2*round(PixelNum_5um_Tran));
PixelNum_5um_Axial = 5/25*256;
x_PixelNum_5um_Axial = 20*ones(1,2*round(PixelNum_5um_Axial));

I0 = StackImgMid;
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
print(h_fig3a_Axial_0,'-dtiffn','-r600','Fig3a_Axial_0_Emi.tif');

%%%%%%%%%%%%%%%%   G16L5 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack5_depth500um.mat');
SurveyImg = squeeze(Fluostack5(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
print(h_fig3a_5,'-dtiffn','-r600','Fig3a_Tran_5_Emi.tif');

%%%%%%%%%%%%%%%%   G16L5 Plot Axial   %%%%%%%%%%%%%%%% 
StackImgFlt = Fluostack5;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

I0 = StackImgMid;
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
print(h_fig3a_Axial_5,'-dtiffn','-r600','Fig3a_Axial_5_Emi.tif');

%%%%%%%%%%%%%%%%   G16L11 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack11_depth500um.mat');
SurveyImg = squeeze(Fluostack11(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
print(h_fig3a_11,'-dtiffn','-r600','Fig3a_Tran_11_Emi.tif');

%%%%%%%%%%%%%%%%   G16L11 Plot Axial   %%%%%%%%%%%%%%%% 
StackImgFlt = Fluostack11;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

I0 = StackImgMid;
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
print(h_fig3a_Axial_11,'-dtiffn','-r600','Fig3a_Axial_11_Emi.tif');

%%%%%%%%%%%%%%%%   G16L15 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack15_depth500um.mat');
SurveyImg = squeeze(Fluostack15(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
print(h_fig3a_15,'-dtiffn','-r600','Fig3a_Tran_15_Emi.tif');

%%%%%%%%%%%%%%%%   G16L15 Plot Axial   %%%%%%%%%%%%%%%% 
StackImgFlt = Fluostack15;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

I0 = StackImgMid;
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
print(h_fig3a_Axial_15,'-dtiffn','-r600','Fig3a_Axial_15_Emi.tif');

%%%%%%%%%%%%%%%%   G16L20 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack20_depth500um.mat');
SurveyImg = squeeze(Fluostack20(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
print(h_fig3a_20,'-dtiffn','-r600','Fig3a_Tran_20_Emi.tif');

%%%%%%%%%%%%%%%%   G16L20 Plot Axial   %%%%%%%%%%%%%%%% 
StackImgFlt = Fluostack20;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

I0 = StackImgMid;
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
print(h_fig3a_Axial_20,'-dtiffn','-r600','Fig3a_Axial_20_Emi.tif');

%%%%%%%%%%%%%%%%   G16L25 Centroid   %%%%%%%%%%%%%%%%
load('Fluostack25_depth500um.mat');
SurveyImg = squeeze(Fluostack25(:,:,FocalPlaneInd));
SurveyImgFlt = SurveyImg;
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
print(h_fig3a_25,'-dtiffn','-r600','Fig3a_Tran_25_Emi.tif');

%%%%%%%%%%%%%%%%   G16L25 Plot Axial   %%%%%%%%%%%%%%%% 
StackImgFlt = Fluostack25;
StackImgMid = squeeze(StackImgFlt(AxialPlaneColInd,:,:));
StackImgMid = StackImgMid';

I0 = StackImgMid;
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
print(h_fig3a_Axial_25,'-dtiffn','-r600','Fig3a_Axial_25_Emi.tif');

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
print(h_fig3a_00,'Fig3a_Colorbar_Emi','-dmeta');  