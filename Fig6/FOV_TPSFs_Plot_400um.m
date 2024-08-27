clear
clc
close all

load pO2_SurveyScan_205856_400um.mat
load PixelIndex_Final_400um
load TPSF_400um_IA_Pixel55

survey_image = fliplr(DaqPara.SurveyScanImg);
survey_max = min(max(survey_image(:)),80);
survey_min = min(survey_image(:));

X_PixelNum = GalvoPara.PixelX;    
FOV = GalvoPara.FOV;
ScaleBarLen = 100/FOV*X_PixelNum;
PixelPos = [GalvoPara.PO2Pos(:,2) X_PixelNum+1-GalvoPara.PO2Pos(:,1)];

pixel_pos = [GalvoPara.PO2Pos(PixelIndex_Final,2) X_PixelNum+1-GalvoPara.PO2Pos(PixelIndex_Final,1)];
pixel_pos_sub = [GalvoPara.PO2Pos(PixelIndex_Final+1,2) X_PixelNum+1-GalvoPara.PO2Pos(PixelIndex_Final+1,1)];
tpsf_pixel_index = [55];

for p = 1:length(tpsf_pixel_index)

tpsf_pixel_pos = [GalvoPara.PO2Pos(tpsf_pixel_index(p),2) X_PixelNum+1-GalvoPara.PO2Pos(tpsf_pixel_index(p),1)];

h1 = figure;
linewidth = 3;
fontsize = 10;
legendfontsize = 8;
h1.Units ='inch';
h1.Position=[2 2 2.3 2.3];
h1.PaperPositionMode='auto';
ax_handle = axes(h1, 'Position', [0 0 1 1]);

imagesc(survey_image,[survey_min survey_max]);
colormap gray
hold on
plot(pixel_pos(:,2),pixel_pos(:,1),'Color',[0.9290 0.6940 0.1250],'LineStyle','none','Marker','.','MarkerSize',10);
plot(pixel_pos_sub(:,2),pixel_pos_sub(:,1),'Color','g','LineStyle','none','Marker','^','MarkerSize',3);
plot(tpsf_pixel_pos(:,2),tpsf_pixel_pos(:,1),'Color','r','LineStyle','none','Marker','o','MarkerSize',6);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
ScaleBarX = linspace(178,178+ScaleBarLen,100);
ScaleBarY = 225*ones(size(ScaleBarX));
plot(ScaleBarX,ScaleBarY,'w-','LineWidth',linewidth);
h_legend = legend('Intravas. Pixel','Ref. Pixel');
if p==4
   set(h_legend,'Fontname', 'Calibri','FontSize',legendfontsize,'Location','northeast','Position',[0.49 0.8560 0.5063 0.1436]);
else
   set(h_legend,'Fontname', 'Calibri','FontSize',legendfontsize,'Location','northeast','Position',[0.065 0.01 0.41 0.14]);
end
print(h1,'-dtiffn','-r600',['FOV_400um_Pixel_',num2str(tpsf_pixel_index(p)),'_v4']);
print(h1,['FOV_400um_Pixel_',num2str(tpsf_pixel_index(p)),'_v4'],'-dmeta'); 

h2 = figure;
h2.Units ='inch';
h2.Position=[8 2 2.62 2.3];
h2.PaperPositionMode='auto';
ax_handle2 = axes(h2, 'Position', [0.16 0.15 0.80 0.80]);

timepts = (0:2:298);
fittedDecayTimePoints_indices = 9:150;
linewidth2 = 1;
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(1,fittedDecayTimePoints_indices),'r-','LineWidth',linewidth2);
hold on
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(2,fittedDecayTimePoints_indices),'g-','LineWidth',linewidth2);
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(3,fittedDecayTimePoints_indices),'b-','LineWidth',linewidth2);
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(4,fittedDecayTimePoints_indices),'c-','LineWidth',linewidth2);
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(5,fittedDecayTimePoints_indices),'m-','LineWidth',linewidth2);
semilogy(timepts(fittedDecayTimePoints_indices),TPSF_400um_IA_Pixel(6,fittedDecayTimePoints_indices),'k-','LineWidth',linewidth2);
grid on
ylim([1 3e3]);

%%% Create legend
plot([120 140],[2e3 2e3],'r-','LineWidth',linewidth2);
text(145, 2e3,'Intravas. pixel (\itm \rm= 0)','Fontname', 'Calibri','FontSize',legendfontsize);
plot([120 140],[1.2e3 1.2e3],'g-','LineWidth',linewidth2);
text(145, 1.2e3,'Ref. pixel (\itm \rm= 0)','Fontname', 'Calibri','FontSize',legendfontsize);
plot([120 140],[720 720],'b-','LineWidth',linewidth2);
text(145, 720,'Intravas. pixel (\itm \rm= 11)','Fontname', 'Calibri','FontSize',legendfontsize);
plot([120 140],[432 432],'c-','LineWidth',linewidth2);
text(145, 432,'Intravas. pixel (\itm \rm= 13)','Fontname', 'Calibri','FontSize',legendfontsize);
plot([120 140],[259 259],'m-','LineWidth',linewidth2);
text(145, 259,'Intravas. pixel (\itm \rm= 15)','Fontname', 'Calibri','FontSize',legendfontsize);
plot([120 140],[155 155],'k-','LineWidth',linewidth2);
text(145, 155,'Intravas. pixel (\itm \rm= 20)','Fontname', 'Calibri','FontSize',legendfontsize);
rectangle('Position',[115 110 182 2700]);

h_labelx = xlabel('Time (\mus)','Fontname', 'Calibri','FontSize',10,'Position',[150.0001  0.47   -1.0]);
h_labely = ylabel('Photon counts (A.U.)','Fontname', 'Calibri','FontSize',10,'Position',[-34  55   -1.0]);
set(ax_handle2,'Fontname', 'Calibri','FontSize',10);
print(h2,['TPSF_400um_Pixel_',num2str(tpsf_pixel_index(p)),'_v4'],'-dmeta');
end