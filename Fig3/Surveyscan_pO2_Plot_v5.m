clear
clc
close all

load pO2_Galvo_TD_192837_0um.mat

SurveyImage = fliplr(CI_SurveyScanImg);
SurveyImage = medfilt2(SurveyImage, [2 2],'symmetric');
SurveyIntMax = max(SurveyImage(:));
SurveyIntMin = min(SurveyImage(:));

X_PixelNum = GalvoPara.PixelX;    
FOV = GalvoPara.FOV;
ScaleBarLen = 100/FOV*X_PixelNum;
PixelPos = [GalvoPara.PO2Pos(:,2) X_PixelNum+1-GalvoPara.PO2Pos(:,1)];

ArteryIndex = [9 13];
VeinIndex = [2 3];

Total_PC_Artery = sum(CI_DataProcessed(ArteryIndex,:),2);
Total_PC_Vein = sum(CI_DataProcessed(VeinIndex,:),2);
Total_PC_Mean = mean([Total_PC_Artery; Total_PC_Vein]); 

Peak_PC_Artery = max(CI_DataProcessed(ArteryIndex,:),[],2);
Peak_PC_Vein = max(CI_DataProcessed(VeinIndex,:),[],2);
Peak_PC_Mean = mean([Peak_PC_Artery; Peak_PC_Vein]);


h1 = figure;
h1.Units ='inch';
h1.Position=[2 2 2.62 2.3];
h1.PaperPositionMode='auto';
ax_handle = axes(h1, 'Position', [0.11 0.10 0.72 0.83]);
imagesc(SurveyImage,[0 80]);
colormap gray
h_bar = colorbar;
set(h_bar,'Position',[0.84   0.10    0.0541    0.83],'Fontname', 'Calibri','FontSize',9);
text(252,-10,'P.C.','Fontname', 'Calibri','FontSize',9);
hold on
plot(PixelPos(ArteryIndex(1),2),PixelPos(ArteryIndex(1),1),'r.','MarkerSize',15);
plot(PixelPos(ArteryIndex(2),2),PixelPos(ArteryIndex(2),1),'m.','MarkerSize',15);
plot(PixelPos(VeinIndex(1),2),PixelPos(VeinIndex(1),1),'c.','MarkerSize',15);
plot(PixelPos(VeinIndex(2),2),PixelPos(VeinIndex(2),1),'b.','MarkerSize',15);
text(PixelPos(ArteryIndex(1),2)-5,PixelPos(ArteryIndex(1),1)-15,'1','Color','r','Fontname', 'Calibri','FontSize',10);
text(PixelPos(ArteryIndex(2),2)-5,PixelPos(ArteryIndex(2),1)+15,'2','Color','m','Fontname', 'Calibri','FontSize',10);
text(PixelPos(VeinIndex(1),2)+10,PixelPos(VeinIndex(1),1)+8,'3','Color','c','Fontname', 'Calibri','FontSize',10);
text(PixelPos(VeinIndex(2),2)+8,PixelPos(VeinIndex(2),1),'4','Color','b','Fontname', 'Calibri','FontSize',10);
ScaleBarX = linspace(20,20+ScaleBarLen,100);
ScaleBarY = 230*ones(size(ScaleBarX));
plot(ScaleBarX,ScaleBarY,'w-','LineWidth',3);

set(ax_handle,'Fontname', 'Calibri','FontSize',9,'XTick',[1 250],'YTick',[1 250]);
print(h1,'Surveyscan_0um_v5','-dmeta');  

h2 = figure;
h2.Units ='inch';
h2.Position=[10 2 2.62 2.3];
h2.PaperPositionMode='auto';
ax_handle2 = axes(h2, 'Position', [0.16 0.15 0.80 0.80]);

load('dataStruct_delay5.mat');
pO2 = dataStruct.pO2;
tau = dataStruct.lifetimes;
fittedDecayTimePoints_indices = dataStruct.fittedDecayTimePoints_indices;
fittedDecayTimePoints_us = dataStruct.fittedDecayTimePoints_us;
originalDecay = dataStruct.summed_lifetime_data;
t = dataStruct.collectionDurationTimePoints_us;
lifetimefitdata = (tau(:,1).*exp(-fittedDecayTimePoints_us./tau(:,2))+tau(:,3)) .* originalDecay (:,fittedDecayTimePoints_indices(1));

semilogy(t(fittedDecayTimePoints_indices),lifetimefitdata(ArteryIndex(1),:),'r-','LineWidth',0.5);
hold on
semilogy(t,originalDecay(ArteryIndex(1),:),'r.','MarkerSize',5,'HandleVisibility','off');
semilogy(t(fittedDecayTimePoints_indices),lifetimefitdata(ArteryIndex(2),:),'m-','LineWidth',0.5);
semilogy(t,originalDecay(ArteryIndex(2),:),'m.','MarkerSize',5,'HandleVisibility','off');

semilogy(t(fittedDecayTimePoints_indices),lifetimefitdata(VeinIndex(1),:),'c-','LineWidth',0.5);
semilogy(t,originalDecay(VeinIndex(1),:),'c.','MarkerSize',5,'HandleVisibility','off');

semilogy(t(fittedDecayTimePoints_indices),lifetimefitdata(VeinIndex(2),:),'b-','LineWidth',0.5);
semilogy(t,originalDecay(VeinIndex(2),:),'b.','MarkerSize',5,'HandleVisibility','off');
axis([0 300 5 2e4]);
grid on
Line_X = linspace(145,170,100); 
Line1_Y = 8e3*ones(size(Line_X));
Line2_Y = 1.7e3*ones(size(Line_X));
Line3_Y = 3.5e2*ones(size(Line_X));
Line4_Y = 7e1*ones(size(Line_X));

plot(Line_X,Line1_Y,'r-','LineWidth',0.5);
plot(Line_X,Line2_Y,'m-','LineWidth',0.5);
plot(Line_X,Line3_Y,'c-','LineWidth',0.5);
plot(Line_X,Line4_Y,'b-','LineWidth',0.5);

text(175,8e3,['\tau_1= ',num2str(tau(ArteryIndex(1),2),'%.1f'),' \mus',newline,'pO_2= ',num2str(pO2(ArteryIndex(1)),'%.1f'),' mmHg'],'Color','k','Fontname', 'Calibri','FontSize',8);
text(175,1.7e3,['\tau_2= ',num2str(tau(ArteryIndex(2),2),'%.1f'),' \mus',newline,'pO_2= ',num2str(pO2(ArteryIndex(2)),'%.1f'),' mmHg'],'Color','k','Fontname', 'Calibri','FontSize',8);
text(175,3.5e2,['\tau_3= ',num2str(tau(VeinIndex(1),2),'%.1f'),' \mus',newline,'pO_2= ',num2str(pO2(VeinIndex(1)),'%.1f'),' mmHg'],'Color','k','Fontname', 'Calibri','FontSize',8);
text(175,7e1,['\tau_4= ',num2str(tau(VeinIndex(2),2),'%.1f'),' \mus',newline,'pO_2= ',num2str(pO2(VeinIndex(2)),'%.1f'),' mmHg'],'Color','k','Fontname', 'Calibri','FontSize',8);

h_labelx = xlabel('Time (\mus)','Fontname', 'Calibri','FontSize',9,'Position',[150 2.2 -1]);
h_labely = ylabel('Photon counts','Fontname', 'Calibri','FontSize',9,'Position',[-35  219 -1]);
set(ax_handle2,'Fontname', 'Calibri','FontSize',9);
print(h2,'TPSF_Fitting_0um_v5','-dmeta');  
