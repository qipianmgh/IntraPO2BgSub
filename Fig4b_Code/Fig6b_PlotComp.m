clear
clc
close all

BeamDiaRaw = zeros(15,5);
load BeamDiaSamllNA_Batch1
BeamDiaRaw(:,1) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_Batch2
BeamDiaRaw(:,2) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_Batch3
BeamDiaRaw(:,3) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_Batch4
BeamDiaRaw(:,4) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_Batch5
BeamDiaRaw(:,5) = BeamDiaSmallNA(:,2);

VortexCharge = BeamDiaSmallNA(:,1);
BeamDiaMean = mean(BeamDiaRaw,2);
BeamDiaStd = std(BeamDiaRaw,[],2);

load('BeamDiaSimu_Emi_NA_0p56.mat');


SaveFlag= 1;
h = figure;
h.Units ='inch';
h.Position=[2 2 2.42 2];
h.PaperPositionMode='auto';
ax_handle = axes(h, 'Position', [0.1672 0.1682 0.8 0.765]);

errorbar(VortexCharge([1:11,13:15]),BeamDiaMean([1:11,13:15]),BeamDiaStd([1:11,13:15]),'r-','linewidth',1);
hold on
plot(BeamDiaSimu_Emi([1:11,13:15],1),BeamDiaSimu_Emi([1:11,13:15],2),'b-^','linewidth',1);
h_xlabel = xlabel('Vortex beam topological charge','Fontname', 'Calibri','FontSize', 10);
h_ylabel = ylabel('Beam diameter (\mum)','Fontname', 'Calibri','FontSize', 10);
set(h_xlabel,'Position',[13.0000   -1.95   -1.0000]);
set(h_ylabel,'Position',[-2.8  9   -1.0000]);
h_legend= legend('Experiment','Simulation');
set(h_legend,'Fontname', 'Calibri','FontSize',9,'Position',[0.1807    0.7294    0.4647    0.1901]);
grid on
axis([0 26 0 18]);
set(gca,'Fontname', 'Calibri','FontSize', 10);
if SaveFlag
   print(h,'Fig3b_VortexBeamDiaMeanStd_Omit20_v8_Final','-dmeta');  
end

BeamDiaEmi_Diff(:,1) = BeamDiaSimu_Emi(:,1);
BeamDiaEmi_Diff(:,2) = BeamDiaMean - BeamDiaSimu_Emi(:,2);