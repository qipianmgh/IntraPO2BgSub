clear
clc
close all

BeamDiaRaw = zeros(15,5);
load BeamDiaSamllNA_123122_Batch1
BeamDia(:,1) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_123122_Batch2
BeamDia(:,2) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_123122_Batch3
BeamDia(:,3) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_123122_Batch4
BeamDia(:,4) = BeamDiaSmallNA(:,2);
load BeamDiaSamllNA_123122_Batch5
BeamDia(:,5) = BeamDiaSmallNA(:,2);

VortexCharge = BeamDiaSmallNA(:,1);

BeamDiaMean = mean(BeamDia,2);
BeamDiaStd = std(BeamDia,[],2);

SaveFlag= 1;
h = figure;
h.Units ='inch';
h.Position=[2 2 2.42 2];
h.PaperPositionMode='auto';
ax_handle = axes(h, 'Position', [0.1672 0.1682 0.8 0.765]);

errorbar(VortexCharge([1:11,13:15]),BeamDiaMean([1:11,13:15]),BeamDiaStd([1:11,13:15]),'r-.','linewidth',1.2);
h_xlabel = xlabel('Vortex beam topological charge','Fontname', 'Calibri','FontSize', 10);
h_ylabel = ylabel('Beam diameter (\mum)','Fontname', 'Calibri','FontSize', 10);
set(h_xlabel,'Position',[13.0000   -1.5   -1.0000]);
set(h_ylabel,'Position',[-2.8   7.0000   -1.0000]);
grid on
axis([0 26 0 14]);
set(gca,'Fontname', 'Calibri','FontSize', 10);
if SaveFlag
   print(h,'-dtiffn','-r600','Fig4b_VortexBeamDiaMeanStd_Omit20.tif');
end

