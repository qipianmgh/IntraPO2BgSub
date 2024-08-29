clear
clc
close all

depth = [50 100 200 300 400 500 600];
vortexcharge = [0 1 5 10 15 20 25];

lambda = 0.95;  %%%% um
n = 1.33;
c = 3e14;  %%%%  um
w0 = 0.538; %%%%  um 0.538 was used for NA = 0.562
zR = n*pi*w0^2/lambda;
mus = 5e-3;   %%%% 5 mm^-1 = 5e-3 um^-1  
g_aniso = 0.9;
z_ext = round(0.5/mus);
z_ext_collect = round(0.5/mus);
z_density = 10;  %%% z points per um
int_scalefactor_z = 1/z_density;  %%% scaling factor for integration
z_focaldepth = round(10*z_density);  %%% ~20um depth of field +-10um

x_limt = 450;
y_limt = 450;
spa_pts = 4095;
x = linspace(-x_limt,x_limt,spa_pts);
y = linspace(-y_limt,y_limt,spa_pts);
[X Y] = meshgrid(x,y);
rho = sqrt(X.^2+Y.^2);

depth_plot = [50 100 200 400 600];
for m=1:length(depth)
load(['FluobSum_',num2str(depth(m)),'um.mat']);    
load(['FluoSum_',num2str(depth(m)),'um.mat']);
z0 = depth(m);
w0_surface = w0*abs(1+i*z0/zR);
z = linspace(0,z0+z_ext,z0+z_ext+1);
z_interp = linspace(0,z0+z_ext,(z0+z_ext)*z_density+1);

FluoSum0_interp = interp1(z,FluoSum0,z_interp,'spline');
FluoSum1_interp = interp1(z,FluoSum1,z_interp,'spline');
FluoSum5_interp = interp1(z,FluoSum5,z_interp,'spline');
FluoSum10_interp = interp1(z,FluoSum10,z_interp,'spline');
FluoSum15_interp = interp1(z,FluoSum15,z_interp,'spline');
FluoSum20_interp = interp1(z,FluoSum20,z_interp,'spline');
FluoSum25_interp = interp1(z,FluoSum25,z_interp,'spline');

FluobSum0_interp = interp1(z,FluobSum0,z_interp,'spline');
FluobSum1_interp = interp1(z,FluobSum1,z_interp,'spline');
FluobSum5_interp = interp1(z,FluobSum5,z_interp,'spline');
FluobSum10_interp = interp1(z,FluobSum10,z_interp,'spline');
FluobSum15_interp = interp1(z,FluobSum15,z_interp,'spline');
FluobSum20_interp = interp1(z,FluobSum20,z_interp,'spline');
FluobSum25_interp = interp1(z,FluobSum25,z_interp,'spline');

Temp_z = abs(z_interp-(z0+z_ext_collect));
zmax_Fluo_index = find(Temp_z == min(Temp_z(:)));
FluoTotal0 = sum(FluoSum0_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal1 = sum(FluoSum1_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal5 = sum(FluoSum5_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal10 = sum(FluoSum10_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal15 = sum(FluoSum15_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal20 = sum(FluoSum20_interp(1:zmax_Fluo_index))*int_scalefactor_z;
FluoTotal25 = sum(FluoSum25_interp(1:zmax_Fluo_index))*int_scalefactor_z;

Temp2_z = abs(z_interp-z0);
zoof_Fluo_index = find(Temp2_z == min(Temp2_z(:)));
FluoOofTotal0 = FluoTotal0 - sum(FluobSum0_interp(zoof_Fluo_index-z_focaldepth:zoof_Fluo_index+z_focaldepth))*int_scalefactor_z;

FluoPerifTotal0 = sum(FluobSum0_interp(zoof_Fluo_index-z_focaldepth:zoof_Fluo_index+z_focaldepth))*int_scalefactor_z;

FluoOofRatio(m,:) = [FluoTotal0 FluoOofTotal0 FluoTotal1 FluoTotal5 FluoTotal10 FluoTotal15 FluoTotal20 FluoTotal25]/FluoTotal0;
FluoOofRatioPercent(m,:) = [FluoTotal0 FluoOofTotal0 FluoTotal1 FluoTotal5 FluoTotal10 FluoTotal15 FluoTotal20 FluoTotal25]/FluoTotal0*100;

FluoPfRatio(m,1) = FluoPerifTotal0/FluoTotal0;
FluoPfRatioPercent(m,1) = FluoPerifTotal0/FluoTotal0*100;

FluoOofTotalEst(m,:) = [FluoOofTotal0 FluoTotal1 FluoTotal5 FluoTotal10 FluoTotal15 FluoTotal20 FluoTotal25];
FluoOofTotalError(m,:) = [FluoTotal1 FluoTotal5 FluoTotal10 FluoTotal15 FluoTotal20 FluoTotal25]-FluoOofTotal0;

RoofError(m,:) = abs([FluoTotal1 FluoTotal5 FluoTotal10 FluoTotal15 FluoTotal20 FluoTotal25]-FluoOofTotal0)/FluoTotal0;

if any(depth_plot == depth(m))
   mm = find(depth_plot == depth(m));
   h = figure;
   linewidth = 2;
   fontsize = 5;
   h.Units ='inch';
   h.Position=[2 2 4 3];
   h.PaperPositionMode='auto';
   ax_handle = axes(h, 'Position', [0.2 0.18 0.755 0.8]);
   semilogy(z_interp,FluoSum0_interp/max(FluoSum0_interp),'r','Linewidth',linewidth);
   hold on
   semilogy(z_interp,FluoSum1_interp/max(FluoSum0_interp),'g','Linewidth',linewidth);
   semilogy(z_interp,FluoSum5_interp/max(FluoSum0_interp),'b','Linewidth',linewidth);
   semilogy(z_interp,FluoSum10_interp/max(FluoSum0_interp),'c','Linewidth',linewidth);
   semilogy(z_interp,FluoSum15_interp/max(FluoSum0_interp),'m','Linewidth',linewidth);
   semilogy(z_interp,FluoSum20_interp/max(FluoSum0_interp),'k','Linewidth',linewidth);
   semilogy(z_interp,FluoSum25_interp/max(FluoSum0_interp),'y','Linewidth',linewidth);
   if depth(m) == 50
      xlim([0 z0+z_ext+60]);
   else
      xlim([-150 z0+z_ext]);
   end
   ylim([0 0.65]);
   yticks([0.001 0.01 0.1]);
   yticklabels({'10^{-3}','10^{-2}','10^{-1}'})
   h_xlabel = xlabel('z (\mum)','Fontname', 'Calibri','Fontsize',16);
   if depth(m) == 50
      set(h_xlabel,'Position', [(z0+z_ext+60)/2 8e-5 -1]);
   elseif depth(m) == 100
      set(h_xlabel,'Position', [z0+z_ext-(z0+z_ext+150)/2 9e-5 -1]); 
   elseif depth(m) == 200
      set(h_xlabel,'Position', [z0+z_ext-(z0+z_ext+150)/2 1e-4 -1]); 
   elseif depth(m) == 400
      set(h_xlabel,'Position', [z0+z_ext-(z0+z_ext+150)/2 1.6e-4 -1]); 
   else
      set(h_xlabel,'Position', [z0+z_ext-(z0+z_ext+150)/2 4.3e-4 -1]); 
   end
   h_ylabel = ylabel('TPEF','Fontname', 'Calibri','Fontsize',16);
   set(gca,'Fontname', 'Calibri','FontSize', 16);
   grid on
   h_legend = legend('Gaussian beam','Vortex charge = 1','Vortex charge = 5',...
                  'Vortex charge = 10','Vortex charge = 15','Vortex charge = 20','Vortex charge = 25');
   if depth(m)< z_ext
      set(h_legend,'Position',[0.53    0.45    0.40    0.52],'Fontname', 'Calibri','Fontsize',11);
   else
      set(h_legend,'Position',[0.225    0.45    0.40    0.52],'Fontname', 'Calibri','Fontsize',11);
   end
   print(h,['LogCompFluoTotal_',num2str(depth(m)),'um'],'-dmeta');  
end

end

save('FluoOofRatio.mat','FluoOofRatio');
save('FluoPfRatio.mat','FluoPfRatio');
save('FluoOofTotalEst','FluoOofTotalEst');
save('FluoOofTotalError','FluoOofTotalError');
save('RoofError','RoofError');

hh = figure;
linewidth = 2;
fontsize = 15;
hh.Units ='inch';
hh.Position=[2 2 4 3];
hh.PaperPositionMode='auto';
ax_handle = axes(hh, 'Position', [0.23 0.18 0.755 0.805]);
plot(depth',RoofError(:,1),'ro','MarkerSize',6);
hold on
plot(depth',RoofError(:,2),'g>','MarkerSize',6);
plot(depth',RoofError(:,3),'bp','MarkerSize',6);
plot(depth',RoofError(:,4),'cd','MarkerSize',6);
plot(depth',RoofError(:,5),'m^','MarkerSize',6);
plot(depth',RoofError(:,6),'ks','MarkerSize',6);
xlim([0 max(depth)+50]);
ylim([0 0.62]);
set(gca,'Fontname', 'Calibri','FontSize', 16);
h_xlabel = xlabel('Focal plane depth (\mum)','Fontname', 'Calibri','Fontsize',16);
set(h_xlabel,'Position', [325 -0.075 -1]);
h_ylabel = ylabel({'Absolute error of ';'R_{oof} estimation'},'Fontname', 'Calibri','Fontsize',16);
set(h_ylabel,'Position', [-65  0.31 -1]);
grid on
% Plot legend manually
plot(50,0.58,'ro','MarkerSize',6);
plot(360,0.58,'g>','MarkerSize',6);
plot(50,0.53,'bp','MarkerSize',6);
plot(360,0.53,'cd','MarkerSize',6);
plot(50,0.48,'m^','MarkerSize',6);
plot(360,0.48,'ks','MarkerSize',6);
text(70,0.58,'Vortex charge = 1','Fontname', 'Calibri','Fontsize',11);
text(380,0.58,'Vortex charge = 5','Fontname', 'Calibri','Fontsize',11);
text(70,0.53,'Vortex charge = 10','Fontname', 'Calibri','Fontsize',11);
text(380,0.53,'Vortex charge = 15','Fontname', 'Calibri','Fontsize',11);
text(70,0.48,'Vortex charge = 20','Fontname', 'Calibri','Fontsize',11);
text(380,0.48,'Vortex charge = 25','Fontname', 'Calibri','Fontsize',11);
rectangle('Position',[30 0.45 610 0.162]);
print(hh,'LogCompFluoTotal_BgError_v2','-dmeta');  