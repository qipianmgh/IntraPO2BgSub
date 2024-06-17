clear
clc
close all

mus = 5e-3;   %%%% 5 mm^-1 = 5e-3 um^-1  
x_limit = 450;
y_limit = 450;
spa_pts = 4096;
z0 = 500;
z_ext = 100;

x_limit_disp = 300;
x_start_disp = round((x_limit-x_limit_disp)/x_limit*(spa_pts-1)/2); 
x_len_disp = round(x_limit_disp/x_limit*(spa_pts-1));

load(['FluomapProj0_depth',num2str(z0),'um.mat']);
load(['FluomapProj20_depth',num2str(z0),'um.mat']);
load(['FluomapSubProj20_depth',num2str(z0),'um.mat']);

NormFactor = max(FluomapProj0(:));
FluomapProj0 = FluomapProj0/NormFactor;
FluomapProj20 = FluomapProj20/NormFactor;
FluomapSubProj20 = FluomapSubProj20/NormFactor;

FluomapProj0 = FluomapProj0(:,x_start_disp+1:x_start_disp+x_len_disp);
FluomapProj20 = FluomapProj20(:,x_start_disp+1:x_start_disp+x_len_disp);
FluomapSubProj20 = FluomapSubProj20(:,x_start_disp+1:x_start_disp+x_len_disp);

h0 = figure;
h0.Units ='inch';
h0.Position=[2 2 1.60 1.7];
h0.PaperPositionMode='auto';
ax_handle = axes(h0, 'Position', [0.08 0.075 0.9 0.847]);

xscalebar_len = 100;
xscalebar_pts = 500;
xscalebar_x = linspace(150,150+xscalebar_len/2/x_limit_disp*x_len_disp,xscalebar_pts);
xscalebar_y = 450*ones(1,xscalebar_pts);

zscalebar_len = 100;
zscalebar_pts = 500;
zscalebar_x = 150*ones(1,zscalebar_pts);
zscalebar_y = linspace(450,450+zscalebar_len/(z0+z_ext-0)*(size(FluomapProj0,1)-1),zscalebar_pts);

imagesc(log(FluomapProj0+eps),[-14 0]);
colormap hot;
hold on
plot(xscalebar_x,xscalebar_y,'w-','LineWidth',2);
plot(zscalebar_x,zscalebar_y,'w-','LineWidth',2);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
h_xlabel = xlabel('x','Fontname', 'Calibri','FontWeight','bold','FontSize',12,'Position',[(x_len_disp+1)/2 580 1]);
h_ylabel = ylabel('z','Fontname', 'Calibri','FontWeight','bold','FontSize',12,'Position',[80 301 1]);
print(h0,'-dtiffn','-r1200','F0_Vortex20_v1.tif');

h1 = figure;
h1.Units ='inch';
h1.Position=[2 2 1.82 1.70];
h1.PaperPositionMode='auto';
ax_handle = axes(h1, 'Position', [0.015 0.075 0.791 0.847]);

imagesc(log(FluomapProj20+eps),[-14 0]);
colormap hot;
h_bar = colorbar;
set(h_bar,'Ticks',[-14 0],'TickLabels',[ ],'Fontname','Calibri','FontWeight','bold','FontSize',10,'Position',[0.815 0.075 0.04 0.847]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
text(x_len_disp+210,10,['0'],'Color','k','Fontname', 'Calibri','FontWeight','bold','FontSize',8.5);
text(x_len_disp+210,600,['-14'],'Color','k','Fontname', 'Calibri','FontWeight','bold','FontSize',8.5);
h_xlabel = xlabel('x','Fontname', 'Calibri','FontWeight','bold','FontSize',12,'Position',[(x_len_disp+1)/2 580 1]);
print(h1,'-dtiffn','-r1200','FB_Vortex20_v1.tif');

h2 = figure;
h2.Units ='inch';
h2.Position=[2 2 1.82 1.70];
h2.PaperPositionMode='auto';
ax_handle = axes(h2, 'Position', [0.015 0.075 0.791 0.847]);

imagesc(FluomapSubProj20,[-1e-5 1e-5]);
cm1a = [32:-1:1]'*[0 1 0]/32; 
cm1a(:,3) = 1;
cm1b = [32:-1:1]'*[0 0 1]/32;
cm2 = hot(64);
colormap([cm1a; cm1b; cm2]);
h_bar = colorbar;
set(h_bar,'Ticks',[-1e-5 1e-5],'TickLabels',[ ],'Fontname','Calibri','FontWeight','bold','FontSize',10,'Position',[0.815 0.075 0.04 0.847]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
text(x_len_disp+200,0,['10^{-5}'],'Color','k','Fontname', 'Calibri','FontWeight','bold','FontSize',8.5);
text(x_len_disp+200,590,['-10^{-5}'],'Color','k','Fontname', 'Calibri','FontWeight','bold','FontSize',8.5);
h_xlabel = xlabel('x','Fontname', 'Calibri','FontWeight','bold','FontSize',12,'Position',[(x_len_disp+1)/2 580 1]);
print(h2,'-dtiffn','-r1200','F0_Sub_Vortex20_v1.tif');
