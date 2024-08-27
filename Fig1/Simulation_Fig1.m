clear
clc
close all

depth = [500];
vortexcharge = [0 1 5 10 15 20 25];

lambda = 0.95;  %%%% um
n = 1.33;
c = 3e14;  %%%%  um
w0 = 0.538; %%%%  um 0.538 was used for NA = 0.562
zR = n*pi*w0^2/lambda;
mus = 5e-3;   %%%% 5 mm^-1 = 5e-3 um^-1  
z_ext = round(0.5/mus);
g_aniso = 0.9;

x_limit = 450;
y_limit = 450;
spa_pts = 4096;
x = linspace(-x_limit,x_limit,spa_pts);
y = linspace(-y_limit,y_limit,spa_pts);
[X Y] = meshgrid(x,y);
rho = sqrt(X.^2+Y.^2);
int_scalefactor = 2*x_limit/(spa_pts-1)*2*y_limit/(spa_pts-1);  %%% scaling factor for 2D integration
int_scalefactor_1d = 2*y_limit/(spa_pts-1);  %%% scaling factor for integration in y dimension

for m = 1:length(depth)
    z0 = depth(m);
    w0_surface = w0*abs(1+i*z0/zR);
    z = linspace(0,z0+z_ext,z0+z_ext+1);
    kesi = 1+ i*(z-z0)/zR;
    
    PSFbSum0 = zeros(1,length(z));
    PSFbSum1 = zeros(1,length(z));
    PSFbSum5 = zeros(1,length(z));
    PSFbSum10 = zeros(1,length(z));
    PSFbSum15 = zeros(1,length(z));
    PSFbSum20 = zeros(1,length(z));
    PSFbSum25 = zeros(1,length(z));
    
    PSFsSum0 = zeros(1,length(z));
    PSFsSum1 = zeros(1,length(z));
    PSFsSum5 = zeros(1,length(z));
    PSFsSum10 = zeros(1,length(z));
    PSFsSum15 = zeros(1,length(z));
    PSFsSum20 = zeros(1,length(z));
    PSFsSum25 = zeros(1,length(z));

    FluobSum0 = zeros(1,length(z));
    FluobSum1 = zeros(1,length(z));
    FluobSum5 = zeros(1,length(z));
    FluobSum10 = zeros(1,length(z));
    FluobSum15 = zeros(1,length(z));
    FluobSum20 = zeros(1,length(z));
    FluobSum25 = zeros(1,length(z));

    FluosSum0 = zeros(1,length(z));
    FluosSum1 = zeros(1,length(z));
    FluosSum5 = zeros(1,length(z));
    FluosSum10 = zeros(1,length(z));
    FluosSum15 = zeros(1,length(z));
    FluosSum20 = zeros(1,length(z));
    FluosSum25 = zeros(1,length(z));

    FluocSum0 = zeros(1,length(z));
    FluocSum1 = zeros(1,length(z));
    FluocSum5 = zeros(1,length(z));
    FluocSum10 = zeros(1,length(z));
    FluocSum15 = zeros(1,length(z));
    FluocSum20 = zeros(1,length(z));
    FluocSum25 = zeros(1,length(z));
    
    FluoSum0_NewMethod = zeros(1,length(z));
    FluoSum1_NewMethod = zeros(1,length(z));
    FluoSum5_NewMethod = zeros(1,length(z));
    FluoSum10_NewMethod = zeros(1,length(z));
    FluoSum15_NewMethod = zeros(1,length(z));
    FluoSum20_NewMethod = zeros(1,length(z));
    FluoSum25_NewMethod = zeros(1,length(z));
       
    FluomapProj0 = zeros(length(z),spa_pts); 
    FluomapProj1 = zeros(length(z),spa_pts); 
    FluomapProj5 = zeros(length(z),spa_pts); 
    FluomapProj10 = zeros(length(z),spa_pts); 
    FluomapProj15 = zeros(length(z),spa_pts); 
    FluomapProj20 = zeros(length(z),spa_pts); 
    FluomapProj25 = zeros(length(z),spa_pts); 
        
    FluomapSubProj1 = zeros(length(z),spa_pts); 
    FluomapSubProj5 = zeros(length(z),spa_pts); 
    FluomapSubProj10 = zeros(length(z),spa_pts); 
    FluomapSubProj15 = zeros(length(z),spa_pts); 
    FluomapSubProj20 = zeros(length(z),spa_pts); 
    FluomapSubProj25 = zeros(length(z),spa_pts); 
    
    for k = 1:length(z)
        Fluomap0 = zeros(spa_pts,spa_pts); 
        Fluomap1 = zeros(spa_pts,spa_pts); 
        Fluomap5 = zeros(spa_pts,spa_pts); 
        Fluomap10 = zeros(spa_pts,spa_pts); 
        Fluomap15 = zeros(spa_pts,spa_pts); 
        Fluomap20 = zeros(spa_pts,spa_pts); 
        Fluomap25 = zeros(spa_pts,spa_pts); 
        
        FluomapSub1 = zeros(spa_pts,spa_pts); 
        FluomapSub5 = zeros(spa_pts,spa_pts); 
        FluomapSub10 = zeros(spa_pts,spa_pts); 
        FluomapSub15 = zeros(spa_pts,spa_pts); 
        FluomapSub20 = zeros(spa_pts,spa_pts); 
        FluomapSub25 = zeros(spa_pts,spa_pts); 
        
        %%% Generate ballistic light intensity at z(k)
        Temp0 = (2/pi/w0^2/abs(kesi(k))^2)*exp(-2*rho.^2/w0^2/abs(kesi(k))^2);
        Temp1 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(0,rho.^2/2/w0^2/kesi(k))-besseli(1,rho.^2/2/w0^2/kesi(k))).^2;
        Temp5 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(2,rho.^2/2/w0^2/kesi(k))-besseli(3,rho.^2/2/w0^2/kesi(k))).^2;
        Temp10 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(4.5,rho.^2/2/w0^2/kesi(k))-besseli(5.5,rho.^2/2/w0^2/kesi(k))).^2;
        Temp15 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(7,rho.^2/2/w0^2/kesi(k))-besseli(8,rho.^2/2/w0^2/kesi(k))).^2;
        Temp20 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(9.5,rho.^2/2/w0^2/kesi(k))-besseli(10.5,rho.^2/2/w0^2/kesi(k))).^2;
        Temp25 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(12,rho.^2/2/w0^2/kesi(k))-besseli(13,rho.^2/2/w0^2/kesi(k))).^2;
        Temp1(isinf(Temp1))= 0;
        Temp1(isnan(Temp1))= 0;
        Temp5(isinf(Temp5))= 0;
        Temp5(isnan(Temp5))= 0;
        Temp10(isinf(Temp10))= 0;
        Temp10(isnan(Temp10))= 0;
        Temp15(isinf(Temp15))= 0;
        Temp15(isnan(Temp15))= 0;
        Temp20(isinf(Temp20))= 0;
        Temp20(isnan(Temp20))= 0;
        Temp25(isinf(Temp25))= 0;
        Temp25(isnan(Temp25))= 0;
        %%% Calculate total fluorescence/phorsopherescence emission at z(k)
        PSF_Ini0 = Temp0;
        PSF_Ini1 = Temp1;
        PSF_Ini5 = Temp5;
        PSF_Ini10 = Temp10;
        PSF_Ini15 = Temp15;
        PSF_Ini20 = Temp20;
        PSF_Ini25 = Temp25;
        
        PSFb0 = PSF_Ini0*exp(-mus*z(k));
        PSFb1 = PSF_Ini1*exp(-mus*z(k));
        PSFb5 = PSF_Ini5*exp(-mus*z(k));
        PSFb10 = PSF_Ini10*exp(-mus*z(k));
        PSFb15 = PSF_Ini15*exp(-mus*z(k));
        PSFb20 = PSF_Ini20*exp(-mus*z(k));
        PSFb25 = PSF_Ini25*exp(-mus*z(k));

        Fluob0 = PSFb0.^2;
        Fluob1 = PSFb1.^2;
        Fluob5 = PSFb5.^2;
        Fluob10 = PSFb10.^2;
        Fluob15 = PSFb15.^2;
        Fluob20 = PSFb20.^2;
        Fluob25 = PSFb25.^2;
        
        PSFbSum0(k) = sum(PSFb0(:))*int_scalefactor;
        PSFbSum1(k) = sum(PSFb1(:))*int_scalefactor;
        PSFbSum5(k) = sum(PSFb5(:))*int_scalefactor;
        PSFbSum10(k) = sum(PSFb10(:))*int_scalefactor;
        PSFbSum15(k) = sum(PSFb15(:))*int_scalefactor;
        PSFbSum20(k) = sum(PSFb20(:))*int_scalefactor;
        PSFbSum25(k) = sum(PSFb25(:))*int_scalefactor;

        FluobSum0(k) = sum(Fluob0(:))*int_scalefactor;
        FluobSum1(k) = sum(Fluob1(:))*int_scalefactor;
        FluobSum5(k) = sum(Fluob5(:))*int_scalefactor;
        FluobSum10(k) = sum(Fluob10(:))*int_scalefactor;
        FluobSum15(k) = sum(Fluob15(:))*int_scalefactor;
        FluobSum20(k) = sum(Fluob20(:))*int_scalefactor;
        FluobSum25(k) = sum(Fluob25(:))*int_scalefactor;
   
        if z(k)==0   
           PSFs0 = zeros(spa_pts,spa_pts);
           PSFs1 = zeros(spa_pts,spa_pts);
           PSFs5 = zeros(spa_pts,spa_pts);
           PSFs10 = zeros(spa_pts,spa_pts);
           PSFs15 = zeros(spa_pts,spa_pts);
           PSFs20 = zeros(spa_pts,spa_pts);
           PSFs25 = zeros(spa_pts,spa_pts);  
        else
           beta = 3*n/mus/2/(1-g_aniso)/(z(k)^3);
           BSF = beta/pi*exp(-beta*rho.^2);
           BSF = BSF/sum(BSF(:));

           PSFs0 = conv2fft(PSF_Ini0*(1-exp(-mus*z(k))),BSF,'same');
           PSFs1 = conv2fft(PSF_Ini1*(1-exp(-mus*z(k))),BSF,'same');
           PSFs5 = conv2fft(PSF_Ini5*(1-exp(-mus*z(k))),BSF,'same');
           PSFs10 = conv2fft(PSF_Ini10*(1-exp(-mus*z(k))),BSF,'same');
           PSFs15 = conv2fft(PSF_Ini15*(1-exp(-mus*z(k))),BSF,'same');
           PSFs20 = conv2fft(PSF_Ini20*(1-exp(-mus*z(k))),BSF,'same');
           PSFs25 = conv2fft(PSF_Ini25*(1-exp(-mus*z(k))),BSF,'same');            
        end

        Fluos0 = PSFs0.^2;
        Fluos1 = PSFs1.^2;
        Fluos5 = PSFs5.^2;
        Fluos10 = PSFs10.^2;
        Fluos15 = PSFs15.^2;
        Fluos20 = PSFs20.^2;
        Fluos25 = PSFs25.^2;

        Fluoc0 = 2*PSFb0.*PSFs0;
        Fluoc1 = 2*PSFb1.*PSFs1;
        Fluoc5 = 2*PSFb5.*PSFs5;       
        Fluoc10 = 2*PSFb10.*PSFs10;       
        Fluoc15 = 2*PSFb15.*PSFs15;       
        Fluoc20 = 2*PSFb20.*PSFs20; 
        Fluoc25 = 2*PSFb25.*PSFs25; 
        
        Fluomap0 = Fluob0 + Fluos0 + Fluoc0; 
        Fluomap1 = Fluob1 + Fluos1 + Fluoc1; 
        Fluomap5 = Fluob5 + Fluos5 + Fluoc5; 
        Fluomap10 = Fluob10 + Fluos10 + Fluoc10; 
        Fluomap15 = Fluob15 + Fluos15 + Fluoc15; 
        Fluomap20 = Fluob20 + Fluos20 + Fluoc20; 
        Fluomap25 = Fluob25 + Fluos25 + Fluoc25; 
        
        PSFsSum0(k) = sum(PSFs0(:))*int_scalefactor;
        PSFsSum1(k) = sum(PSFs1(:))*int_scalefactor;
        PSFsSum5(k) = sum(PSFs5(:))*int_scalefactor;
        PSFsSum10(k) = sum(PSFs10(:))*int_scalefactor;
        PSFsSum15(k) = sum(PSFs15(:))*int_scalefactor;
        PSFsSum20(k) = sum(PSFs20(:))*int_scalefactor;
        PSFsSum25(k) = sum(PSFs25(:))*int_scalefactor;
        
        FluosSum0(k) = sum(Fluos0(:))*int_scalefactor;
        FluosSum1(k) = sum(Fluos1(:))*int_scalefactor;
        FluosSum5(k) = sum(Fluos5(:))*int_scalefactor;
        FluosSum10(k) = sum(Fluos10(:))*int_scalefactor;
        FluosSum15(k) = sum(Fluos15(:))*int_scalefactor;
        FluosSum20(k) = sum(Fluos20(:))*int_scalefactor;
        FluosSum25(k) = sum(Fluos25(:))*int_scalefactor;

        FluocSum0(k) = sum(Fluoc0(:))*int_scalefactor;
        FluocSum1(k) = sum(Fluoc1(:))*int_scalefactor;
        FluocSum5(k) = sum(Fluoc5(:))*int_scalefactor;
        FluocSum10(k) = sum(Fluoc10(:))*int_scalefactor;
        FluocSum15(k) = sum(Fluoc15(:))*int_scalefactor;
        FluocSum20(k) = sum(Fluoc20(:))*int_scalefactor;
        FluocSum25(k) = sum(Fluoc25(:))*int_scalefactor;

        FluoSum0_NewMethod(k) = sum(Fluomap0(:))*int_scalefactor;
        FluoSum1_NewMethod(k) = sum(Fluomap1(:))*int_scalefactor;
        FluoSum5_NewMethod(k) = sum(Fluomap5(:))*int_scalefactor;
        FluoSum10_NewMethod(k) = sum(Fluomap10(:))*int_scalefactor;
        FluoSum15_NewMethod(k) = sum(Fluomap15(:))*int_scalefactor;
        FluoSum20_NewMethod(k) = sum(Fluomap20(:))*int_scalefactor;
        FluoSum25_NewMethod(k) = sum(Fluomap25(:))*int_scalefactor;
        
        FluomapProj0(k,:) = sum(Fluomap0,1)*int_scalefactor_1d;
        FluomapProj1(k,:) = sum(Fluomap1,1)*int_scalefactor_1d; 
        FluomapProj5(k,:) = sum(Fluomap5,1)*int_scalefactor_1d; 
        FluomapProj10(k,:) = sum(Fluomap10,1)*int_scalefactor_1d; 
        FluomapProj15(k,:) = sum(Fluomap15,1)*int_scalefactor_1d; 
        FluomapProj20(k,:) = sum(Fluomap20,1)*int_scalefactor_1d; 
        FluomapProj25(k,:) = sum(Fluomap25,1)*int_scalefactor_1d; 

        FluomapSub1 = Fluomap0-Fluomap1; 
        FluomapSub5 = Fluomap0-Fluomap5; 
        FluomapSub10 = Fluomap0-Fluomap10; 
        FluomapSub15 = Fluomap0-Fluomap15; 
        FluomapSub20 = Fluomap0-Fluomap20; 
        FluomapSub25 = Fluomap0-Fluomap25; 

        FluomapSubProj1(k,:) = sum(FluomapSub1,1)*int_scalefactor_1d;
        FluomapSubProj5(k,:) = sum(FluomapSub5,1)*int_scalefactor_1d;
        FluomapSubProj10(k,:) = sum(FluomapSub10,1)*int_scalefactor_1d;
        FluomapSubProj15(k,:) = sum(FluomapSub15,1)*int_scalefactor_1d;
        FluomapSubProj20(k,:) = sum(FluomapSub20,1)*int_scalefactor_1d;
        FluomapSubProj25(k,:) = sum(FluomapSub25,1)*int_scalefactor_1d;
        m
        k
    end
       
    FluoSum0 = FluobSum0+FluosSum0+FluocSum0;
    FluoSum1 = FluobSum1+FluosSum1+FluocSum1;
    FluoSum5 = FluobSum5+FluosSum5+FluocSum5;
    FluoSum10 = FluobSum10+FluosSum10+FluocSum10;
    FluoSum15 = FluobSum15+FluosSum15+FluocSum15;
    FluoSum20 = FluobSum20+FluosSum20+FluocSum20;
    FluoSum25 = FluobSum25+FluosSum25+FluocSum25;
    
    save(['FluobSum_',num2str(depth(m)),'um.mat'],'FluobSum0','FluobSum1','FluobSum5','FluobSum10','FluobSum15','FluobSum20','FluobSum25');
    save(['FluosSum_',num2str(depth(m)),'um.mat'],'FluosSum0','FluosSum1','FluosSum5','FluosSum10','FluosSum15','FluosSum20','FluosSum25');
    save(['FluocSum_',num2str(depth(m)),'um.mat'],'FluocSum0','FluocSum1','FluocSum5','FluocSum10','FluocSum15','FluocSum20','FluocSum25');
    save(['FluoSum_',num2str(depth(m)),'um.mat'],'FluoSum0','FluoSum1','FluoSum5','FluoSum10','FluoSum15','FluoSum20','FluoSum25');
    save(['PSFbSum_',num2str(depth(m)),'um.mat'],'PSFbSum0','PSFbSum1','PSFbSum5','PSFbSum10','PSFbSum15','PSFbSum20','PSFbSum25');
    save(['PSFsSum_',num2str(depth(m)),'um.mat'],'PSFsSum0','PSFsSum1','PSFsSum5','PSFsSum10','PSFsSum15','PSFsSum20','PSFsSum25');

    save(['FluoSum_NewMethod_',num2str(depth(m)),'um.mat'],'FluoSum0_NewMethod','FluoSum1_NewMethod','FluoSum5_NewMethod','FluoSum10_NewMethod','FluoSum15_NewMethod','FluoSum20_NewMethod','FluoSum25_NewMethod');
    
    save(['FluomapProj0_depth',num2str(z0),'um.mat'],'FluomapProj0');
    save(['FluomapProj1_depth',num2str(z0),'um.mat'],'FluomapProj1');
    save(['FluomapProj5_depth',num2str(z0),'um.mat'],'FluomapProj5');
    save(['FluomapProj10_depth',num2str(z0),'um.mat'],'FluomapProj10');
    save(['FluomapProj15_depth',num2str(z0),'um.mat'],'FluomapProj15');
    save(['FluomapProj20_depth',num2str(z0),'um.mat'],'FluomapProj20');
    save(['FluomapProj25_depth',num2str(z0),'um.mat'],'FluomapProj25');

    save(['FluomapSubProj1_depth',num2str(z0),'um.mat'],'FluomapSubProj1');
    save(['FluomapSubProj5_depth',num2str(z0),'um.mat'],'FluomapSubProj5');
    save(['FluomapSubProj10_depth',num2str(z0),'um.mat'],'FluomapSubProj10');
    save(['FluomapSubProj15_depth',num2str(z0),'um.mat'],'FluomapSubProj15');
    save(['FluomapSubProj20_depth',num2str(z0),'um.mat'],'FluomapSubProj20');
    save(['FluomapSubProj25_depth',num2str(z0),'um.mat'],'FluomapSubProj25');
end
