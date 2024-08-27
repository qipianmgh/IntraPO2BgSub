clear
clc
close all

depth = 500;
vortexcharge = [0 1 5 11 15 20 25];

lambda = 0.95;  %%%% um
n = 1.33;
c = 3e14;  %%%%  um
w0 = 0.538; %%%%  um 0.538 was used for NA = 0.562
zR = n*pi*w0^2/lambda;
mus = 0e-3;   %%%% non-scattering medium  
g_aniso = 0.9;

x_limt = 12.5;
y_limt = 12.5;
spa_pts = 256;
x = linspace(-x_limt,x_limt,spa_pts);
y = linspace(-y_limt,y_limt,spa_pts);
[X Y] = meshgrid(x,y);
rho = sqrt(X.^2+Y.^2);
z_limit = 12;
stackdepth = linspace(depth-z_limit,depth+z_limit,spa_pts);
stacksize = length(stackdepth);

for m = length(depth)
    z0 = depth;
    z = stackdepth;
    kesi = 1+ i*(z-z0)/zR;
    
    PSFstack0 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack1 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack5 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack11 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack15 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack20 = zeros(spa_pts,spa_pts,stacksize); 
    PSFstack25 = zeros(spa_pts,spa_pts,stacksize);
    
    Fluostack0 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack1 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack5 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack11 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack15 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack20 = zeros(spa_pts,spa_pts,stacksize); 
    Fluostack25 = zeros(spa_pts,spa_pts,stacksize);  
    
    for k = 1:length(stackdepth)
        Temp0 = (2/pi/w0^2/abs(kesi(k))^2)*exp(-2*rho.^2/w0^2/abs(kesi(k))^2);
        Temp1 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(0,rho.^2/2/w0^2/kesi(k))-besseli(1,rho.^2/2/w0^2/kesi(k))).^2;
        Temp5 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(2,rho.^2/2/w0^2/kesi(k))-besseli(3,rho.^2/2/w0^2/kesi(k))).^2;
        Temp11 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(5,rho.^2/2/w0^2/kesi(k))-besseli(6,rho.^2/2/w0^2/kesi(k))).^2;
        Temp15 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(7,rho.^2/2/w0^2/kesi(k))-besseli(8,rho.^2/2/w0^2/kesi(k))).^2;
        Temp20 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(9.5,rho.^2/2/w0^2/kesi(k))-besseli(10.5,rho.^2/2/w0^2/kesi(k))).^2;
        Temp25 = (rho.^2/2/w0^4/abs(kesi(k))^3).*exp(-rho.^2/w0^2/abs(kesi(k))^2).*abs(besseli(12,rho.^2/2/w0^2/kesi(k))-besseli(13,rho.^2/2/w0^2/kesi(k))).^2;
        Temp1(isinf(Temp1))= 0;
        Temp1(isnan(Temp1))= 0;
        Temp5(isinf(Temp5))= 0;
        Temp5(isnan(Temp5))= 0;
        Temp11(isinf(Temp11))= 0;
        Temp11(isnan(Temp11))= 0;
        Temp15(isinf(Temp15))= 0;
        Temp15(isnan(Temp15))= 0;
        Temp20(isinf(Temp20))= 0;
        Temp20(isnan(Temp20))= 0;
        Temp25(isinf(Temp25))= 0;
        Temp25(isnan(Temp25))= 0;
        %%% Calculate two-photon excitation and emission at z(k)
        PSF_Ini0 = Temp0/sum(Temp0(:));
        PSF_Ini1 = Temp1/sum(Temp1(:));
        PSF_Ini5 = Temp5/sum(Temp5(:));
        PSF_Ini11 = Temp11/sum(Temp11(:));
        PSF_Ini15 = Temp15/sum(Temp15(:));
        PSF_Ini20 = Temp20/sum(Temp20(:));
        PSF_Ini25 = Temp25/sum(Temp25(:));
        
        PSFb0 = PSF_Ini0;
        PSFb1 = PSF_Ini1;
        PSFb5 = PSF_Ini5;
        PSFb11 = PSF_Ini11;
        PSFb15 = PSF_Ini15;
        PSFb20 = PSF_Ini20;
        PSFb25 = PSF_Ini25;

        Fluob0 = PSFb0.^2;
        Fluob1 = PSFb1.^2;
        Fluob5 = PSFb5.^2;
        Fluob11 = PSFb11.^2;
        Fluob15 = PSFb15.^2;
        Fluob20 = PSFb20.^2;
        Fluob25 = PSFb25.^2;
        
        PSFstack0(:,:,k) = PSFb0; 
        PSFstack1(:,:,k) = PSFb1; 
        PSFstack5(:,:,k) = PSFb5; 
        PSFstack11(:,:,k) = PSFb11; 
        PSFstack15(:,:,k) = PSFb15; 
        PSFstack20(:,:,k) = PSFb20; 
        PSFstack25(:,:,k) = PSFb25; 

        Fluostack0(:,:,k) = Fluob0; 
        Fluostack1(:,:,k) = Fluob1; 
        Fluostack5(:,:,k) = Fluob5; 
        Fluostack11(:,:,k) = Fluob11; 
        Fluostack15(:,:,k) = Fluob15; 
        Fluostack20(:,:,k) = Fluob20; 
        Fluostack25(:,:,k) = Fluob25; 
        m
        k
    end
    
    save(['Fluostack0_depth',num2str(z0),'um.mat'],'Fluostack0');
    save(['Fluostack1_depth',num2str(z0),'um.mat'],'Fluostack1');
    save(['Fluostack5_depth',num2str(z0),'um.mat'],'Fluostack5');
    save(['Fluostack11_depth',num2str(z0),'um.mat'],'Fluostack11');
    save(['Fluostack15_depth',num2str(z0),'um.mat'],'Fluostack15');
    save(['Fluostack20_depth',num2str(z0),'um.mat'],'Fluostack20');
    save(['Fluostack25_depth',num2str(z0),'um.mat'],'Fluostack25');
 
end
