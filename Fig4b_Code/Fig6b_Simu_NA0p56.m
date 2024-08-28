clear
clc
close all

%%% Medium parameters
lambda = 0.95;  %%%% um
n = 1.33;
c = 3e14;  %%%%  um
w0 = 0.538; %%%%  um 0.538 was used for NA = 0.562
zR = n*pi*w0^2/lambda;

%%% Spatial coordinates
x = linspace(-30,30,601);
y = linspace(-30,30,601);
[X Y] = meshgrid(x,y);
rho = sqrt(X.^2+Y.^2);
z = 0;
kesi = 1+ i*z/zR;
FOV = max(x)-min(x);   %%%  um
vortexcharge = [0 1:2:19 20 21:2:25];
BeamDist_Emi = zeros(length(x),length(y),length(vortexcharge));
BeamDia_Emi = zeros(length(vortexcharge),1);

%%% Simulation
for k = 1:length(vortexcharge)
    if vortexcharge(k) == 0 
       Temp = (2/pi/w0^2/abs(kesi)^2)*exp(-2*rho.^2/w0^2/abs(kesi)^2); 
       BeamDist_Emi(:,:,k) = Temp.^2;

       I_emi = BeamDist_Emi(:,:,k); 
       I_emi_norm_max = I_emi/max(I_emi(:));
       I_emi_mask = I_emi_norm_max >= 0.1353;
       props_emi = regionprops(I_emi_mask, 'Area', 'EquivDiameter', 'Centroid');
       dia_pixel_emi = props_emi.EquivDiameter;
       BeamDia_Emi(k) = dia_pixel_emi*FOV/(length(x)-1);
    else
       Temp = (rho.^2/2/w0^4/abs(kesi)^3).*exp(-rho.^2/w0^2/abs(kesi)^2).*abs(besseli((vortexcharge(k)-1)/2,rho.^2/2/w0^2/kesi)-besseli((vortexcharge(k)+1)/2,rho.^2/2/w0^2/kesi)).^2;
       Temp(isinf(Temp))= 0;
       Temp(isnan(Temp))= 0;
       BeamDist_Emi(:,:,k) = Temp.^2;
       
       I_emi = BeamDist_Emi(:,:,k); 
       I_emi_norm_max = I_emi/max(I_emi(:));
       I_emi_mask = I_emi_norm_max >= 0.01;
       I_emi_mask = bwareafilt(I_emi_mask, 1, 8);
       I_emi = I_emi.*I_emi_mask;
       I_emi_norm_sum = I_emi/sum(I_emi(:));
       [row,col] = size(I_emi_norm_sum);
       [XX,YY] = ndgrid(1:row,1:col);
       row_centroid = sum(sum(XX.*I_emi_norm_sum));
       col_centroid = sum(sum(YY.*I_emi_norm_sum));
       ind = find(I_emi_mask == 1);
       [row_beam col_beam]= ind2sub(size(I_emi_mask),ind);
       BeamDia_Emi(k) = 2*FOV/(length(x)-1)*sum(I_emi_norm_sum(ind).*sqrt((row_beam-row_centroid).^2+(col_beam-col_centroid).^2));
    end
end

BeamDiaSimu_Emi = [vortexcharge' BeamDia_Emi];
save('BeamDiaSimu_Emi_NA_0p56','BeamDiaSimu_Emi');