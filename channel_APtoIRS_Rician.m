function [G] = channel_APtoIRS_Rician( Nx,Ny,Mx,My,UAV_location,IRS_location_center,lambda  )
%% 生成IRStoUser的LoS信道
%% 需输入横纵反射单元数Nx\Ny,横纵天线数Mx\My，反射单元中心坐标，无人机坐标
dis_AP_IRS = sqrt( sum( abs(UAV_location - IRS_location_center ).^2, 1 ) )';
%lambda = 0.1; % the wavelength of the carrier frequency. 0.1 m

M = Mx*My;
AAR_UR_1 = zeros(1,Mx);
AAR_UR_2 = zeros(1,My);
AAR_UR = zeros(1,M);
d_antenna = lambda/2;

N = Nx*Ny;
AAR_RU_1 = zeros(1,Nx);
AAR_RU_2 = zeros(1,Ny);
AAR_RU = zeros(1,N);
d_elements = lambda/2;

sin_vertical_AOD = abs(UAV_location(3)-IRS_location_center(3))/dis_AP_IRS;
cos_horizontal_AOD = abs(UAV_location(2)-IRS_location_center(2))/sqrt( (UAV_location(1)-IRS_location_center(1))^2+(UAV_location(2)-IRS_location_center(2))^2 );
sin_horizontal_AOD = abs(UAV_location(1)-IRS_location_center(1))/sqrt( (UAV_location(1)-IRS_location_center(1))^2+(UAV_location(2)-IRS_location_center(2))^2 );

racian_factor = 2; % 3 dB
sigma = 1;         % the variance of CSCG

for j = 1 : Mx
    arg = -2*pi*d_antenna*(j-1)/lambda*sin_vertical_AOD*cos_horizontal_AOD;
    AAR_UR_1(j) = exp(1i.*arg);
end
for k = 1 : My
    arg = -2*pi*d_antenna*(k-1)/lambda*sin_vertical_AOD*sin_horizontal_AOD;
    AAR_UR_2(k) = exp(1i.*arg);
end
AAR_UR = kron(AAR_UR_1',AAR_UR_2');

for j = 1 : Nx
    arg = -2*pi*d_elements*(j-1)/lambda*sin_vertical_AOD*cos_horizontal_AOD;
    AAR_RU_1(j) = exp(1i.*arg);
end
for k = 1 : Ny
    arg = -2*pi*d_elements*(k-1)/lambda*sin_vertical_AOD*sin_horizontal_AOD;
    AAR_RU_2(k) = exp(1i.*arg);
end
AAR_RU = kron(AAR_RU_1',AAR_RU_2');

G = kron(AAR_UR,AAR_RU').';
NLOS_average = 0;
for rand_num = 1 :200
    NLOS_average = NLOS_average + (randn(1,M)+1i*randn(1,M));
end
NLOS_average = NLOS_average/200;

G = sqrt(racian_factor/(racian_factor+1))*kron(AAR_UR,AAR_RU').' + sqrt(1/(racian_factor+1))*NLOS_average*sigma/sqrt(2)   ;
end

