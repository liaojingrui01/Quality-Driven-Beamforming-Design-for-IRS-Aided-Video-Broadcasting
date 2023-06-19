function [Hi_d] = channel_UAVtoUser( Mx,My,Ki,UAV_location,IU_location,lambda  )
%% 生成UAVtoUser的Racia信道
%% 需输入横纵天线数Mx\My,用户数Ki，无人机坐标，用户坐标
dis_IU_AP = sqrt( sum( abs(UAV_location - IU_location ).^2, 1 ) )';

M = Mx*My;
AAR_1 = zeros(Ki,Mx);
AAR_2 = zeros(Ki,My);
AAR = zeros(Ki,M);
Hi_d = zeros(Ki,M);
sin_vertical_AOD = zeros(1,Ki);
cos_horizontal_AOD = zeros(1,Ki);
sin_horizontal_AOD = zeros(1,Ki);
d_antenna = lambda/2;
%lambda = 0.1; % the wavelength of the carrier frequency. 0.1 m

racian_factor = 1.26; % 2dB
sigma = 1;         % the variance of CSCG

for i = 1 :Ki
    sin_vertical_AOD(i) = abs(UAV_location(3)-IU_location(3,i))/dis_IU_AP(i);
    cos_horizontal_AOD(i) = abs(UAV_location(2)-IU_location(2,i))/sqrt( (UAV_location(1)-IU_location(1,i))^2+(UAV_location(2)-IU_location(2,i))^2 );
    sin_horizontal_AOD(i) = abs(UAV_location(1)-IU_location(1,i))/sqrt( (UAV_location(1)-IU_location(1,i))^2+(UAV_location(2)-IU_location(2,i))^2 );
    for j = 1 : Mx
        arg = -2*pi*d_antenna*(j-1)/lambda*sin_vertical_AOD(i)*cos_horizontal_AOD(i);
        AAR_1(i,j) = exp(1i.*arg);
    end
    for k = 1 : My
        arg = -2*pi*d_antenna*(k-1)/lambda*sin_vertical_AOD(i)*sin_horizontal_AOD(i);
        AAR_2(i,k) = exp(1i.*arg);
    end    
    AAR(i,:) = kron(AAR_1(i,:),AAR_2(i,:));

    NLOS_average = 0;
    for rand_num = 1 :200
        NLOS_average = NLOS_average + (randn(1,M)+1i*randn(1,M));
    end
    NLOS_average = NLOS_average/200;

    Hi_d(i,:) = sqrt(racian_factor/(racian_factor+1))*AAR(i,:) + sqrt(1/(racian_factor+1))*(NLOS_average)*sigma/sqrt(2)   ;
end

Hi_d = Hi_d.';

end

