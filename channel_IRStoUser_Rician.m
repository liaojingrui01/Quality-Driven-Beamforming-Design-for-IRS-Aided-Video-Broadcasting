function [Hi_r] = channel_IRStoUser( Nx,Ny,Ki,IRS_location_center,IU_location,lambda  )
%% 生成IRStoUser的Racia信道
%% 需输入横纵反射单元数Nx\Ny,用户数Ki，反射单元中心坐标，用户坐标
dis_IU_IRS = sqrt( sum( abs(IRS_location_center - IU_location ).^2, 1 ) )';

N = Nx*Ny;
AAR_1 = zeros(Ki,Nx);
AAR_2 = zeros(Ki,Ny);
AAR = zeros(Ki,N);
Hi_r = zeros(Ki,N);
sin_vertical_AOD = zeros(1,Ki);
cos_horizontal_AOD = zeros(1,Ki);
sin_horizontal_AOD = zeros(1,Ki);
d_elements = lambda/2;
%lambda = 0.1; % the wavelength of the carrier frequency. 0.1 m

racian_factor = 1.26; % 3dBm
sigma = 1;         % the variance of CSCG

for i = 1 :Ki
    sin_vertical_AOD(i) = abs(IRS_location_center(3)-IU_location(3,i))/dis_IU_IRS(i);
    cos_horizontal_AOD(i) = abs(IRS_location_center(2)-IU_location(2,i))/sqrt( (IRS_location_center(1)-IU_location(1,i))^2+(IRS_location_center(2)-IU_location(2,i))^2 );
    sin_horizontal_AOD(i) = abs(IRS_location_center(1)-IU_location(1,i))/sqrt( (IRS_location_center(1)-IU_location(1,i))^2+(IRS_location_center(2)-IU_location(2,i))^2 );
    for j = 1 : Nx
        arg = -2*pi*d_elements*(j-1)/lambda*sin_vertical_AOD(i)*cos_horizontal_AOD(i);
        AAR_1(i,j) = exp(1i.*arg);
    end
    for k = 1 : Ny
        arg = -2*pi*d_elements*(k-1)/lambda*sin_vertical_AOD(i)*sin_horizontal_AOD(i);
        AAR_2(i,k) = exp(1i.*arg);
    end
    AAR(i,:) = kron(AAR_1(i,:)',AAR_2(i,:)');
    
    NLOS_average = 0;
    for rand_num = 1 :200
        NLOS_average = NLOS_average + (randn(1,N)+1i*randn(1,N));
    end
    NLOS_average = NLOS_average/200;
    
    Hi_r(i,:) = sqrt(racian_factor/(racian_factor+1))*AAR(i,:) + sqrt(1/(racian_factor+1))*(NLOS_average)*sigma/sqrt(2)   ;
end
Hi_r = Hi_r';
end

