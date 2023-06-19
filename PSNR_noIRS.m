%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mx = 2;
My = 5;
M = Mx*My;           %%  number of transmit antennas at the AP
Nx = 10;
Ny = 10;
N = Nx*Ny;

noise_pow = 0.01; %10^(-10);  %%   noise power: -160 dBm/Hz
power_max = 1.5;
% frequency = 2*1e9;      %% carrrier center frequncy 2GHz
frequency = 7.5494e+08;
lambda = 3e8/frequency; %% carrier wavelength = speed of light / frequency
path_loss = (lambda/(4*pi))^2;
% path_loss = 10^(-3);
penalty = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����IRS����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRS_location_center = [200 0 25]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����û�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
IU_location = [    
    168  180  205  210  164
     25   41   59   42   53
      1    1    1    1    1];  
Ki = size(IU_location,2); %%�û�������
K = Ki;
AP_location = [0 40 10]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ŵ�ʵ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
channel_max = 1000;
Hi_d_all = zeros(M,K,channel_max);
Hi_r_all = zeros(N,K,channel_max);
G_all = zeros(N,M,channel_max);
beta_r2 = zeros(1,K);
beta_d = zeros(1,K);
dis_AP_IRS = sqrt( sum( abs(AP_location - IRS_location_center ).^2, 1 ) )';   % distance between the AP and the IRS, equal for all users
dis_IU_IRS = sqrt( sum( abs( IU_location( : , : , : ) - IRS_location_center ).^2, 1 ) )';  % distance between the IRS and users
dis_IU_AP = sqrt( sum( abs(AP_location - IU_location( : , : , : ) ).^2, 1 ) )'; % distance between the AP and users
for user_index = 1: K
    beta_r2(user_index) = 2*path_loss*path_loss*dis_AP_IRS^(-2.2)*dis_IU_IRS(user_index)^(-2.8)./10^(-10);
    beta_d(user_index)  = path_loss*dis_IU_AP(user_index)^(-3.8)/10^(-10);
end

muu_noIRS_all = zeros(1,channel_max);
Wi_noIRS_all = zeros(M,channel_max);


for channal_realization = 1 : channel_max

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ŵ���Ϣ%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    G = channel_APtoIRS_Rician( Nx,Ny,Mx,My,AP_location,IRS_location_center,lambda  );
    G_all(:,:,channal_realization) = G;

    H_r_rayleigh = channel_IRStoUser_Rician( Nx,Ny,Ki,IRS_location_center,IU_location,lambda );
    Hi_r  = H_r_rayleigh*diag( sqrt(beta_r2) );        %%  IRS-user reflecting link channel   N x K
    Hi_r_all(:,:,channal_realization) = Hi_r;

    H_d_rayleigh =  channel_APtoUser_Rician(Mx,My,Ki,AP_location,IU_location,lambda);
    Hi_d  = H_d_rayleigh*diag( sqrt(beta_d) );         %%  AP-user direct link channel   M x K
    Hi_d_all(:,:,channal_realization) = Hi_d;

%     G = G_all(:,:,channal_realization);
%     Hi_r = Hi_r_all(:,:,channal_realization);
%     Hi_d = Hi_d_all(:,:,channal_realization);

    Wi_original =  ones( M, 1, 1 ); %ones(M,Ki);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%no IRS%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Hi_original  = ( Hi_d )';
    [muu_noIRS, Wi_original_noIRS] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow);
    muu_noIRS
    
    muu_noIRS_all(channal_realization) = muu_noIRS;
    Wi_noIRS_all(:,channal_realization) = Wi_original_noIRS;
end

average_muu = sum(muu_noIRS_all)/channel_max

alpha = 0;
c1 = 0.905;
c2 = 1.34;
B = 200;

theta1 = 13870;
beta1 = 493.2;
Q1=-10*log10( theta1/( c1*B*log2(1+average_muu/c2)-beta1 ) -alpha ) + 20*log10(255)

theta2 = 2876;
beta2 = 23.6;
Q2=-10*log10( theta2/( c1*B*log2(1+average_muu/c2)-beta2 ) -alpha ) + 20*log10(255)

save('no_IRS_p15','Wi_noIRS_all','G_all','Hi_r_all','Hi_d_all','muu_noIRS_all','average_muu','Q1','Q2','power_max','path_loss');

