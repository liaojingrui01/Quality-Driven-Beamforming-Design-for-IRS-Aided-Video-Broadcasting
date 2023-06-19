%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mx = 2;
My = 5;
M = Mx*My;           %%  number of transmit antennas at the AP
Nx = 40;
Ny = 10;
N = Nx*Ny;

noise_pow = 0.01; %10^(-10);  %%   noise power: -160 dBm/Hz
power_max = 1.5;
frequency = 2.4*1e9;      %% carrrier center frequncy 2GHz
lambda = 3e8/frequency; %% carrier wavelength = speed of light / frequency
path_loss = (lambda/(4*pi))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成IRS坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRS_location_center = [200 0 25]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成用户坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
IU_location = [    
    168  180  205  210  164
     25   41   59   42   53
      1    1    1    1    1];  
Ki = size(IU_location,2); %%用户坐标数
K = Ki;
AP_location = [0 40 10]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%信道实现%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
    beta_r2(user_index) = 2*10^(-6)*dis_AP_IRS^(-2.2)*dis_IU_IRS(user_index)^(-2.8)./10^(-10);
    beta_d(user_index)  = 10^(-3)*dis_IU_AP(user_index)^(-3.8)/10^(-10);
end
load('no_IRS_p15','Wi_noIRS_all','Hi_d_all','muu_noIRS_all');

muu_randIRS_all = zeros(1,channel_max);
Wi_randIRS_all = zeros(M,channel_max);
v_randIRS_all = zeros(N,channel_max);

for channal_realization = 1 : channel_max

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%信道信息%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    Hi_d = Hi_d_all(:,:,channal_realization);
    
    G = channel_APtoIRS_Rician( Nx,Ny,Mx,My,AP_location,IRS_location_center,lambda  );
    G_all(:,:,channal_realization) = G;

    H_r_rayleigh = channel_IRStoUser_Rician( Nx,Ny,Ki,IRS_location_center,IU_location,lambda );
    Hi_r  = H_r_rayleigh*diag( sqrt(beta_r2) );        %%  IRS-user reflecting link channel   N x K
    Hi_r_all(:,:,channal_realization) = Hi_r;

    Wi_original =  Wi_noIRS_all(:,channal_realization); %ones(M,Ki);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%rand IRS%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    muu_randIRS = 0;
    muu_noIRS = muu_noIRS_all(channal_realization);
    while ( muu_randIRS < (muu_noIRS - 0.1) )
        phz_ini = 0 + (2*pi-0).*rand(N,1);
        phz = phz_ini;
        v_original = exp(1i.*phz);
        Hi_original  = (G'*diag(v_original)'*Hi_r + Hi_d )';
        [muu_randIRS, Wi_original_randIRS] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow);
    end
    muu_randIRS
    
    muu_randIRS_all(channal_realization) = muu_randIRS;
    Wi_randIRS_all(:,channal_realization) = Wi_original_randIRS;
    v_randIRS_all(:,channal_realization) = v_original;
end

average_muu = sum(muu_randIRS_all)/channel_max

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

save('rand_IRS_p15_N400_initial','Wi_randIRS_all','G_all','Hi_r_all','Hi_d_all','muu_randIRS_all','v_randIRS_all','average_muu','Q1','Q2','power_max','N', '-v7.3');

