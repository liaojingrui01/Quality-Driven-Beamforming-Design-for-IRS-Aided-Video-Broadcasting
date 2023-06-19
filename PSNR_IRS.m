%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mx = 2;
My = 5;
M = Mx*My;           %%  number of transmit antennas at the AP
Nx = 40;
Ny = 10;
N = Nx*Ny;

noise_pow = 0.01; %10^(-10);  %%   noise power: -150 dBm/Hz
power_max = 1.4;
%frequency = 2.4*1e9;      %% carrrier center frequncy 
frequency = 7.5494e+08;
lambda = 3e8/frequency;   %% carrier wavelength = speed of light / frequency
path_loss = (lambda/(4*pi))^2;
penalty = 0.001;
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
channel_max = 100;

snr_IRS_all = zeros(1,channel_max);
Wi_IRS_all = zeros(M,channel_max);
v_all = zeros(N,channel_max);
SNR_list = zeros(50,channel_max);

muu_b1_all = zeros(1,channel_max);
Wi_b1_all = zeros(M,channel_max);
v_b1_all = zeros(N,channel_max);
muu_b2_all = zeros(1,channel_max);
Wi_b2_all = zeros(M,channel_max);
v_b2_all = zeros(N,channel_max);
muu_b3_all = zeros(1,channel_max);
Wi_b3_all = zeros(M,channel_max);
v_b3_all = zeros(N,channel_max);

load('rand_IRS_p14_N400_initial.mat','Wi_randIRS_all','G_all','Hi_r_all','Hi_d_all','v_randIRS_all','power_max','N');

for channal_realization = 1 : channel_max
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%信道信息%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    Hi_d = Hi_d_all(:,:,channal_realization);
    G = G_all(:,:,channal_realization);
    Hi_r = Hi_r_all(:,:,channal_realization);

    Wi_original =  Wi_randIRS_all(:,channal_realization); %ones(M,Ki);
    v_original = v_randIRS_all(:,channal_realization);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%求连续反射角%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    muu_last = 0;
    iteration_num = 0;
    %% 开始SCA
    while(1)
        Hi_original  = (G'*diag(v_original)'*Hi_r + Hi_d)';
        B = zeros(M,Ki);
        A = zeros(Ki,1);
        for i = 1:Ki
            A(i) = Hi_original(i,:)*Wi_original;
            B(:,i) = A(i)*Hi_original(i,:)'+Wi_original;
        end
        
        cvx_begin
        variable Wi_bf(M,1) complex
        variable v_ps(N,1) complex
        variable muu
        expression lower_bound(Ki,1)
        Hi_channel  = (G'*diag(v_ps)'*Hi_r + Hi_d)';
        for i = 1:Ki
            lower_bound(i) = real( B(:,i)'*( A(i)* Hi_channel(i,:)' + Wi_bf ) ) - 0.5*norms( B(:,i))^2 - 0.5*norm( A(i)*Hi_channel(i,:)' - Wi_bf , 2 ) - norm(A(i))^2;
        end
        maximize (  muu + penalty *( 2*real(v_original'*v_ps) - norm(v_original)^2 )  )
        subject  to
        for i = 1:Ki
            lower_bound(i) >= noise_pow*muu;
        end
        for n = 1:N
            norms(v_ps(n),2) <= 1;
        end
        norms(Wi_bf,2) <= sqrt(power_max);
        cvx_end
        cvx_status;
        
        Wi_original =  Wi_bf; %ones(M,Ki);
        v_original = v_ps;
        Hi_original  = (G'*diag(v_original)'*Hi_r + Hi_d)';
        for i = 1:Ki
            snr_item(i) = abs(Hi_original(i,:) *  Wi_bf)^2 / noise_pow;
        end
        min_snr = min(snr_item);
        iteration_num = iteration_num + 1;
        SNR_list(iteration_num,channal_realization) = min_snr;
        
        if abs(muu_last-muu)/abs(muu) <= 10^(-3)
            break;
        end
        muu_last = muu;
    end
    
    snr_IRS_all(channal_realization) = min_snr;
    Wi_all(:,channal_realization) = Wi_original;
    v_all(:,channal_realization) = v_original;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% b = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    b = 1; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end

    Hi_original = (G'*diag(v_discrete)'*Hi_r + Hi_d)';
    [muu_IRSb1, Wi_original_IRSb1] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow);
    
    v_b1_all(:,channal_realization) = v_discrete;
    muu_b1_all(channal_realization) = muu_IRSb1;
    Wi_b1_all(:,channal_realization) = Wi_original_IRSb1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% b = 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    b = 2; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end

    Hi_original = (G'*diag(v_discrete)'*Hi_r + Hi_d)';
    [muu_IRSb2, Wi_original_IRSb2] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow);
    
    v_b2_all(:,channal_realization) = v_discrete;
    muu_b2_all(channal_realization) = muu_IRSb2;
    Wi_b2_all(:,channal_realization) = Wi_original_IRSb2; 
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% b = 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    b = 3; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end

    Hi_original = (G'*diag(v_discrete)'*Hi_r + Hi_d)';
    [muu_IRSb3, Wi_original_IRSb3] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow);
    
    v_b3_all(:,channal_realization) = v_discrete;
    muu_b3_all(channal_realization) = muu_IRSb3;
    Wi_b3_all(:,channal_realization) = Wi_original_IRSb3; 
    
    min_snr
    muu_IRSb3
    muu_IRSb2
    muu_IRSb1
end

average_muu = sum(snr_IRS_all)/channel_max
average_muu_b1 = sum(muu_b1_all)/channel_max
average_muu_b2 = sum(muu_b2_all)/channel_max
average_muu_b3 = sum(muu_b3_all)/channel_max

% alpha = 0;
% c1 = 0.905;
% c2 = 1.34;
% B = 200;
% 
% theta1 = 13870;
% beta1 = 493.2;
% theta2 = 2876;
% beta2 = 23.6;
% 
% Q1=-10*log10( theta1/( c1*B*log2(1+average_muu/c2)-beta1 ) -alpha ) + 20*log10(255)
% Q2=-10*log10( theta2/( c1*B*log2(1+average_muu/c2)-beta2 ) -alpha ) + 20*log10(255)

save('new_IRS_N400_p14','Hi_d_all','Hi_r_all','G_all','N','power_max'...
    ,'snr_IRS_all','Wi_all','v_all','SNR_list','average_muu'...
    ,'v_b1_all','muu_b1_all','Wi_b1_all','average_muu_b1'...
    ,'v_b2_all','muu_b2_all','Wi_b2_all','average_muu_b2'...
    ,'v_b3_all','muu_b3_all','Wi_b3_all','average_muu_b3');


