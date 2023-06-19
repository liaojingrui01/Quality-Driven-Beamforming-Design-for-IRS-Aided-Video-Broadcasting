clc; clear all;tic;
%% 2019/08/20
%% compare pdd and SDR, after discussion with Mingmin Zhao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mx = 2;
My = 5;
M = Mx*My;           %%  number of transmit antennas at the AP
Nx = 6;
Ny = 10;
N = Nx*Ny;

frequency = 2.4*1e9;    %% carrrier center frequncy 750MHz
lambda = 3e8/frequency; %% carrier wavelength = speed of light / frequency
path_loss=(lambda/(4*pi));
noise_pow = 10^(-2);    %% back ground noise power  10^(-13)/10^(-10);
% user_SINR_dB  = 20;
% user_SINR     = 10^( user_SINR_dB/10 );
%time_delay = 3.2 * 10^(-7);
time_delay = 2.5 * 10^(-6);
bandwidth = 0.1*10^6;
user_SINR_dB = 2^(1/(time_delay*bandwidth))-1;
user_SINR     = 10^( user_SINR_dB/10 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成IRS坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRS_location_center = [400 700 20]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成用户坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
IU_location = [    
    58  151  145   81  200
   148  137  183  190  200
     0    0    0    0   0]/250*800;  
Ki = size(IU_location,2); %%用户坐标数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成无人机轨迹%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('trajectory2');
T = ceil(sum(delta_r));                  % the number of time slots  
Q1 = zeros(2,T+1);
for t = 1 : T
    for i = 1 : L-1
        if sum(delta_r(1:i)) <= t
            Q1(:,t+1) = Q1(:,t+1) + delta_r(i)*V1_r(:,i);
        elseif sum(delta_r(1:i-1)) <= t && sum(delta_r(1:i)) > t
            %disp([num2str( t-sum(delta_r(1:i-1)) )]);
            Q1(:,t+1) = Q1(:,t+1) + (t-sum(delta_r(1:i-1)))*V1_r(:,i);
        else
            break;
        end
    end
end


H = 80;
T = T + 1;
Q1_r = H*ones(3,T);
Q1_r(1:2,:) = Q1;

T = 1;
channel_realization_max_num = T;
total_pow_IRS_qq = zeros(1,T);
total_pow_randomIRS = zeros(1,T);
total_pow_noIRS = zeros(1,T);
total_pow_discreteIRS_b1 = zeros(1,T);
total_pow_discreteIRS_b2 = zeros(1,T);
total_pow_discreteIRS_b3 = zeros(1,T);
total_pow_discreteIRS_b4 = zeros(1,T);

total_Hi_d = zeros(M,Ki,T);
total_Hi_r = zeros(N,Ki,T);
total_G = zeros(N,M,T);
total_X_ini = zeros(Ki,Ki,T);
total_Wi_bf = zeros(M,Ki,T);
total_v_original = zeros(N,T);

outer_layer_obj = zeros(T,500);
all_inner_obj = zeros(T,1000);
for channel_num = 1: channel_realization_max_num  %% average the fading channel
    
    AP_location = Q1_r(:,channel_num);
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['the channel realization number : ', num2str( channel_num )]);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
   %% 求AP-user信道
    Hi_d_AAR =  channel_UAVtoUser_Rician( Mx,My,Ki,AP_location,IU_location,lambda );
    Hi_d_AAR =  channel_cof2(M,Ki,1);
    dis_IU_AP = sqrt( sum( abs(AP_location - IU_location( : , : ) ).^2, 1 ) )';
    beta_i_d = path_loss/10^(-10)*dis_IU_AP.^(-3.5);      %%  the distance-dependent path loss UAV-user
    Hi_d  = Hi_d_AAR*diag( sqrt(beta_i_d) )*1;            %%  AP-user direct link channel   M x K
    
    %% 求IRS-user信道
    Hi_r_AAR =  channel_IRStoUser_Rician( Nx,Ny,Ki,IRS_location_center,IU_location,lambda  );
    dis_AP_IRS = sqrt( sum( abs(AP_location - IRS_location_center ).^2, 1 ) )';
    dis_IU_IRS  = zeros(1,Ki);
    beta_i_r2 = zeros(1,Ki);                                                                                 %% the distance-dependent path loss UAV-AP-user
    for user_index = 1: Ki   
        dis_IU_IRS( user_index ) = sqrt( sum( abs( IU_location( user_index ) - IRS_location_center ).^2, 1 ) )';
        beta_i_r2( user_index )  = 2*path_loss*path_loss/10^(-10)*dis_AP_IRS^(-3)*dis_IU_IRS( user_index ).^(-3);    %% 10^(-6)/10^(-4);
        Hi_r(:,user_index )  = diag( sqrt(beta_i_r2( user_index )) )*Hi_r_AAR(:,user_index );               %% IRS-user reflecting link channel  N x K
    end    
    
    %% 求AP-IRS信道
    G = channel_UAVtoIRS_LoS( Nx,Ny,Mx,My,AP_location,IRS_location_center,lambda  );
    
    Hi_d;
    Hi_r;
    G;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%              Initialize the double loop algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho =  0.01;    %% penalty coefficient
    rho_ini = rho;
    c = 0.95;       %% scaling coefficient of rho
    obj = 0;
    outer_layer_iteration_num = 0;
    all_inner_iteration_num = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% initialize BCD variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_ini = channel_cof2( Ki, Ki, 1 ); %ones( Ki, Ki, 1 ); 
    X = X_ini;
    phz_ini = 0 + (2*pi-0).*rand(N,1);
    phz = phz_ini;
    v_original = exp(1i.*phz);
    Wi_bf =  channel_cof2( M, Ki, 1 ); %ones(M,Ki); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%              Start the double loop algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(1) %% double loop iterative algorithm
        inner_iteration_num = 0;
        inner_layer_obj_decrease = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(1)  %% inner loop: using AO/BCD to alternatly optimize the W, v,and X.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp('%%%%%%%%%%%%%%%%%  inner layer: setp 1  %%%%%%%%%%%%%%%%%%%%%%%');
            inner_iteration_num = inner_iteration_num + 1;
            
            %%% Step 1: optimize the transmit precoding vectors, with given { X, v }
            %%% optimize { Wi_bf }, by using closed-form expression
            Hi_eff_chl  = G'*diag(v_original)'*Hi_r + Hi_d;
            A_bf_quad_cof  = eye(M);
            
            for i = 1:Ki     %%  caculate the quadratic term coefficient for both IDRs and EHRs
                A_bf_quad_cof  = A_bf_quad_cof + Hi_eff_chl(:,i)*Hi_eff_chl(:,i)'/(2*rho);
            end
            
            bi_bf_linr_cof = zeros(M,Ki);
            for i = 1:Ki     %% caculate the linear term coefficient for IDRs
                for k = 1:Ki
                    bi_bf_linr_cof( :, i ) = bi_bf_linr_cof( :, i ) +  Hi_eff_chl(:,k)*(  X(k,i)   );
                end
            end
            
            for i = 1:Ki     %% caculate the information precoder
                Wi_bf(:,i) = ( A_bf_quad_cof\bi_bf_linr_cof(:,i) )/(2*rho);
            end

            %%% caculate the AL objective in this step1
            penalty_term = 0;
            %load bf_data
            for i = 1:Ki   % caculate the penalty term in the augmented lagrangian function for only IDRs
                for k = 1:Ki
                    penalty_term = penalty_term + abs( Hi_eff_chl(:,i)'*Wi_bf(:,k) - X(i,k) )^2;
                end
            end
            transmit_pow_step1  =  real(  trace( Wi_bf*Wi_bf' ) );
            obj_AL_step1 = transmit_pow_step1 + 1/(2*rho)*(  penalty_term );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp('%%%%%%%%%%%%%%%%%  inner layer: setp 2  %%%%%%%%%%%%%%%%%%%%%%%');
            inner_iteration_num = inner_iteration_num + 1;
            
            %%% Step 2: optimize IRS phase shifts, given { Wi_bf, X }
            %%% optimize phz_shift, v, by using one-iteration BCD
            X_b = zeros(Ki,Ki);
            X_a = zeros(N,Ki,Ki);
            A_phz_quad_cof = zeros(N,N);
            b_phz_linr_cof = zeros(N,1);
            v = conj( v_original );
            
            for i = 1:Ki   % construct one-iteration BCD coefficient as in the TWC paper
                for k = 1:Ki
                    X_b(i,k)   =  Hi_d(:,i)'*Wi_bf(:,k) - X(i,k);
                    X_a(:,i,k) =  diag( Hi_r(:,i)' )*G*Wi_bf(:,k);
                end
            end
            
            %%% caculate the QCQP coefficient for the final expression of phz shift optimization
            for i = 1:Ki   %
                for k = 1:Ki
                    A_phz_quad_cof = A_phz_quad_cof + X_a(:,i,k)*X_a(:,i,k)';
                    b_phz_linr_cof = b_phz_linr_cof + X_a(:,i,k)*X_b(i,k)';
                end
            end
            
            %%% caculate phase shifts one by one
            obj_phz_last = 0;
            phz_iteration_num = 0;
            phz_obj_converge = [];
            while (1)
                for n = 1:N
                    n;
                    phz_iteration_num = phz_iteration_num + 1;
                    AX = A_phz_quad_cof*v;
                    coefficient_temp = AX(n) - A_phz_quad_cof(n,n)*v(n) + b_phz_linr_cof(n);
                    if coefficient_temp == 0
                        v_n = 1;
                    else
                        v_n =  -coefficient_temp/abs( coefficient_temp );    %% minimize AL function, so take additional pi phase shift
                        v(n) = v_n;
                    end
                    obj_phz = real( v'*A_phz_quad_cof*v + 2*real( v'*b_phz_linr_cof ) );
                    phz_obj_converge(phz_iteration_num ) = obj_phz;
                end
                obj_phz = real( v'*A_phz_quad_cof*v + 2*real( v'*b_phz_linr_cof ) );
                if real( obj_phz - obj_phz_last ) > -10^(-4)  && real( obj_phz - obj_phz_last ) <= 10^(-3)
                    %disp(['%%%%%%  One-iteration BCD for phase shift converges!  %%%%%%%']);
                    v;
                    break
                else
                    obj_phz_last = obj_phz;
                end
                
            end
            v_original = conj( v );    %%  column original theta vector
            
            %%% caculate the AL objective in this step2
            Hi_eff_chl  = G'*diag(v_original)'*Hi_r + Hi_d;
            penalty_term = 0;
            %load bf_data
            for i = 1:Ki   % caculate the penalty term in the augmented lagrangian function for only IDRs
                for k = 1:Ki
                    penalty_term = penalty_term + abs( Hi_eff_chl(:,i)'*Wi_bf(:,k) - X(i,k) )^2;
                end
            end
            
            transmit_pow_step2  =  real(  trace( Wi_bf*Wi_bf' ) );
            obj_AL_step2 = transmit_pow_step2 + 1/(2*rho)*(  penalty_term  );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp('%%%%%%%%%%%%%%%%%  inner layer: setp 3  %%%%%%%%%%%%%%%%%%%%%%%');
            inner_iteration_num = inner_iteration_num + 1;
            
            %%% Step 2: optimize the newly introduced auxiliary variables, given { Wi_bf, v }
            %%% optimize { X } by using QCQP with one single constraint
            X_backup = X;
            X = zeros(Ki,Ki);
            C_x = zeros(Ki,Ki);
            
            for i = 1:Ki   % construct constants in QCQP with one single constraint
                for k = 1:Ki
                    C_x(i,k) = (  Hi_r(:,i)'*diag(v_original)*G + Hi_d(:,i)' )*Wi_bf(:,k) ;
                end
            end
            
            for i = 1:Ki   % caculate newly introduced variables { X }
                [ X(i,:) ]  =  QCQP_with_SingleConstra_X( i, C_x, user_SINR, noise_pow );
            end
            
            %%% caculate the AL objective in this step3
            Hi_eff_chl  = G'*diag(v_original)'*Hi_r + Hi_d;
            penalty_term = 0;
            %load bf_data
            for i = 1:Ki   % caculate the penalty term in the augmented lagrangian function for only IDRs
                for k = 1:Ki
                    penalty_term = penalty_term + abs( Hi_eff_chl(:,i)'*Wi_bf(:,k) - X(i,k) )^2;
                end
            end
            
            transmit_pow_step3   =  real(  trace( Wi_bf*Wi_bf' ) );
            obj_AL_step3 = transmit_pow_step3 + 1/(2*rho)*(  penalty_term  );
            
            all_inner_iteration_num = all_inner_iteration_num + 1;
            all_inner_obj(channel_num,all_inner_iteration_num) = transmit_pow_step3;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp('%%%%%%%%%%%%%%%%%  inner layer: overview AL obj  %%%%%%%%%%%%%%%%%%%%%%%');
            %%% to check if the inner loop BCD converges to the pre-defined accuracy
            if inner_iteration_num>=2 & abs(obj_AL_step3 - obj_AL_step1)<=1e-3*obj_AL_step1
                break;
            end
            if inner_iteration_num >= 20
                break;
            end
        end  %%% end of inner loop: using AO/BCD to alternatly optimize the W,theta, and X.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %outer loop: decrease the penalty term rho in PDD so as to force the equality constraint
        %disp('%%%%%%%%%%%%%  decrease rho to increase the penalty term value  %%%%%%%%%%%%%% ');
        %constraint violation for newly introduced equalities
        X_violation = max( max(  abs( X - C_x  )  ) );
        max_violation = max( [ X_violation ] );
        rho = c*rho;
        outer_layer_iteration_num  =  outer_layer_iteration_num + 1;
        outer_layer_obj(channel_num,outer_layer_iteration_num) = transmit_pow_step3;
        if max_violation < 1e-6
            %disp('%%%%  the outer layer optimization converges  %%%%');
            rho;
            rho_ini;
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%%%      caculate the random IRS upper bound by using SDR       %%%%%%%%%%%%%%%%%%%%'])
    phz_sum = zeros(N,1);
    for i = 1:200
        phz_rand = 0 + (2*pi-0).*rand(N,1);
        phz_sum = phz_sum + phz_rand;
    end
    phz_ave = phz_sum/200;
    v1 = exp(1i.*phz_ave);
    Hi_eff_chl_randomIRS  =  G'*diag(v1)'*Hi_r + Hi_d;
    pow_randomIRS = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_randomIRS, M, Ki);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%      caculate the no IRS upper bound by using SDR      %%%%%%%%%%%%%%%%%%%']) 
    Hi_eff_chl_noIRS  =  Hi_d;
    pow_noIRS = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_noIRS, M, Ki);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%      caculate the discrete IRS b=1 upper bound by using SDR      %%%%%%%%%%%%%%%%%%%']) 
    b = 1; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end
    Hi_eff_chl_discreteIRS_b1  =  G'*diag(v_discrete)'*Hi_r + Hi_d;
    pow_discreteIRS_b1 = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_discreteIRS_b1, M, Ki);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%      caculate the discrete IRS b=2 upper bound by using SDR      %%%%%%%%%%%%%%%%%%%']) 
    b = 2; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end
    Hi_eff_chl_discreteIRS_b2  =  G'*diag(v_discrete)'*Hi_r + Hi_d;
    pow_discreteIRS_b2 = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_discreteIRS_b2, M, Ki);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%      caculate the discrete IRS b=3 upper bound by using SDR      %%%%%%%%%%%%%%%%%%%']) 
    b = 3; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end
    Hi_eff_chl_discreteIRS_b3  =  G'*diag(v_discrete)'*Hi_r + Hi_d;
    pow_discreteIRS_b3 = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_discreteIRS_b3, M, Ki);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disp(['%%%%%%%%%%%%%%%%%%      caculate the discrete IRS b=4 upper bound by using SDR      %%%%%%%%%%%%%%%%%%%']) 
    b = 4; %%离散个数2^b
    phase_shift = linspace(0,2*pi,2^b);
    v_discrete = ones(N,1);
    for i = 1:N
        [m,index]=min( abs(exp(1i.*phase_shift)-v_original(i)) );
        v_discrete(i)=exp(1i.*phase_shift(index));
    end
    Hi_eff_chl_discreteIRS_b4  =  G'*diag(v_discrete)'*Hi_r + Hi_d;
    pow_discreteIRS_b4 = SDR_benchmark( user_SINR, noise_pow, Hi_eff_chl_discreteIRS_b4, M, Ki);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    transmit_pow_step3
    pow_randomIRS
    pow_noIRS
    pow_discreteIRS_b1
    pow_discreteIRS_b2
    pow_discreteIRS_b3
    pow_discreteIRS_b4
    
    total_pow_randomIRS(channel_num) = pow_randomIRS;
    total_pow_noIRS(channel_num) = pow_noIRS;
    total_pow_discreteIRS_b1(channel_num) = pow_discreteIRS_b1;
    total_pow_discreteIRS_b2(channel_num) = pow_discreteIRS_b2;
    total_pow_discreteIRS_b3(channel_num) = pow_discreteIRS_b3;
    total_pow_discreteIRS_b4(channel_num) = pow_discreteIRS_b4;
    total_pow_IRS_qq(channel_num) = transmit_pow_step3;

    total_Hi_d(:,:,channel_num) = Hi_d;
    total_Hi_r(:,:,channel_num) = Hi_r;
    total_G(:,:,channel_num) = G;
    total_X_ini(:,:,channel_num) = X;
    total_Wi_bf(:,:,channel_num) = Wi_bf;
    total_v_original(:,channel_num) = v_original;

end  %%  for channel_num = 1: channel_realization_max_num  %% average the fading channel

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
sum(total_pow_IRS_qq)
sum(total_pow_randomIRS)
sum(total_pow_noIRS)
sum(total_pow_discreteIRS_b1)
sum(total_pow_discreteIRS_b2)
sum(total_pow_discreteIRS_b3)
sum(total_pow_discreteIRS_b4)

save('initial_M10_delay25_N60.mat','all_inner_obj','outer_layer_obj','total_pow_IRS_qq','total_pow_randomIRS','total_pow_noIRS',...
    'total_pow_discreteIRS_b1','total_pow_discreteIRS_b2','total_pow_discreteIRS_b3','total_pow_discreteIRS_b4',...
    'total_X_ini','total_v_original','total_Wi_bf','total_Hi_d','total_Hi_r','total_G','N','T','M','time_delay','bandwidth');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


a