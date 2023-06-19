function [min_snr, Wi_original] = SCA_beamforming(Hi_original, Wi_original, M, Ki, power_max, noise_pow)

muu_last = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%å¼?§‹SCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)

    A = zeros(Ki,1);
    for i = 1:Ki
        A(i) = norm(Hi_original(i,:)*Wi_original)^2
    end

    cvx_begin
    variable Wi_bf(M,1) complex
    variable muu
    maximize (  muu  )
    subject  to
    for i = 1:Ki
        A(i)*( 2*real( Wi_original'*Wi_bf ) - norm(Wi_original)^2 ) >= noise_pow*muu;
    end
    norms(Wi_bf,2) <= sqrt(power_max);
    cvx_end
    cvx_status;
    Wi_original =  Wi_bf; %ones(M,Ki);

    if abs(muu_last-muu)/abs(muu) <= 10^(-3)
        break;
    end
    muu_last = muu;
end

        for i = 1:Ki
            snr_item(i) = abs(Hi_original(i,:) *  Wi_bf)^2 / noise_pow;
        end
        min_snr = min(snr_item);

end

