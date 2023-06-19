function [ penalty_term_no_rho ] = Caculate_Augmented_Lagran_dual( Ki, Ke, Hi_eff_chl, He_eff_chl, Wi_bf, We_bf,X, T,S, rho,...
                                       lambda_X, lambda_S, lambda_T  )
%CACULATE_AUGMLAGRAN 
% the function caculate the augmented lagrangian function




penalty_term = 0;
for i = 1:Ki   % caculate the penalty term in the augmented lagrangian function for only IDRs
    for k = 1:Ki
        penalty_term = penalty_term + abs( Hi_eff_chl(:,i)'*Wi_bf(:,k) - X(i,k) + rho*lambda_X(i,k)  )^2;
    end
end
for j = 1:Ke   % caculate the penalty term in the augmented lagrangian function for only EHRs
    for m = 1:Ke
        penalty_term = penalty_term + abs( He_eff_chl(:,j)'*We_bf(:,m) - T(j,m) + rho*lambda_T(j,m)  )^2;
    end
end
% for i = 1:Ki   % caculate the penalty term in the augmented lagrangian function for corss terms
%     for j = 1:Ke
%         penalty_term = penalty_term + abs( Hi_eff_chl(:,i)'*We_bf(:,j) - Y(i,j) + rho*lambda_Y(i,j)  )^2;
%     end
% end

for j = 1:Ke   % caculate the penalty term in the augmented lagrangian function for corss terms
    for i = 1:Ki
        penalty_term = penalty_term + abs( He_eff_chl(:,j)'*Wi_bf(:,i) - S(j,i) + rho*lambda_S(j,i)  )^2;
    end
end

penalty_term_no_rho = penalty_term;





end

