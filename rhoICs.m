function [rho_vec_Current,rho_net_vec_Current] = rhoICs(var_rho,K,V_start,lsPoints)

rho_add = -V_start*K; % why the division by 2?
rho_vec_Current = ones(lsPoints,1)*var_rho + rho_add;
rho_net_vec_Current = rho_vec_Current - var_rho; % creates net charge density vector

end