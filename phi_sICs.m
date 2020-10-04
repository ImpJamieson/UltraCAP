function [phi_s_vec_Current] = phi_sICs(rho_net_vec_Current,K)

% Use specific capacitance equation

phi_s_vec_Current = rho_net_vec_Current/K;

end