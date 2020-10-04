function [phi_s_vec_Next] = SolvePhiS(rho_net_vec_Next,K)

% Use specific capacitance equation

phi_s_vec_Next = rho_net_vec_Next/K;

end