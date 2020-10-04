function  [rho_vec_Next,rho_net_vec_Next] = SolveRho(rho_vec_Current,R_1,A_1,var_rho,rho_t,connection_type)

switch connection_type
   case 'ConstantCurrent'
        rho_vec_Next = R_1*rho_vec_Current + A_1;
        rho_net_vec_Next = rho_vec_Next - var_rho; % creates net charge density vector
   case 'ConstantVoltage'
        rho_vec_Current = rho_vec_Current(2:end); % trims rho current vector for use
        rho_vec_Next = R_1*rho_vec_Current + A_1;
        rho_vec_Next = vertcat(rho_t,rho_vec_Next); % concatenates vectors 
        rho_net_vec_Next = rho_vec_Next - var_rho; % creates net charge density vector
   otherwise
      disp('ERROR: Specified connection type is unknown.')
end

end