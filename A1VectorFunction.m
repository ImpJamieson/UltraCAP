function [A_1] = A1VectorFunction(Sigma_vec,var_rho,Upsilon,Theta,RlsPoints,connection_type)

%% Initialise blank vector for addition

A_1 = ones(RlsPoints,1)*(-2*Upsilon);

%% Enforce boundary conditions through A_1 vector if necessary

switch connection_type
   case 'ConstantCurrent'
        A_1(end,1) = - 2*Upsilon; % sets final value
   case 'ConstantVoltage'
        A_1(1,1) = - 2*Upsilon + var_rho*(Theta/4 + Sigma_vec(2) - Sigma_vec(3)/4); % sets first value
   otherwise
      disp('ERROR: Specified connection type is unknown.')
      disp('Cannot correctly populate A1 Vector.')
      error('An error occured. Simulation halted. See the above messages for information.')
end