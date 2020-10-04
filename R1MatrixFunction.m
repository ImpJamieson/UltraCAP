function [R_1] = R1MatrixFunction(Sigma_vec,Theta,e,sigma,k_B,T,j_in,lsStep,RlsPoints,connection_type)

%% Initialise blank matrix

R_1 = zeros(RlsPoints,RlsPoints);

%% Create Xi value

Xi = (sigma*k_B*T)/((sigma*k_B*T/Sigma_vec(2))-(2*lsStep*j_in*e/Sigma_vec(1)));

%% Populate interior of R_5 matrix

switch connection_type
    
case 'ConstantCurrent'
    
   % Inner values

    for i = 2:(RlsPoints-1)
    
    R_1(i,i-1) = Sigma_vec(i-1)/4 + Sigma_vec(i) - Sigma_vec(i+1)/4; % sets above diagonals
    R_1(i,i) = 1; % sets diagonals
    R_1(i,i+1) = - Sigma_vec(i-1)/4 + Sigma_vec(i) + Sigma_vec(i+1)/4; % sets below diagonals
    
    end
    
    % Boundary values
   
    R_1(1,1) = 1 - (2*lsStep*j_in*e/sigma/k_B/T)*(Xi/4 + Sigma_vec(1) - Sigma_vec(2)/4); % sets top corner
    R_1(1,2) = 2*Sigma_vec(1); % sets top corner plus one to right
    R_1(end,end) = 1; % sets below top corner
    R_1(end,end-1) = 2*Sigma_vec(end); % sets below top corner plus one to left
        
   case 'ConstantVoltage'
       
    % Inner values

    for i = 2:(RlsPoints-1)
    
    R_1(i,i-1) = Sigma_vec(i)/4 + Sigma_vec(i+1) - Sigma_vec(i+2)/4; % sets above diagonals
    R_1(i,i) = 1; % sets diagonals
    R_1(i,i+1) = - Sigma_vec(i)/4 + Sigma_vec(i+1) + Sigma_vec(i+2)/4; % sets below diagonals
    
    end
    
    % Boundary values
    
    R_1(1,1) = 1; % sets top corner
    R_1(1,2) = -Theta/4 + Sigma_vec(2) + Sigma_vec(3)/4; % sets top corner plus one to right
    R_1(end,end) = 1; % sets below top corner
    R_1(end,end-1) = 2*Sigma_vec(end); % sets below top corner plus one to left
    
   otherwise
      disp('ERROR: Specified connection type is unknown.')
      disp('Cannot correctly populate R1 Matrix.')
end

end