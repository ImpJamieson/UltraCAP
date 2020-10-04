function [c_vec_Next] = SolveC(C_1,C_1_var,alpha,beta,c_vec_Current,c_0,phi_s_vec_Current,phi_s_vec_Next,phi_l_vec_Current,phi_l_vec_Next,PlsPoints,PllPoints)

% SOLVE FOR C IN THE SOLID PHASE REGION

% Trim vectors to account only for solid region

c_vec_Current_Solid = c_vec_Current(1:PlsPoints+1);
phi_s_vec_Current_Solid = phi_s_vec_Current(1:PlsPoints+1);
phi_s_vec_Next_Solid = phi_s_vec_Next(1:PlsPoints+1);
phi_l_vec_Current_Solid = phi_l_vec_Current(1:PlsPoints+1);
phi_l_vec_Next_Solid = phi_l_vec_Next(1:PlsPoints+1);

c_vec_Next_Solid = C_1*c_vec_Current_Solid + alpha*phi_s_vec_Next_Solid - alpha*phi_s_vec_Current_Solid - alpha*phi_l_vec_Next_Solid + alpha*phi_l_vec_Current_Solid;

% SOLVE FOR C IN THE LIQUID PHASE REGION

% Trim vector to account only for liquid region

c_vec_Current_Liquid = c_vec_Current(PlsPoints+2:end-1);

% Generate Dirichlet condition and add to A_3 Vector

c_SolidEnd = c_vec_Next_Solid(end);

[A_3] = A3VectorFunction(beta,c_0,c_SolidEnd,PlsPoints,PllPoints);

c_vec_Next_Liquid = C_1_var*c_vec_Current_Liquid + A_3;

c_vec_Next_Liquid = vertcat(c_vec_Next_Liquid,c_0);

% CONCATENATE VECTORS

c_vec_Next = vertcat(c_vec_Next_Solid,c_vec_Next_Liquid);

end