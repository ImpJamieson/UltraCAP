function [phi_l_vec_Next] = SolvePhiL(B_1,phi_l_vec_Current,phi_s_vec_Current,phi_s_vec_Next,lsEnd,llEnd,llVec,PlsPoints,PllPoints)

% SOLVE FOR PHI_L IN THE SOLID PHASE REGION

% Trim vectors to account only for solid region

phi_l_vec_Current_Solid = phi_l_vec_Current(1:PlsPoints+1);

phi_l_vec_Next_Solid = B_1*phi_l_vec_Current_Solid + phi_s_vec_Next - phi_s_vec_Current;

% SOLVE FOR PHI_L IN THE LIQUID PHASE REGION

% Linear potential gradient across seperator

phiL_SolidEnd = phi_l_vec_Next_Solid(end);
phiL_LiquidEnd = 0;

phi_l_vec_Next_Liquid = (phiL_LiquidEnd-phiL_SolidEnd)/(llEnd-lsEnd)*(llVec(PlsPoints+1:PllPoints+1) - lsEnd) + phiL_SolidEnd;
phi_l_vec_Next_Liquid = transpose(phi_l_vec_Next_Liquid(2:end));

% CONCATENATE VECTORS

phi_l_vec_Next = vertcat(phi_l_vec_Next_Solid,phi_l_vec_Next_Liquid);

end