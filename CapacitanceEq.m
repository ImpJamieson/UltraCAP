function [beth,daleth] = CapacitanceEq(phi_s_vec_Current,I_terminal,C_total,tStep)

beth = phi_s_vec_Current(1);

daleth = I_terminal*tStep/C_total + beth;

end