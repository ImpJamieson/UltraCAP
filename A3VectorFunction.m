function [A_3] = A3VectorFunction(beta,c_0,c_SolidEnd,PlsPoints,PllPoints)

%% Initialise blank vector for addition

A_3 = zeros(PllPoints-PlsPoints-1,1);

%% Enforce boundary conditions through A_5 vector if necessary

        A_3(1,1) = beta*c_SolidEnd; % sets first value
        A_3(end,1) = beta*c_0; % sets final value

end