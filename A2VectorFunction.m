function [A_2] = A2VectorFunction(beth,daleth,Lambda,PlsPoints)

%% Initialise blank vector for addition

A_2 = zeros(PlsPoints,1);

%% Enforce boundary conditions through A_2 vector

A_2(1,1) = beth - daleth*(1+Lambda); % sets first value

end