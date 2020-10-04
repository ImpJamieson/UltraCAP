function [P_2] = P2MatrixFunction(PlsPoints)

%% Initialise blank matrices

P_2 = zeros(PlsPoints,PlsPoints);

%% Populate interior of P_2 matrix

% Inner values

for i = 2:PlsPoints-1
    
    P_2(i,i-1) = 1; % sets below diagonals
    P_2(i,i) = - 2; % sets diagonals
    P_2(i,i+1) = 1; % sets above diagonals
    
end

% Boundary values

        P_2(1,1) = - 2; % sets top corner
        P_2(1,2) = 1; % sets top corner plus one to right
        P_2(end,end) = - 2; % sets bottom corner
        P_2(end,end-1) = 2; % sets bottom corner plus one to left


end