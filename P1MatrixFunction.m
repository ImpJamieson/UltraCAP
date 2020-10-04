function [P_1_inv] = P1MatrixFunction(Lambda,PlsPoints)

%% Initialise blank matrices

P_1 = zeros(PlsPoints,PlsPoints);

%% Populate interior of P_1 matrix

% Inner values

for i = 2:PlsPoints-1
    
    P_1(i,i-1) = 1+Lambda; % sets below diagonals
    P_1(i,i) = - 2*(1+Lambda); % sets diagonals
    P_1(i,i+1) = 1+Lambda; % sets above diagonals
    
end

% Boundary values

        P_1(1,1) = - 2*(1+Lambda); % sets top corner
        P_1(1,2) = 1+Lambda; % sets top corner plus one to right
        P_1(end,end) = - 2*(1+Lambda); % sets bottom corner
        P_1(end,end-1) = 2*(1+Lambda); % sets bottom corner plus one to left

%% Compute inverse of P_1 matrix

P_1_inv = inv(P_1);

end