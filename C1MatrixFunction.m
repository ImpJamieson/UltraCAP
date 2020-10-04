function [C_1,C_1_var] = C1MatrixFunction(beta,PlsPoints,PllPoints)

%% Initialise blank matrix

C_1 = zeros(PlsPoints+1,PlsPoints+1);

%% Populate interior of C_1 matrix

% Inner values

for i = 2:PlsPoints
    
    C_1(i,i-1) = beta; % sets above diagonals
    C_1(i,i) = 1 - 2*beta; % sets diagonals
    C_1(i,i+1) = beta; % sets below diagonals
    
end

% Boundary values

        C_1(1,1) = 1 - 2*beta; % sets top corner
        C_1(1,2) = 2*beta; % sets top corner plus one to right
        C_1(end,end) = 1 + (35/12*beta); % sets below top corner
        C_1(end,end-1) = -104/12*beta; % sets bottom corner plus one to left
        C_1(end,end-2) = 114/12*beta; % sets bottom corner plus two to left
        C_1(end,end-3) = -56/12*beta; % sets bottom corner plus three to left
        C_1(end,end-4) = 11/12*beta; % sets bottom corner plus four to left

%% Initialise blank matrix

C_1_var = zeros(PllPoints-PlsPoints-1,PllPoints-PlsPoints-1);

%% Populate interior of C_1 Var matrix

% Inner values

for i = 2:(PllPoints-PlsPoints-2)
    
    C_1_var(i,i-1) = beta; % sets above diagonals
    C_1_var(i,i) = 1 - 2*beta; % sets diagonals
    C_1_var(i,i+1) = beta; % sets below diagonals
    
end

% Boundary values

        C_1_var(1,1) = 1 - 2*beta; % sets top corner
        C_1_var(1,2) = beta; % sets top corner plus one to right
        C_1_var(end,end) = 1 - 2*beta; % sets below top corner
        C_1_var(end,end-1) = beta; % sets below top corner plus one to right

end