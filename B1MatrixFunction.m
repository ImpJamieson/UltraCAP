function [B_1] = B1MatrixFunction(varkappa_vec,PlsPoints)

%% Initialise blank matrix

B_1 = zeros(PlsPoints+1,PlsPoints+1);

%% Populate interior of R_5 matrix

% Inner values

for i = 2:PlsPoints
    
    B_1(i,i-1) = varkappa_vec(i); % sets above diagonals
    B_1(i,i) = 1 - 2*varkappa_vec(i); % sets diagonals
    B_1(i,i+1) = varkappa_vec(i); % sets below diagonals
    
end

% Boundary values

        B_1(1,1) = 1 - 2*varkappa_vec(1); % sets top corner
        B_1(1,2) = 2*varkappa_vec(1); % sets top corner plus one to right
        B_1(end,end) = 1 + (35/12*varkappa_vec(PlsPoints+1)); % sets below top corner
        B_1(end,end-1) = -104/12*varkappa_vec(PlsPoints+1); % sets bottom corner plus one to left
        B_1(end,end-2) = 114/12*varkappa_vec(PlsPoints+1); % sets bottom corner plus two to left
        B_1(end,end-3) = -56/12*varkappa_vec(PlsPoints+1); % sets bottom corner plus three to left
        B_1(end,end-4) = 11/12*varkappa_vec(PlsPoints+1); % sets bottom corner plus four to left

end