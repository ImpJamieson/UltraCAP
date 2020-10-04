clear all
close all
clc

RlPoints = 11;
lStep = 0.5;

F_hat = 1;
sigma_hat = 2;
varepsilon_hat = 4;
sigma = 1;

j_in = 1000; % dummy value
var_rho = -2; % dirichlet condition of rho
bar_rho = 0; % average rho in the domain
e = -1; % dummy value
k_B = 2; %dummy value
T = 1; %dummy value
xi = 1;

Lambda = 0.5;
Upsilon = 1;
Theta = 4;

Pi = (sigma*k_B*T)/(e*var_rho); % value of D at the dirichlet condition

PllPoints = 10;

varkappa_vec = ones(PllPoints,1)*1;

connection_type = 'ConstantCurrent';

[B_1] = B1MatrixFunction(varkappa_vec,PllPoints,connection_type);
