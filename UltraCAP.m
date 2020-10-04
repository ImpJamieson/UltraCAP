%% Clear data for initialisation

close all
clear variables
clc

%% Input parameters

% TIME PARAMETERS

tBeg = 0; % time always starts at zero [s]
tEnd = 60000; % simulation end time, equal to maximum permissible charge duration [s]
tSteps = 6000000; % number of intervals for time to be computed at
tStep = (tEnd-tBeg)/tSteps; % time step between each calculation [s]
tVec = tBeg:tStep:tEnd; % temporal vector
tPoints = length(tVec); % number of time points to be computed

% GEOMETRICAL PARAMETERS

% Solid phase domain

lsBeg = 0; % solid phase length starts at zero [m]
lsEnd = 50e-6; % solid phase end length [m]
lsSteps = 1000; % number of intervals for solid phase points to be computed at
lsStep = (lsEnd - lsBeg)/lsSteps; % spatial step between points [m]
lsVec = lsBeg:lsStep:lsEnd; % solid phase spatial vector
lsPoints = length(lsVec); % number of solid phase spatial points to be computed

% Liquid phase domain

llBeg = 0; % liquid phase length starts at zero [m]
llEnd = 62.5e-6; % liquid phase end length [m]
llStep = lsStep; % spatial step must be the same in both domains
llVec = llBeg:llStep:llEnd; % liquid phase spatial vector
llPoints = length(llVec); % number of liquid phase spatial points to be computed

% VARYING PHYSICAL PARAMETERS

sigma = 0.0521; % electrode conductivity [S m^{-1}]
Free_E = 1.0479e29; % free electron density of conductor [m^-3]
D_0 = 3.5e-11; % Diffusion coefficient (free solution) in liquid phase [m^{2}s^{?1}]
c_0 = 930; % base ionic concentration in liquid phase [mol m^{?3}]
T = 298.15; % ambient absolute temperature [K]
I_in = 1; % current supply in [A, C s^-1]
V_t = -2.0; % target voltage across terminals [V]
A_xc = 1.5; % electrode cross sectional area [m^{2}] chosen to match Verbrugge
lambda_s = 0.3e-09; % Stern layer thickness [m]
varepsilon_l = 36.6; % relative permittivity of solvent in liquid phase
a_d = 3.89e07; % specific surface area [m^{2} m^{-3}]
epsilon = 0.67; % void fraction
K = 2.1e+09; % effective specific capacitance 2.3333e+09
V_start = 1.5; % capacitor starting voltage (V)

% CONSTANT PHYSICAL PARAMETERS

k_B = 1.380649e-23; % Boltzmann's constant [J K^{-1}]
N_A = 6.02214076e23; % Avogadro constant [mol^{-1}]
e = -1.6e-19; % electron elementary charge [C]
varepsilon_0 = 8.85418782e-12; % permittivity of free space [F m^{-1}]
varepsilon = 0.1750; % permittivity of electrode [F m^{-1}]

% CAPACITANCE MODEL CALCULATION

C_D = varepsilon_l*varepsilon_0/lambda_s; % Helmholtz prediction of capacitance

% BOUNDARY CONDITION SPECIFICATION

% Choose either 'ConstantCurrent' or 'ConstantVoltage'

connection_type = 'ConstantCurrent';

%% Precision initialisation and conversion

% CHOOSE FLOATING POINT NUMBER TYPE

% Select either 'DoublePrecision' or 'VariablePrecision'

number_type = 'DoublePrecision';

switch number_type
    case 'VariablePrecision'
        disp('Variable precision floating point arithmetic chosen.')
        
% SET PRECISION LEVEL

digits(32) % number of bits in floating point calculations

        disp('Precision is:')
        disp(digits)

% CONVERT RELEVANT INPUT PARAMETERS TO SYMBOLIC WITH VARIABLE PRECISION

tStep = vpa(tStep);
tVec = vpa(tVec);
lsStep = vpa(lsStep);
lsVec = vpa(lsVec);
sigma = vpa(sigma);
Free_E = vpa(Free_E);
T = vpa(T);
k_B = vpa(k_B);
e = vpa(e);
varepsilon_0 = vpa(varepsilon_0);
I_in = vpa(I_in);
V_sup = vpa(V_sup);
d_xc = vpa(d_xc);

    case 'DoublePrecision'
        disp('Double precision floating point arithmetic chosen.')
end

%% Derived parameters

% DERIVED SIMULATION PARAMETERS

j_in = I_in/A_xc; % current supply in [C s^-1 m^-2]
rho_0 = Free_E*e; % base level of charge concentration [C m^-3]
var_rho = rho_0; % value of rho at i=N [C m^-3]
Gamma = tStep/(lsStep)^2; % Fourier number
D_prime = sigma*k_B*T/e; % diffusivity paramter collection of coefficients
Lambda = sigma*tStep/varepsilon; % dimensionless number
Upsilon = sigma*k_B*T/e*Gamma; % collection of values
Theta = sigma*k_B*T/e/var_rho*Gamma; % collection of values
xi = varepsilon_0/(lsStep)^2; % collection of values
alpha = a_d*C_D/2/abs(e)/N_A/epsilon;
beta = D_0*Gamma/epsilon;
tStep_Critical = varepsilon_0/sigma; % maximum time step for numerical stability
C_total = a_d*C_D*A_xc*lsEnd; % total effective capacitance

rho_t = 0; % TEST VALUE FOR NOW, NEEDS TO BE CHANGED FOR CONSTANT VOLTAGE

% OTHER DERVIED PARAMETERS

R = 1/sigma*lsEnd; % effective resistance (Ohm)
C = varepsilon/lsEnd; % effective capacitance (Farad)
TheoreticalTimeConstant = R*C; % effective time constant (s)
ComputationalTimeConstant = varepsilon/sigma; % computational time constant (s)

%% Preparations to save output data

tLinesNumber = 10; % number of time evolution lines to plot
tInterval = (tPoints-1)/tLinesNumber; % time intervals to plot at

for i = 1:tLinesNumber+1
    
    if i ==1
        
        tEvolutionPoints = [];
        
    end

tSelect = round(i*tInterval - (tInterval-1));
tEvolutionPoints = vertcat(tEvolutionPoints,tSelect);

end

phi_s_mat = zeros(length(tEvolutionPoints),length(lsVec));
rho_mat = zeros(length(tEvolutionPoints),length(lsVec));
ls_mat = zeros(length(tEvolutionPoints),length(lsVec));
phi_l_mat = zeros(length(tEvolutionPoints),length(llVec));
c_mat = zeros(length(tEvolutionPoints),length(llVec));
ll_mat = zeros(length(tEvolutionPoints),length(llVec));

%% Initial conditions

% CALL INITIAL CONDITION FUNCTIONS

[rho_vec_Current,rho_net_vec_Current] = rhoICs(var_rho,K,V_start,lsPoints); % creates rho initial conditions
[phi_s_vec_Current] = phi_sICs(rho_net_vec_Current,K); % creates phi_s initial conditions
[phi_l_vec_Current] = phi_lICs(llPoints); % creates phi_l initial conditions
[c_vec_Current,c_net_vec_Current] = cICs(c_0,llPoints); % creates c initial conditions
D_vec_Current = D_prime./rho_vec_Current; % sets initial condition for the diffusivity parameter at all x
Sigma_vec = D_vec_Current*Gamma;
kappa_vec = 2*e^2*N_A/k_B/T*D_0*c_vec_Current; % vector of intitial values of kappa collection of values
varkappa_vec = kappa_vec*Gamma/a_d/C_D; % initial values of varkappa collection of values
v_terminal = phi_s_vec_Current(1); % initial value of v_terminal

switch connection_type
   case 'ConstantCurrent'
       RlsPoints = lsPoints; % trims number of length points for solid phase rho matrices
       PlsPoints = lsPoints-1; % trims number of length points for solid phase phi matrices
       PllPoints = llPoints-1; % trims number of length points for solid phase phi matrices
       var_phi1 = 0; % sets end value of phi1 to be 0
   case 'ConstantVoltage'
       RlsPoints = lsPoints-1; % trims number of length points for solid phase rho matrices
       PlsPoints = lsPoints-1; % trims number of length points for solid phase phi matrices
       PllPoints = llPoints-1; % trims number of length points for solid phase phi matrices
       var_phi1 = 0; % sets end value of phi1 to be 0
   otherwise
     disp('ERROR: Specified connection type is unknown.')
     disp('Cannot correctly trim matrices.')
     error('An error occured. Simulation halted. See the above messages for information.')
end

err = 1; % sets initial error (dummy value to allow computation to start)
tol = 1e-27; % sets acceptable tolerance limit
j=1; % sets up iteration counter
t=0; % sets up time counting

%% Display information in command window

disp('Time step [s] is:')
disp(tStep)
disp('Dimensionless Lambda is:')
disp(Lambda)
disp('Maximum dimensionless varkappa is:')
disp(max(varkappa_vec))
disp('Dimensionless Beta is:')
disp(beta)

if Lambda > 1
    disp('Lambda stability criterion not satisfied.')
    disp('Reduce time step [s] to below')
    disp(tStep_Critical)
else
    disp('Lambda stability criterion satisfied.')
end 

disp('UltraCAP initialisation complete.')

%% Main computation

% CALL FUNCTIONS TO GENERATE TIME-INVARIANT MATRIX AND VECTOR FUNCTIONS

[P_1_inv] = P1MatrixFunction(Lambda,PlsPoints);
[P_2] = P2MatrixFunction(PlsPoints);
[C_1,C_1_var] = C1MatrixFunction(beta,PlsPoints,PllPoints);

while V_t <= v_terminal % DUMMY FOR NOW
    
% DETERMINE VALUE OF LAMBDA AT THIS TIME STEP

lambda = rho_net_vec_Current(end);
lambda = 0; % TEST VALUE

% CALL FUNCTIONS TO GENERATE TIME-VARAIANT MATRIX AND VECTOR FUNCTIONS

[R_1] = R1MatrixFunction(Sigma_vec,Theta,e,sigma,k_B,T,j_in,lsStep,RlsPoints,connection_type);
[A_1] = A1VectorFunction(Sigma_vec,var_rho,Upsilon,Theta,RlsPoints,connection_type);
% [B_1] = B1MatrixFunction(varkappa_vec,PlsPoints);

% SOLVE FOR SOLID PHASE ELECTRIC CHARGE DENSITY AT NEXT TIME LEVEL

[rho_vec_Next,rho_net_vec_Next] = SolveRho(rho_vec_Current,R_1,A_1,var_rho,rho_t,connection_type);

% RECORD AND SAVE TERMINAL VOLTAGE, CHARGE AND CURRENT INFORMATION

if j == 1
    I_vec = [];
    j_vec = [];
    Q_vec = [];
    rho_edge_vec = [];
end 

j_terminal = D_prime./rho_vec_Next(1)/lsStep*(- 1/2*rho_vec_Next(3) + 2*rho_vec_Next(2) - 3/2*rho_vec_Next(1)); % records estimate of specific current at current collector 
I_terminal = j_terminal*A_xc; % records estimate of current at current collector
rho_terminal = rho_net_vec_Next(1); % records net charge density at current collector

j_vec = horzcat(j_vec,j_terminal);
I_vec = horzcat(I_vec,I_terminal);
rho_edge_vec = horzcat(rho_edge_vec,rho_terminal);

% CALL FUNCTIONS TO GENERATE TIME-VARAIANT A2 VECTOR

[A_2] = A2VectorFunction(beth,daleth,Lambda,PlsPoints);

% SOLVE FOR SOLID PHASE ELECTRIC POTENTIAL AT NEXT TIME LEVEL

[phi_s_vec_Next] = SolvePhiS(rho_net_vec_Next,K);

% RECORD AND SAVE TERMINAL VOLTAGE AND CURRENT INFORMATION

if j == 1
    v_vec = [];
end 

v_terminal = phi_s_vec_Current(1); % records potential at current collector 

v_vec = horzcat(v_vec,v_terminal);

% SOLVE FOR LIQUID PHASE ELECTRIC POTENTIAL AT NEXT TIME LEVEL

[phi_l_vec_Next] = SolvePhiL(B_1,phi_l_vec_Current,phi_s_vec_Current,phi_s_vec_Next,lsEnd,llEnd,llVec,PlsPoints,PllPoints);

% SOLVE FOR LIQUID PHASE IONIC CONCENTRATION AT NEXT TIME LEVEL

[c_vec_Next] = SolveC(C_1,C_1_var,alpha,beta,c_vec_Current,c_0,phi_s_vec_Current,phi_s_vec_Next,phi_l_vec_Current,phi_l_vec_Next,PlsPoints,PllPoints);

% SAVE OUTPUT DATA AT REGULAR TIME INTERVALS

tMatch = find(j == tEvolutionPoints);

if isempty(tMatch) == 0

phi_s_mat(tMatch,:) = phi_s_vec_Next;
rho_mat(tMatch,:) = rho_net_vec_Next;
ls_mat(tMatch,:) = lsVec;
phi_l_mat(tMatch,:) = phi_l_vec_Next;
c_mat(tMatch,:) = c_vec_Next;
ll_mat(tMatch,:) = llVec;

end

j = j + 1; % update step count
t = t + tStep; % update time count
err = rho_vec_Current(1) - rho_vec_Next(1); % provisional error calculation
rho_vec_Current = rho_vec_Next;
rho_net_vec_Current = rho_net_vec_Next;
phi_s_vec_Current = phi_s_vec_Next;
phi_l_vec_Current = phi_l_vec_Next;
c_vec_Current = c_vec_Next;
D_vec_Current = D_prime./rho_vec_Current;
Sigma_vec = D_vec_Current*Gamma;

if t > tEnd
        disp('Simulation time limit reached.')
        break
end

end

%% Test output data

tVec = tVec(1:length(I_vec)); % trims time vector to correct length

figure(9)
plot(tVec,I_vec,'k');
title('Terminal current')
xlabel('Time ,s') 
ylabel('Current, A')

vVecAbs = abs(v_vec);

figure(10)
plot(tVec,vVecAbs,'k');
title('Terminal potential')
xlabel('Time, s') 
ylabel('Potential, V')

save('INSERTFILENAME.mat','variable_name_to_save')