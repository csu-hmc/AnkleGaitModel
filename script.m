%opten optim
% script.m
% Runs a single optimization using optim.m

% This is an example of how to use optim.m. 
% It runs an optimization using walk data from Gait_Data_Sub1

clear all

% define the movement data to be tracked: Excel filename and sheet nasme
problem.gaitdata = {'Gait_Data_Sub1.xls' 'Slow Walk'};

% Uncomment one of these structure elements so the program will use either
% a random initial guess, or results from the prevoius optimization. Only
% use result.mat if the program converged and the results were okay.
problem.initialguess = 'random';
%problem.initialguess = 'result.mat';

% general optimization settings
problem.MaxIterations 	= 10000;
problem.print = 1;
problem.solver = 'IPOPT';		% use this solver

problem.N = 50;					% use this number of nodes (optim will reduce this to number of gait samples, if needed

% lower and upper bounds on some of the variables we are optimizing (set lower=upper if variable should not be optimized)
problem.L_u1 = 0.001;		problem.U_u1 = 1.0; 		% uncomment this line to keep valve 1 closed
problem.L_u2 = 0.001;		problem.U_u2 = 1.0; 		% valve 2 control
problem.L_k  = 0.1;			problem.U_k  = 12;			% spring stiffness in MPa/ml
problem.L_u3 = -50.0;       problem.U_u3 = 50.0;        %(NM-added u3 lower and upper bounds to the problem statement)

% weights for cost function
problem.w1 				= 10.0;			% weight for angle tracking term
problem.w2				= 1.0;			% weight for moment tracking term
problem.w3 				= 0.1;			% weight for valve 1 control accelerations
problem.w4 				= 0.1;			% weight for valve 2 control accelerations
problem.w5              = 0.05;          % weight for the pump controller (NM)
problem.wreg 			= 0.0;			% weight for regularization (smoo thness) term
problem.Wdisc 			= 0.0;			% weighting to encourage discrete control levels
problem.Ndisc 			= 2.0;			% number of discrete control levels (dont touch)
problem.phisd = 5*pi/180;				% we want to be within x deg when tracking
problem.Msd = 5;						% we want to be within x Nm when tracking

% model parameters for the hydraulic hardware
problem.G 			= 0.137;		% units: cm^-3
problem.C1max		= 30;			% Cmax for valve 1
problem.C2max		= 30;			% Cmax for valve 2
problem.B1 			= 0.01;			% viscous drag coefficient for valve 1 & its tubing
problem.B2 			= 0.01;			% viscous drag coefficient for valve 2 & its tubing
problem.B3          = 0.01;         % viscous drag coefficient for pump and hardware
problem.C           = 0;            % coefficient for parallel accumulator
	
% do the optimization	
result = optim(problem);
result
if (result.info ~= 0)
	disp('THE OPTIMIZATION DID NOT CONVERGE, USE RESULTS WITH CAUTION');
end





