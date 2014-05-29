% script.m
% Runs a single optimization using optim.m

% This is an example of how to use optim.m. 
% It runs an optimization using walk data from gaitdata_template.xls

clear all

% define the movement data to be tracked: Excel filename and sheet name
problem.gaitdata = {'gaitdata_template.xls' 'Slow Walk'};
problem.initialguess = 'result.mat';

% general optimization settings
problem.MaxIterations 	= 50000;
problem.print = 1;
problem.solver = 'IPOPT';		% use this solver

%=================================NM======================================
problem.N = 50;					% use this number of nodes (optim will reduce this to number of gait samples, if needed)
%=================================NM======================================

% lower and upper bounds on some of the variables we are optimizing (set lower=upper if variable should not be optimized)
problem.L_u1 = 0.001;		problem.U_u1 = 1.0; 		% uncomment this line to keep valve 1 closed
problem.L_u2 = 0.001;		problem.U_u2 = 1.0; 		% valve 2 control
problem.L_k  = 0.1;			problem.U_k  = 12;			% spring stiffness in MPa/ml

% weights for cost function
problem.w1 				= 1;			% weight for angle tracking term
problem.w2				= 1;			% weight for moment tracking term
problem.w3 				= .1;			% weight for valve 1 control accelerations
problem.w4 				= .1;			% weight for valve 2 control accelerations
problem.wreg 			= 0;			% weight for regularization (smoothness) term
problem.Wdisc 			= 0;			% weighting to encourage discrete control levels
problem.Ndisc 			= 2;			% number of discrete control levels
problem.phisd = 5*pi/180;				% we want to be within x deg when tracking
problem.Msd = 5;						% we want to be within x Nm when tracking

% model parameters for the hydraulic hardware
problem.G 			= 0.137;		% units: cm^-3
problem.C1max		= 30;			% Cmax for valve 1
problem.C2max		= 30;			% Cmax for valve 2
problem.B1 			= 0.01;			% viscous drag coefficient for valve 1 & its tubing
problem.B2 			= 0.01;			% viscous drag coefficient for valve 2 & its tubing
	
% do the optimization	
result = optim(problem);
result
if (result.info ~= 0)
	disp('THE OPTIMIZATION DID NOT CONVERGE, USE RESULTS WITH CAUTION');
end
