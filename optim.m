
function [result] = optim(problem);

% Solve the optimal control problem for the rotary hydraulic knee
% For model and methods, see RotaryKneeModel.doc

	global model mass
	model.eval = 0;
	
	tic
	
	% the following were taken out of optimsettings, because they do not need to be user modifiable
	checkderivatives 		= 1;			% check the derivatives (1) or not (0)
	modifyinitialguess		= 'none';		% gaitdata: phi and M in initial guess are replaced by gait data
	FeasibilityTolerance 	= 1e-5;         
	OptimalityTolerance 	= 1e-4;
	model.plot 				= 1;			% plots on screen during optimization (1) or not (0)
	model.pause 			= 0;			% pause after every cost function evaluation (1) or not (0)

	% some copying and transformations of problem settings
	solver			= problem.solver;
	MaxIterations 	= problem.MaxIterations;   
	model.print		= problem.print;
	if isfield(problem,'N')     
		N = problem.N;
	else
		N = 75;		% default is huge number, will be reduced to number of gait data samples
	end
	if ~isfield(problem,'prescribe_kinematics')
		problem.prescribe_kinematics = 0;		% default is not to prescribe kinematics
	end
	gaitdata 		= problem.gaitdata;
	initialguess 	= problem.initialguess;
	L_u1 = problem.L_u2;	U_u2 = problem.U_u1;    % Lower/Upper bounds for pump signal
	L_u2 = problem.L_u3;	U_u3 = problem.U_u2;    % Lower/Upper bounds for valve U2
    L_u3 = problem.L_u1;    U_u1 = problem.U_u3;    % Lower/Upper bounds for valve U3
	L_k = problem.L_k;		U_k = problem.U_k;      % spring stiffness in MPa/ml
	L_CPA = problem.L_CPA;	U_CPA = problem.U_CPA;  % compliance of second accumulator (ml/MPa)
	model.w1 		= problem.w1;                   % weight for angle tracking term
	model.w2 		= problem.w2;                   % weight for moment tracking term
	model.w3 		= problem.w3;                   % weight for valve 1 control accelerations
	model.w4 		= problem.w4;                   % weight for valve 2 control accelerations
    model.w5        = problem.w5;                   % weight for the pump controller
	model.wreg 		= problem.wreg;                 % weight for regularization (smoothness) term
	model.Wdisc 	= problem.Wdisc;                % weighting to encourage discrete control levels
	model.Ndisc 	= problem.Ndisc;                % number of discrete control levels
	model.gait.phisd 	= problem.phisd;            % we want to be within x deg when tracking
	model.gait.Msd 		= problem.Msd;              % we want to be within x Nm when tracking
	model.C1maxsquared 	= problem.C1max^2;			% C1max from document, squared
	model.C2maxsquared 	= problem.C2max^2;			% C2max, squared
	model.G = problem.G;
	model.B1 = problem.B1;
	model.B2 = problem.B2;
	model.datafile = char(gaitdata(1,1));
	model.movement = char(gaitdata(1,2));
	
	% load and store gait data
	ndata = size(gaitdata,1);
	if (ndata > 1)
		error('Current version of optim2.m can not track more than one movement');
	end
	warning off MATLAB:xlsread:Mode
	for i=1:ndata
        % read one sheet (one movement) from an XLS file with new format
        data = readgaitdata(char(gaitdata(i,1)), char(gaitdata(i,2)));
        tdiff = diff(data.t);
        if (std(tdiff) > 1e-6)
            error('Time interval on gait data file is not constant');
        end
        model.gait.T = max(data.t) + mean(tdiff);
        model.gait.phi = data.phi_a;
        model.gait.M = data.M_a;
        model.gait.P = data.P_a;    %(NM-power for the ankle)
    end
    

	% resample the gait data into N time points, if needed
	if (N > size(model.gait.phi,1))							% do not do more than gait data
		disp('WARNING: N was decreased to number of data samples.');
		N = size(model.gait.phi,1);
	elseif (N < size(model.gait.phi,1))
		oldtimes = [data.t ; model.gait.T];							% add a time point at start of next cycle
		model.gait.phi = [model.gait.phi ; model.gait.phi(1)];		% and the corresponding angle
		model.gait.M   = [model.gait.M   ; model.gait.M(1)  ];		% and the corresponding moment
        model.gait.P   = [model.gait.P   ; model.gait.P(1)  ];      % corresponding ankle power
		newtimes = (0:N-1)'/N*model.gait.T;							% the resample times
		model.gait.phi = interp1(oldtimes, model.gait.phi, newtimes);
		model.gait.M   = interp1(oldtimes, model.gait.M  , newtimes);
        model.gait.P   = interp1(oldtimes, model.gait.P  , newtimes);
    end
	
  
	N = size(model.gait.phi,1);			% number of samples in gait data
	model.N = N;
	h = model.gait.T/N;					% time step from the gait data, will also be the time step for direct collocation
	model.h = h;
	fprintf('Gait data loaded: %d samples per gait cycle.\n',N); 
    
	
	% collocation grid and unknowns
	Nvarpernode = 8;			% number of unknowns per node: u2,u3,s,v1,v2,phi,M, u1
	model.Nconpernode = 4;		% number of constraint equations per node
	model.Nvar = model.N * Nvarpernode;			% total number of unknowns
	model.Ncon = model.N * model.Nconpernode;			% total number of constraints
	model.Nvarpernode = Nvarpernode;
    
	% precalculate some indices for X array, to speed up the calculations
	iu1 = zeros(N,3);		% index to u1 control at three successive nodes (pump)
	iu2 = zeros(N,3);		% index to u2 control at three successive nodes (valve U2)
    iu3 = zeros(N,3);       % index to u3 control at three successive nodes (valve U3)
	for i=1:N
		if (i == N)
			iu1(i,:) = [N-1 0 1]*Nvarpernode + 1;	% u1 is the first variable		
		elseif (i == N-1)
			iu1(i,:) = [N-2 N-1 0]*Nvarpernode + 1;		
		else
			iu1(i,:) = [i-1 i i+1]*Nvarpernode + 1;		
		end	
	end
	iu2 = iu1 + 1;							% u2 is stored immediately after u1
    iu3 = iu2 + 1;                          % u3 is stored immediately after u2  (NM-line added for storing u3)
	iM = (0:N-1)*Nvarpernode + 7;			% M is 7th variable at each node
	iphi = (0:N-1)*Nvarpernode + 6;			% phi is 6th variable
	model.iu1 = iu1;
	model.iu2 = iu2;
    model.iu3 = iu3;                        
	model.iu = [iu1(:,2) ; iu2(:,2) ; iu3(:,2)];		% simply a list of all controls within X
	model.iphi = iphi;
	model.iM = iM;
	model.iCPA = model.Nvar+1;      model.Nvar = model.Nvar+1;
	model.ik =   model.Nvar+1;      model.Nvar = model.Nvar+1;
	
	% precalculate the Hessian of objective function (since it is constant for our least squares objective)
	Htrack = spalloc(model.Nvar,model.Nvar,1);
	Hcontrol = spalloc(model.Nvar,model.Nvar,1);
	e = ones(N,1);
	% term 1: angle tracking
	Htrack(iphi, iphi) = model.w1 * 2 * spdiags(e,0,N,N)/model.gait.phisd^2/N;
	% term 2: moment tracking
	Htrack(iM, iM) = model.w2 * 2 * spdiags(e,0,N,N)/model.gait.Msd^2/N;
	% term 3: acceleration of valve u2 control
	% general pattern is Hu for finite difference accelerations with periodicity
	Hu = spdiags([2*e -8*e 12*e -8*e 2*e],-2:2,N,N)/N/h^2;
	Hu(1,N-1:N) = [2 -8]/N/h^2;
	Hu(2,N) = 2/N/h^2;
	Hu(N-1:N,1) = [2 -8]'/N/h^2;
	Hu(N,2) = 2/N/h^2;
	Hcontrol(iu2(:,1),iu2(:,1)) = model.w3*Hu;
	% term 4: acceleration of valve 2 control
	Hcontrol(iu3(:,1),iu3(:,1)) = model.w4*Hu;

	model.Htrack = Htrack;
	model.Hcontrol = Hcontrol;
	model.Hreg = model.wreg * 2*spdiags(ones(model.Nvar,1),0,model.Nvar,model.Nvar);
	Hnnz = nnz(Htrack + Hcontrol);
	fprintf('Hessian sparsity:  %d nonzero elements out of %d (%8.3f %%)\n', Hnnz, model.Nvar^2, Hnnz/model.Nvar^2);
		
	% set lower and upper bounds
	Lnode = [L_u1	L_u2 -100 -100 -100 -150*pi/180 -350  L_u3]';          
	Unode = [U_u1	U_u2  100  100  100  150*pi/180  350  U_u3]';  
	L = [repmat(Lnode,N,1)];
	U = [repmat(Unode,N,1)];
	L(model.ik) = L_k;
	U(model.ik) = U_k;
	L(model.iCPA) = L_CPA;				
	U(model.iCPA) = U_CPA;				
    
	if problem.prescribe_kinematics		% constrain kinematics to be equal to gait data
		L(model.iphi) = model.gait.phi;
		U(model.iphi) = model.gait.phi;
	end
	
	if strcmp(initialguess,'mid')
		X0 = L + 0.5*(U-L);
	elseif strcmp(initialguess,'random')
		X0 = L + (U-L).*rand(size(L));
	elseif strcmp(initialguess,'meas');
		X0 = (L+U)/2;						% start with midpoint for all unknowns
		X0(model.iphi) = model.gait.phi;	% insert measured joint angles
		X0(model.iM) = model.gait.M;		% insert measured joint moments
	elseif strcmp(initialguess,'zeros');
		X0 = zeros(size(L));
	else
		load(initialguess);
		X0 = result.X;
		
		% extract the time series of all 8 variables, and k and P0
		N0 = (size(X0,1)-2)/8;             
        if (N0 ~= round(N0))
			error('N0 is not a whole number');
		end
		x0 = reshape(X0(1:end-2),8, N0)';
		P0 = X0(end-1);
		k = X0(end);
		
		% add one node so we can interpolate with periodicity
		x0 = [x0 ; x0(1,:)];
		
		% interpolate to the current number of nodes
		x0 = interp1((0:N0)'/N0,x0,(0:N-1)'/N,'linear','extrap');
		X0 = reshape(x0',8*N,1);    
		X0 = [X0 ; P0 ; k];

	end

	if numel(strfind(modifyinitialguess,'gaitdata')) > 0
		% replace the phi and M unknowns with the corresponding gait data
		X0(model.iphi) = model.gait.phi;
		X0(model.iM) = model.gait.M;
	end
	if numel(strfind(modifyinitialguess,'randomize')) > 0
		% add small random numbers to initial guess
		X0 = X0 + 0.001*randn(size(X0));
	end

	
	% find the Jacobian pattern
	X = L + (U-L).*rand(size(L));		% a random vector of unknowns
    model.Jnnz = 50*model.N;            % first guess of Jacobian size
	J = conjac(X);
	model.Jnnz = nnz(J);                % actual Jacobian size
	fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n', ...
		model.Jnnz, model.Ncon*model.Nvar, 100*model.Jnnz/(model.Ncon*model.Nvar));
	model.Jpattern = double(J~=0);
	
	% check the derivatives at initial guess X0
	if (checkderivatives)
		model.FDvar = 1;
		hh = 1e-7;
		X = X0;
		f = objfun(X);
		grad = objgrad(X);
		hess = objhess(X);
		c = confun(X);
		cjac = conjac(X);
		cjac_num = zeros(model.Ncon,model.Nvar);
		grad_num = zeros(model.Nvar,1);
		hess_num = zeros(model.Nvar,model.Nvar);
		for i=1:model.Nvar
			fprintf('checking derivatives for unknown %4d of %4d\n',i,model.Nvar);
			Xisave = X(i);
			X(i) = X(i) + hh;
			cjac_num(:,i) = (confun(X) - c)/hh;
			grad_num(i) =   (objfun(X) - f)/hh;
			hess_num(:,i) = (objgrad(X) - grad)/hh; 
			X(i) = Xisave;
		end

		% find the max difference in constraint jacobian and objective gradient
		[maxerr,irow] = max(abs(cjac-cjac_num));
		[maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d\n', maxerr, irow(icol), icol);
        fprintf('    Conjac result: %f\n', full(cjac(irow(icol), icol)));
        fprintf('    Fin Diff appr: %f\n', full(cjac_num(irow(icol), icol)));
		d = 2*abs(cjac - cjac_num)./(cjac + cjac_num);		% relative error
		fprintf('Max. relative error in constraint jacobian: %8.5f\n', max(max(d)));
		[maxerr,irow] = max(abs(grad_num-grad));
		fprintf('Max.error in objective gradient:  %8.5f at %d\n', maxerr, irow);
		[maxerr,irow] = max(abs(hess-hess_num));
		[maxerr,icol] = max(maxerr);
		fprintf('Max.error in objective Hessian:   %8.5f at %d %d\n', maxerr, irow(icol), icol);
		d = 2*abs(hess - hess_num)./(hess + hess_num);		% relative error
		fprintf('Max. relative error in objective Hessian: %8.5f\n', max(max(d)));
		
		% find errors in jacobian pattern, are there nonzeros in actual Jacobian that are not in pattern?
		Jpattern_num = sparse(cjac_num ~= 0);
		Jdiff = Jpattern_num - model.Jpattern;				% any +1 value in the diff matrix indicates an error
		[ierr,jerr] = find(Jdiff>0);
		fprintf('Errors in constraint jacobian pattern: \n');
		[ierr jerr]
		model.FDvar = 0;
		
	end

	% report something about the initial guess, unless we're not even optimizing
	model.FDvar = 1;
	if ~strcmp(solver,'none')
		fprintf('Initial guess cost function value: %8.5f\n',objfun(X0));
		model.normc = norm(confun(X0));
		fprintf('Initial guess constr. violation:   %8.5f\n',model.normc);
		report(X0);
		model.FDvar = 0;
		if (model.print)
			disp('Hit ENTER to start optimization');
			pause
		end
	end
		
	starttime = cputime;
	if strcmp(solver, 'SNOPT')
		% do the optimization using SNOPT
		Prob = conAssign(@objfun, @objgrad, @objhess, [], L, U, 'Rotary Hydraulic Knee', X0, ...
					[], 0, ...
					[], [], [], @confun, @conjacSNOPT, [], model.Jpattern, ...
					zeros(model.Ncon,1), zeros(model.Ncon,1), ...
					[], [], [],[]);
		Prob.SOL.optPar(1)= 11;			% uncomment this to get snoptsum.txt output file
		Prob.SOL.optPar(9) = FeasibilityTolerance;
		Prob.SOL.optPar(10) = OptimalityTolerance;
		Prob.SOL.optPar(11) = 1e-6; % Minor feasibility tolerance (1e-6)
		Prob.SOL.optPar(30) = 1000000; % maximal sum of minor iterations (max(10000,20*m))
		Prob.SOL.optPar(35) = MaxIterations;
		Prob.SOL.optPar(36) = 40000; % maximal number of minor iterations in the solution of the QP problem (500)
		Prob.SOL.PrintFile = 0;
		Prob.SOL.SummFile = 0;
		Prob.SOL.moremem = 10000000; % increase internal memory
		
		Result = tomRun('snopt',Prob);
		X = Result.x_k;
		disp('--------------------------------------');
		disp('SNOPT finished.')
		disp(Result.ExitText);
		fprintf('Number of iterations: %d\n',Result.Iter);
		info = Result.ExitFlag;
        
	elseif strcmp(solver,'IPOPT')
		% do the optimization using IPOPT
        options.lb = L;
        options.ub = U;
        options.cl = zeros(model.Ncon,1);
        options.cu = zeros(model.Ncon,1);
        funcs.objective = @objfun;
        funcs.gradient = @objgrad;
        funcs.constraints = @confun;
        funcs.jacobian = @conjac;
        funcs.jacobianstructure = @conjacstructure;
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.print_level           = 0;       
            options.ipopt.max_iter = MaxIterations;
            options.ipopt.tol = OptimalityTolerance;
            options.ipopt.acceptable_constr_viol_tol = FeasibilityTolerance;
        
		[X, info] = ipopt(X0, funcs, options);
        
		disp('--------------------------------------');
		disp('IPOPT finished.');
        info = info.status;
		disp(IPOPTstatus(info));
	else
		disp('Solver name not recognized, reproducing initial guess.');
		X = X0;
		info = 0;
    end
	
    
    % display the results
    disp(['Total time used: ' num2str(cputime-starttime) ' seconds.']);
    disp('')
    report(X,1);
    disp('')
    % display model parameters
    fprintf('Optimal parameter values:\n');
    fprintf('    k   = %8.4f MPa/ml       (stiffness of spring loaded reservoir)\n', X(model.ik));
    fprintf('    CPA = %8.4f ml/MPa       (compliance of second accumulator)\n', X(model.iCPA));
    disp('')
    
    % If using discritization with pump for continuous flow, print the value of this signal
    if model.Wdisc ~0
    fprintf('    U1  = %8.4f ml/sec       (flow from the pump)\n', X(model.iu1));
    end

%     % Choosing to save the results on file
%     saveresult = input('Would you like to save the results? (y/n):   ','s');
%     if strcmpi(saveresult, 'y');
%         file = input('Enter filename:  ','s');
%         disp('')
%         clear Result result 
%         result.info = info;
%         result.X = X;
%         result.N = N;			% in case it was reduced to number of data samples
%         result.RMSang = 180/pi*sqrt( mean( (X(model.iphi) - model.gait.phi).^2 ) );
%         result.RMSmom = sqrt( mean( (X(model.iM)   - model.gait.M).^2   ) );
%         result.costfun = objfun(X);
%         result.k = X(end);
%         result.w5 = model.w5;
%         save(file,'result');
%         disp('')
%         disp('Result of optimization is saved')
%         disp('Please change initial guess in readgaitdata.m to optimize these results')
%     elseif strcmpi(saveresult, 'n');
%         return
%     elseif isempty(saveresult)
%         disp('Saving results...')
%         disp('')
%         clear Result result 
%         result.info = info;
%         result.X = X;
%         result.N = N;			% in case it was reduced to number of data samples
%         result.RMSang = 180/pi*sqrt( mean( (X(model.iphi) - model.gait.phi).^2 ) );
%         result.RMSmom = 		sqrt( mean( (X(model.iM)   - model.gait.M).^2   ) );
%         result.costfun = objfun(X);
%         result.k = X(end);
%         result.w5 = model.w5;
%         save('resultX.mat','result');
%         disp('')
%         disp('Result of optimization is saved in resultX.mat.')
%     else
%         return
%     end

    
    % Autosaving instead of prompt save. This will save file anyway
    disp('Saving results...')
    disp('')
    clear Result result 
    result.info = info;
    result.X = X;
    result.N = N;			% in case it was reduced to number of data samples
    result.RMSang = 180/pi*sqrt( mean( (X(model.iphi) - model.gait.phi).^2 ) );
    result.RMSmom = 		sqrt( mean( (X(model.iM)   - model.gait.M).^2   ) );
    result.costfun = objfun(X);
    result.k = X(end);
    result.w5 = model.w5;
    save('resultX.mat','result');
    disp('')
    disp('Result of optimization is saved in resultX.mat.')
    
    
    
end
%===============================================================================
function report(X, powerreport)

	global model 
    
	phi = X(model.iphi)*180/pi;
	M =   X(model.iM);
	u1 =  X(model.iu1(:,2));
	u2 =  X(model.iu2(:,2));
    u3 =  X(model.iu3(:,2));    
	P1 = X(model.ik) * [X(3:model.Nvarpernode:end) ; X(3)]'; % Pressure in the accumulator
	v1 = [X(4:model.Nvarpernode:end) ; X(4)]';
	v2 = [X(5:model.Nvarpernode:end) ; X(5)]';
	tperc = 100*(0:model.N)'/model.N;
	phi = [phi ; phi(1)];
	M =   [M ; M(1)]';
	u1 =  [u1 ; u1(1)]';
	u2 =  [u2 ; u2(1)]';
    u3 =  [u3 ; u3(1)]';    
	P = -M * model.G; %Pressure in actuator
	ppump = u1.*P1;
	gaitphi = model.gait.phi*180/pi;
	gaitM = model.gait.M;
    gaitP = (model.gait.P);    
	gaitphi = [gaitphi ; gaitphi(1)];
	gaitM = [gaitM ; gaitM(1)];
    gaitP = [gaitP ; gaitP(1)];
	figure(1);
	clf;
	
	subplot(3,3,1)  
	plot(tperc,-gaitphi,'b',tperc,-phi,'r');
	ylabel('Ankle Angle (deg)');
	title(model.datafile);
	
	subplot(3,3,4)  
	plot(tperc,gaitM,'b',tperc,M,'r');
	ylabel('Ankle moment (Nm)');
	legend('desired','prosthesis');
	xlabel('Time (% of cycle)');
	
	subplot(3,3,2)  
	plot(tperc,u2);
	ylabel('Valve U2 control (a.u.)');
	title(model.movement);
	
	subplot(3,3,5)  
	plot(tperc,u3);
	ylabel('Valve U3 control (a.u.)');
	xlabel('Time (% of cycle)');

	subplot(3,3,3); 
	plot(tperc, P, tperc, P1);
	ylabel('pressure (MPa)');
	legend('P','P1');
	title(['N = ' num2str(model.N)]);

	subplot(3,3,6); 
	plot(tperc,v1,tperc,v2, tperc,u1);
	xlabel('Time (% of cycle)');
	ylabel('flow (ml/s)');
	legend('v1','v2','U1');
    
    
    subplot(3,3,7); 
    plot(tperc,ppump,'r',tperc,gaitP,'b');
    xlabel('Time (% of cycle)');
    ylabel('Power (Watts)');
    legend('Pump','Desired');
    
	
	if (nargin > 1 && powerreport == 1)
		figure(2)
		pspring = P1.*v1;
		pvalve1 = (P-P1).*v1;
		pvalve2 = P.*v2;
        ppump = u1.*P1;
		ptotal = pspring + pvalve1 + pvalve2 + ppump;
		plot(tperc,ptotal,tperc,pspring,tperc,pvalve1,tperc,pvalve2,tperc,ppump);
		xlabel('Time (% of gait cycle');
		ylabel('Power (W)');
		legend('total','spring','valve1','valve2','pump');
	end
	
end
%====================================================================
function [c] = confun(X)
%Becuase there are constraints on our variables, it is necessary to utilize
%this otherwise optional funcion

	global model
	
	h = model.h;
	iCPA = model.iCPA;
	ik = model.ik;
 	
	c = zeros(model.Ncon,1);
	irow = 1;
	ix1 = 1:model.Nvarpernode;
	for i=1:model.N
		% extract variables from successive nodes
		x1 = X(ix1);
		if (i < model.N)
			ix2 = ix1 + model.Nvarpernode;
		else
			ix2 = 1:model.Nvarpernode;
		end
		x2 = X(ix2);
		
		% generate the four constraints
		% the eight variables are: u1,u2,s,v1,v2,phi,M,u3
		
		% ds/dt + v1 + u1 = 0:
		c(irow) = (x2(3) - x1(3))/h  + (x1(4) + x2(4))/2.0  +  (x1(1) + x2(1))/2;   
		
		% u2^2 * C1max * (k s + M G - B1 v1) - v1 * |v1| = 0
		c(irow+1) = x1(2)^2*model.C1maxsquared * (X(ik) * x1(3) + x1(7) * model.G - model.B1 * x1(4) ) - x1(4)*abs(x1(4));
		
		% dphi/dt - G*(v1 + v2 - G*dM/dt*CPA) = 0
		c(irow+2) = (x2(6)-x1(6))/h - model.G*(x1(4)+x2(4) + x1(5)+x2(5))/2.0 - model.G^2 * X(iCPA) * (x2(7)-x1(7))/h;
        
		% u3^2 * C2max * (M * G - B2 * v2) - v2 * |v2|  = 0
		c(irow+3) = x1(8)^2 * model.C2maxsquared * (x1(7) * model.G - model.B2 * x1(5) ) - x1(5)*abs(x1(5)); 
        
        
				
		%  advance ix1 and irow to next node
		ix1 = ix1 + model.Nvarpernode;
		irow = irow + model.Nconpernode;
	end
	
	if model.FDvar
		return
	end
	
	model.normc = norm(c);

end
%==========================================================================================
function J = conjacSNOPT(X);
	% returns constraint Jacobian matrix, for SNOPT
	J = conjac(X,0);
end
%==========================================================================================
function [Jstruct] = conjacstructure(X)
    global model
	Jstruct = model.Jpattern;
end
%==========================================================================================
function [J] = conjac(X)
	global model

	h = model.h;		% time step size
	iCPA = model.iCPA;
	ik = model.ik;

	J = spalloc(model.Ncon,model.Nvar, model.Jnnz);		% allocate memory space
	irow = 1;
	ix1 = 1:model.Nvarpernode;
	for i=1:model.N
		% extract variables from successive nodes
		x1 = X(ix1);
		if (i < model.N)
			ix2 = ix1 + model.Nvarpernode;      % use the x variables from the node after ix1
		else
			ix2 = 1:model.Nvarpernode;          % use the x variables from node 1
		end
		x2 = X(ix2);
		
		% generate the constraint derivatives
		% the variables are: u1,u2,s,v1,v2,phi,M,u3
		
		% ds/dt + v1 + u1 = 0:
		% c(irow) = (x2(3) - x1(3))/h  + (x1(4) + x2(4))/2.0  +  (x1(1) + x2(1))/2;
		J(irow,ix1(3)) = -1/h;
		J(irow,ix1(4)) = 0.5;
        J(irow,ix1(1)) = 0.5;
		J(irow,ix2(3)) = 1/h;
		J(irow,ix2(4)) = 0.5;
        J(irow,ix2(1)) = 0.5;
        
        
		% u2^2 * C1max * (k s + M G - B1 v1) - v1 * |v1| = 0
		% c(irow+1) = x1(1)^2*model.C1max * (X(ik) * x1(3) + x1(7) * model.G - model.B1 * x1(4) ) - x1(4)*abs(x1(4));
		J(irow+1,ix1(2)) = 2*x1(2)*model.C1maxsquared * (X(ik) * x1(3) + x1(7) * model.G - model.B1 * x1(4) );
		J(irow+1,ix1(3)) = x1(2)^2*model.C1maxsquared * X(ik);
		J(irow+1,ix1(4)) = -x1(2)^2*model.C1maxsquared * model.B1 - 2*abs(x1(4));
		J(irow+1,ix1(7)) = x1(2)^2*model.C1maxsquared * model.G;
		J(irow+1,ik)     = x1(2)^2*model.C1maxsquared * x1(3);
		
      
        % dphi/dt - G*(v1 + v2 - G*dM/dt*CPA) = 0
		%c(irow+2) = (x2(6)-x1(6))/h - model.G*(x1(4)+x2(4) + x1(5)+x2(5))/2.0 - model.G^2 * X(iCPA) * (x2(7)-x1(7))/h;  
		J(irow+2,ix1(4)) = -model.G/2.0;
		J(irow+2,ix1(5)) = -model.G/2.0;
		J(irow+2,ix1(6)) = -1/h;
        J(irow+2,ix1(7)) = model.G^2*X(iCPA)/h;
		J(irow+2,ix2(4)) = -model.G/2.0;
		J(irow+2,ix2(5)) = -model.G/2.0;
		J(irow+2,ix2(6)) = 1/h;
        J(irow+2,ix2(7)) = -model.G^2*X(iCPA)/h;
        J(irow+2,iCPA)   = -model.G^2 * (x2(7)-x1(7))/h;
		
		% u3^2 * C2max * (P0 + M * G - B2 * v2) - v2 * |v2|  = 0
		%c(irow+3) = x1(2)^2 * model.C2max * (P0 + x1(7) * model.G - model.B2 * x1(5) ) - x1(5)*abs(x1(5);
		J(irow+3,ix1(8)) =  2*x1(8)*model.C2maxsquared * (x1(7) * model.G - model.B2 * x1(5) );
		J(irow+3,ix1(5)) = -x1(8)^2 * model.C2maxsquared * model.B2 - 2*abs(x1(5));
		J(irow+3,ix1(7)) =  x1(8)^2 * model.C2maxsquared * model.G;

		%  advance ix1 and irow to next node
		ix1 = ix1 + model.Nvarpernode;
		irow = irow + model.Nconpernode;
        
       
    end
end
%====================================================================
function [f] = objfun(X, Prob);
	% objective function for the optimization
	
	global model
	model.eval = model.eval+1;

	h = model.h;		% time step size
	
	% term 1: angle tracking
	f1 = model.w1 * mean( ( (X(model.iphi) - model.gait.phi)/model.gait.phisd ).^2 );
	
	% term 2: moment tracking
	f2 = model.w2 * mean( ( (X(model.iM) - model.gait.M)/model.gait.Msd ).^2 );
	
	% term 3: acceleration of valve U2 control
	f3 = model.w3 * mean( (X(model.iu2(:,1)) - 2*X(model.iu2(:,2)) + X(model.iu2(:,3)) ).^2 )/h^2;
    
	% term 4: acceleration of valve U3 control
	f4 = model.w4 * mean( (X(model.iu3(:,1)) - 2*X(model.iu3(:,2)) + X(model.iu3(:,3)) ).^2 )/h^2;
	
    % term 5: Square of the flow rate for pump U1
    f5 = model.w5 * mean( (X(model.iu1(:,2))) ).^2;
    
	% regularization term
	freg = model.wreg * sum(X.^2);
	
	% discretization term
	if model.Wdisc ~= 0
		fdisc = model.Wdisc * sum(sin(pi*(model.Ndisc-1)*X(model.iu)).^2);  % Should there be an integer by the third model element?
    else
		fdisc = 0;
    end
    
	% add them up
	f = f1 + f2 + f3 + f4 + f5 + freg + fdisc;

	if model.FDvar
		return
	end

	if (toc < 1.0)
		return
	end
	tic;
	
	if (model.print)
        disp('')
		fprintf('%d: Normc: %9.5f   ', model.eval, model.normc);
		fprintf('Objfun: %8.4f=%8.4f(ang)+%8.4f(mom)+%8.4f(u1)+%8.4f(u2)+%8.4f(u3)+%8.4f(reg)+%8.4f(dis)\n', f,f1,f2,f3,f4,f5,freg, fdisc);
        disp('')
	end

	if (model.plot)
		report(X);
		drawnow;
	end
	if (model.pause)
		pause
	end
	
end
%====================================================================
function [g] = objgrad(X);
	% gradient of objective function

	global model
	h = model.h;
	g = zeros(model.Nvar,1);
	
	% term 1: angle tracking
	g(model.iphi) = model.w1 * 2 * (X(model.iphi) - model.gait.phi)/model.gait.phisd^2/model.N;

	% term 2: moment tracking
	g(model.iM) = model.w2 * 2 * (X(model.iM) - model.gait.M)/model.gait.Msd^2/model.N;
	
	% term 3: acceleration of valve U2 control
	g(model.iu2(:,1)) = model.w3 * 2 * (X(model.iu2(:,1)) - 2*X(model.iu2(:,2)) + X(model.iu2(:,3)) )/model.N/h^2;
	g(model.iu2(:,2)) = g(model.iu2(:,2)) - model.w3 * 4 * (X(model.iu2(:,1)) - 2*X(model.iu2(:,2)) + X(model.iu2(:,3)) )/model.N/h^2;
	g(model.iu2(:,3)) = g(model.iu2(:,3)) + model.w3 * 2 * (X(model.iu2(:,1)) - 2*X(model.iu2(:,2)) + X(model.iu2(:,3)) )/model.N/h^2;
	
	% term 4: acceleration of valve U3 control
	g(model.iu3(:,1)) = model.w4 * 2 * (X(model.iu3(:,1)) - 2*X(model.iu3(:,2)) + X(model.iu3(:,3)) )/model.N/h^2;
	g(model.iu3(:,2)) = g(model.iu3(:,2)) - model.w4 * 4 * (X(model.iu3(:,1)) - 2*X(model.iu3(:,2)) + X(model.iu3(:,3)) )/model.N/h^2;
	g(model.iu3(:,3)) = g(model.iu3(:,3)) + model.w4 * 2 * (X(model.iu3(:,1)) - 2*X(model.iu3(:,2)) + X(model.iu3(:,3)) )/model.N/h^2;

	% regularization term
	g = g + model.wreg*2*X;
	
	% discretization term
	if model.Wdisc ~= 0
		a = pi*(model.Ndisc-1)*X(model.iu);
		g(model.iu) = g(model.iu) + model.Wdisc * pi * (model.Ndisc-1)*sin(2*a);
	end
	

end
%====================================================================
function [H] = objhess(X);
	% hessian of objective function
	
	global model
	
	% hessian is constant, so we just pull it out of the model struct
	H = model.Htrack + model.Hcontrol + model.Hreg;
	
	% discretization term Hessian is not constant:
	if model.Wdisc ~= 0
		a = pi*(model.Ndisc-1)*X(model.iu);
		Hdiag = zeros(model.Nvar,1);
		Hdiag(model.iu) = model.Wdisc * 2 * pi^2 * (model.Ndisc-1)^2 * cos(2*a);
		H = H + spdiags(Hdiag,0,model.Nvar,model.Nvar);
	end

end
%=====================================================================
function [s] = IPOPTstatus(code);
	% translates IPOPT status code to a string
	
	if code==0
		s = 'solved';
	elseif code==1
		s = 'solved to acceptable level';
	elseif code==2
		s = 'infeasible problem detected';
	elseif code==3
		s = 'search direction becomes too small';
	elseif code==4
		s = 'diverging iterates';
	elseif code==5
		s = 'user requested stop';
	elseif code==-1
		s = 'IPOPT failed: maximum number of iterations exceeded';
	elseif code==-2
		s = 'IPOPT failed: restoration phase failed';
	elseif code==-3
		s = 'IPOPT failed: error in step computation';
	elseif code==-10
		s = 'IPOPT failed: not enough degrees of freedom';
	elseif code==-11
		s = 'IPOPT failed: invalid problem definition';
	elseif code==-12
		s = 'IPOPT failed: invalid option';
	elseif code==-13
		s = 'IPOPT failed: invalid number detected';

	elseif code==-100
		s = 'IPOPT failed: unrecoverable exception'; % divide by zero
	elseif code==-101
		s = 'IPOPT failed: non-IPOPT exception thrown';
	elseif code==-102
		s = 'IPOPT failed: insufficient memory';
	elseif code==-199
		s = 'IPOPT failed: internal error';
	else
		s = ['IPOPT status code ' num2str(code) ' not recognized'];
	end
	
end
%=====================================================================
function iterfunc(T,F);
	% function called by IPOPT after each iteration
	global model
	% model.eval = T;
end
