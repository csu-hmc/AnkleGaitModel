
function [data] = readgaitdata(filename, movement)
% Reads Orthotrak gait data from XLS file and performs the necessary conversions.

% Revision history
% 2010-08-18	Extracting data.F_rot_ankle from Excel file (Rick Rarick)
% 2011-01-18	Corrected column numbers for data files sent by Sergey on 1/17
%				NOTE only corrected those needed for optim program, check others too!
%				Added abs() to detecting gait from hip motion, so also works when subject walks with decreasing X
% 2014-05-22    Modifying code to output values for hydraulic prosthetic
%               ankle system. (Nick Mavros)

% Output is a struct containing time series data:
%	t				time in seconds
%	t_perc			time in percentage of the total movement time (or gait cycle)
%   F_rot_ankle     ankle force rot. (N)
%	phi_1			thigh angle
%	phi_knee		knee angle
%	x_h, y_h		hip position
%	x_k, y_k		knee position
%	x_a, y_a		ankle position
%	Fx_h, Fy_h		force applied to top of thigh by other leg (global reference frame)
%	Fx_k, Fy_k		force applied to knee by shank (shank reference frame)
%	M_h				hip moment
%   M_k				knee moment

% These variables are all defined in Prosthesis_Geometry_2010_06_01.doc (except knee and ankle position)

% Output struct also contains scalar data
%	subjectID		subject ID number
%	mass			subject mass (kg)
%	cadence			cadence (steps/min), if present on file, otherwise NaN.
%	cycle			total movement time (s)
%	L1				mean distance from hip to knee
%	L2				mean distance from knee to ankle

	% constants
	g = 9.81;			% gravity

	% first read the subject data sheet
	d = xlsread(filename, 'Subject', '', 'basic');
	data.subjectID = d(3,1);
	data.mass = d(4,1);

	% now read the movement data
	d = xlsread(filename, movement, '', 'basic');
	keyboard
	data.cadence = d(1,2);
	d = d(5:end,:);			% remove the four header lines

	% time
	data.t_perc = d(:,2);
	data.t = d(:,3);
	if (data.t_perc(1) ~= 0)
		error('Time(%) does not start at zero');
	elseif (data.t_perc(end) == 100)
		error('Time(%) ends at 100');
	elseif (std(diff(data.t_perc)) > 0.001)
		error('Time(%) is not equally spaced.');
	end
	data.cycle = max(data.t) + data.t(end) - data.t(end-1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Added by Rick Rarick 8/18/2010
    
    data.F_rot_ankle = d(:, 9)*data.mass*g;  % convert to N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% hip, knee, ankle XY position data, convert from mm to m
	data.x_h = d(:,25)/1000;
	data.y_h = d(:,27)/1000;
	data.x_k = d(:,28)/1000;
	data.y_k = d(:,30)/1000;
	data.x_a = d(:,31)/1000;
	data.y_a = d(:,33)/1000;
	
	% if there was more than 0.5 m of horizontal motion, we assume it was gait
	if (abs(data.x_h(end) - data.x_h(1)) > 0.5)
		fprintf('readgaitdata: it looks like %s is gait.\n', movement);
		fprintf('it will be assumed that gait is symmetrical with 50%% phase shift between legs.\n');
		gait = 1;
	else
		fprintf('readgaitdata: it looks like %s is not gait.\n', movement);
		fprintf('it will be assumed that movement is symmetrical with 0%% phase shift between legs.\n');
		gait = 0;
    end
	
	% calculate thigh, shank, and ankle angles
	data.phi_1 = atan2(data.x_k - data.x_h, data.y_h - data.y_k);   %Thigh
	data.phi_2 = atan2(data.x_a - data.x_k, data.y_k - data.y_a);   %Shank
    
    %=================================NM====================================
	% joint angles and moments
	data.phi_k = -d(:,4)*pi/180;				% convert to ankle angle in rad, negative when flexed
	data.M_h = d(:,20)*data.mass;				% convert to Nm, positive for dorsiflexion moment (flexion)
	data.M_k = d(:,19)*data.mass;				% convert to Nm, positive for plantarflexion moment (extension)
	%=================================NM====================================
    
    %=================================NM====================================
	% joint angles and moments
	% joint force at ankle
	Fx = d(:,8)*data.mass*g;
	Fy = d(:,9)*data.mass*g;
	% transform from thigh coordinate system (Orthotrak) to global coordinates
	% Orthotrak: Fx,Fy>0 when leg pushes pelvis posteriorly and up (see, e.g sit-stand-sit results)
	Fx = -Fx.*cos(data.phi_2) - Fy.*sin(data.phi_2);
	Fy = Fy.*cos(data.phi_2) - Fx.*sin(data.phi_2);
     %=================================NM====================================
	
	if (gait)
		% resample the data at a 50% phase shift, so we get force at other leg
		% BUT ONLY IF GAIT IS SYMMETRIC!
		% this will cause a discontinuity when the data is not perfectly periodic
		F = [Fx Fy];					% put Fx and Fy data together in matrix
		F = [F ; F];					% make two gait cycles
		t = [data.t_perc ; data.t_perc + 100];	% and corresponding original percentage times
		F = interp1(t, F, data.t_perc  + 50);	% resample from 50% to 150% of the original gait cycle
		Fx = F(:,1);					% extract Fx and Fy again
		Fy = F(:,2);
	end
	data.Fx_h = Fx;
	data.Fy_h = Fy;
    
	%=================================NM (not used)========================
	% joint force at knee
	% we leave this in shank reference frame, this is how the sensors see it
	Fx = d(:,8)*data.mass*g;
	Fy = d(:,9)*data.mass*g;
	% transform from thigh coordinate system (Orthotrak) to shank XY coordinates
	% Orthotrak: Fx,Fy>0 when shank pushes thigh posteriorly and up (see, e.g sit-stand-sit results)
	data.Fx_k = -Fx;
	data.Fy_k = Fy;
	%=================================NM====================================
    
    %=================================NM====================================
	% link lengths for thigh (L1) and shank (L2)
	data.L1 = mean(sqrt((data.x_k-data.x_h).^2 + (data.y_k-data.y_h).^2));
	data.L2 = mean(sqrt((data.x_a-data.x_k).^2 + (data.y_a-data.y_k).^2));
    %=================================NM====================================
end
