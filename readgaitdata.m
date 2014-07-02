
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
%   P_a             ankle power
%   M_a             ankle moment
%   phi_a           ankle angle

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

    % Do you want the program to use the gait data from an Excel file or from the Motek system.
    Excel   = 0;        % use the Excel data (1) or not (0gait.p)
    Motek   = 1;        % use the Motek data (1) or not (0)
    

% If the Excel data is chosen to be used
if (Excel)
	% first read the subject data sheet
	d = xlsread(filename, 'Subject', '', 'basic');
	data.subjectID = d(3,1);
	data.mass = d(4,1);
    
    

	% now read the movement data
	d = xlsread(filename, movement, '', 'basic');
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
    
    keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Added by Rick Rarick 8/18/2010
    
    data.F_rot_ankle = d(:, 9)*data.mass*g;  % convert to N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% hip, knee, ankle XY position data, convert from mm to m
	data.x_h = d(:,25)/1000; data.y_h = d(:,27)/1000;
	data.x_k = d(:,28)/1000; data.y_k = d(:,30)/1000;
	data.x_a = d(:,31)/1000; data.y_a = d(:,33)/1000;
	
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
    
    
	% joint angles and moments and power output for the ankle
	data.phi_a = d(:,4)*pi/180;				    % convert to ankle angle in rad, positive for dorsiflexion
	data.M_a = -d(:,19)*data.mass;				% convert to Nm, positive for dorsiflexion
    data.P_a = d(:,22)*data.mass;               % power data for ankle, converted to Watts
    
    
    
    % Joint angular velocities from manual calculation
    nsamples = size(d,1);
    omega = NaN(nsamples,1);
    for i=2:nsamples-1                         % Loop starts at 2 becuase angle analysis does not start at 0
        omega(i) = (data.phi_a(i+1) - data.phi_a(i-1))...
               /(data.t(i+1) - data.t(i-1));
    end
    data.P_a = omega .* data.M_a;
    plot(data.t_perc,data.P_a)

end



% If the motek data is chosen to be used
if (Motek)
    % Use the Motek data instead of the Excel data
    motek = load('T017_S12');
    
    % Rename the motek structure variable to 'data'. The purupse is so we have less code to modify
    % The apostrophie is to convert the motek data's row vectors into column vectors
    data = motek;
    data.phi_a      = motek.Right_Ankle_PlantarFlexion_Angle_Mean';    % Ankle Angle
    data.M_a        = motek.Right_Ankle_PlantarFlexion_Moment_Mean';   % Ankle Moment
    data.omega_a    = motek.Right_Ankle_PlantarFlexion_Rate_Mean';     % Ankle angular velocity
    data.t_perc     = motek.percent_gait_cycle';
    data.t          = (motek.mean_cycle_duration*(0:99)/100)';
    keyboard
    % Error checking the time data
    if (data.t_perc(1) ~= 0)
		error('Time(%) does not start at zero');
	elseif (data.t_perc(end) == 100)
		error('Time(%) ends at 100');
	elseif (std(diff(data.t_perc)) > 0.001)
		error('Time(%) is not equally spaced.');
	end
	data.cycle = max(data.t) + data.t(end) - data.t(end-1);
    
    % Multiplying angular velocity and moment to get power. 
    % The function 'times' is a element-by-element multiplication function for matrices
    data.P_a        = times(data.M_a, data.omega_a);                  %Ankle power
    plot(data.t_perc,data.P_a)

    %keyboard
end

end