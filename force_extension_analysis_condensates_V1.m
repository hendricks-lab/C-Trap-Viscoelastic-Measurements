% Ben G. 2023
% Force extension analysis for single bead oscillations on a condensate.
% Produced a time trace of zeroed force signal and trap position.
% as well as a force extention curve corresponding to droplet indentation.

clear all

file = "droplet_0.5um_1Hz_15%.h5"; 
%h5disp(file);
% The number of indivual data points in the force_1x file is determined

% Extracting force, position, time and sampling data from the h5 files obtained via the
% ctrap bluelake software
force_1x = h5read(file, '/Force HF/Force 1x');
%force_1y = h5read(file, '/Force HF/Force 1y');
%force_2x = h5read(file, '/Force HF/Force 2x');
%force_2y = h5read(file, '/Force HF/Force 2y');

% determintion of time array
num_points = length(force_1x);
sampling_rate = h5readatt(file,'/Force HF/Force 2x','Sample rate (Hz)');
time_step = 1/sampling_rate;
t = [0:time_step:0+(num_points-1)*time_step];
T = transpose(t);

% Only Trap 1 position can be exported. 
xAOD = h5read(file, '/Trap position/1X'); 
yAOD = h5read(file, '/Trap position/1Y');
offset = h5readatt(file, '/Calibration/2_experiment calibration/Force 1x','Offset (pN)');

% The following parameters are obtained from the bluelake calibration - these will be used to uncalibrate the force measurements. 
% RF is the force response (pN/V). RD is the distance response (um/V) 

%RF_2x = h5readatt(file, '/Calibration/2_experiment calibration/Force 2x','Response (pN/V)'); 
RF_1x = h5readatt(file, '/Calibration/2_experiment calibration/Force 1x','Response (pN/V)'); 
%RD_1x = h5readatt(file,'/Calibration/1/Force 1x','Rd (um/V)');
%RF_2y = h5readatt(file, '/Calibration/2_experiment calibration/Force 2y','Response (pN/V)'); 
RF_1y = h5readatt(file, '/Calibration/2_experiment calibration/Force 1y','Response (pN/V)'); 
%RD_1y = h5readatt(file,'/Calibration/1/Force 1y','Rd (um/V)');

% The following code is used to decalibrate the force data.
%volts_x = (force_1x)/RF_1x; %xQUAD (V)
%volts_y = (force_1y)/RF_1y; %yQUAD (V)
bead_1_volts_x = (force_1x - offset)/RF_1x; %xQUAD (V)
%bead_1_volts_y = (force_1y - offset)/RF_1y; %yQUAD (V)
%bead_2_volts_x = (force_2x - offset)/RF_2x; %xQUAD (V)
%bead_2_volts_y = (force_2y - offset)/RF_2y; %yQUAD (V)

desired_calibration = "Calibration 4_15%_near surface.h5"; % load the desired calibration profile that will be used to decalibrate the data. 
h5disp(desired_calibration);
% extract the desired calibration parameters
RD_calibrated = h5readatt(desired_calibration, '/Calibration/4_15%_near surface/Force 1x' ,'Rd (um/V)')*1000; % *1000 puts the distance reponse in nm.
RF_calibrated = h5readatt(desired_calibration, '/Calibration/4_15%_near surface/Force 1x' ,'Rf (pN/V)'); %force response
TK_calibrated = h5readatt(desired_calibration, '/Calibration/4_15%_near surface/Force 1x' ,'kappa (pN/nm)'); % trap stifness. 

%recalibrate the force data using the new calibrated parameters
calibrated_force1x = bead_1_volts_x * RF_calibrated;

%zero the force and AOD data using the average values prior to the start of trap oscillation 
zeroing_average = mean(calibrated_force1x(1:1000));
zeroed_calibrated_force1x = calibrated_force1x - zeroing_average;
% make trap zeroing is making sense
xAOD_min = min(xAOD);
xAOD_max = max(xAOD);
zeroing_position = (xAOD_max+xAOD_min)/2; % finds the average positions about the oscillation
zeroed_position = xAOD - zeroing_position;

figure()
yyaxis left
plot(t, zeroed_calibrated_force1x);
ylabel('Force (pN)');
yyaxis right
plot(t, zeroed_position);
ylabel('Trap position (um)');
xlabel("Time (s)")
title("Zeroed Force and Trap Position");

%Using Hookes law and the trap's force signal - determine the bead's
%displacement from the trap during the oscillation. 

center_disp = zeroed_calibrated_force1x/TK_calibrated; % (nm)
droplet_indentation = (zeroed_position*1000) - center_disp;

% Determine the transition point betweenn linear regimes.
% Ben G 2023 - this first iteration is done empirically. In future
% iterations, include a weighing function that determines the breakpoint
% automatically for each data set. 
breakpoint = 200; % nm
coefficients = polyfit(droplet_indentation(breakpoint:end), zeroed_calibrated_force1x(breakpoint:end), 1);
slope = coefficients(1);
intercept = coefficients(2);

% Generate fitted line
xFit = min(droplet_indentation(breakpoint:end)):0.1:max(droplet_indentation(breakpoint:end));
yFit = slope * xFit + intercept;

% Plot the data and the fitted line
figure;
plot(droplet_indentation, zeroed_calibrated_force1x); % Original data points
hold on;
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Fitted line
xlabel('Droplet Indentation (nm)');
ylabel('Force (pN)');
legend('Data', 'Fitted Line');
title('Linear Model Fitting');
% Display the parameter values
parametersText = sprintf('Slope: %.2f\nIntercept: %.2f', slope, intercept);
annotation('textbox', [0.15, 0.7, 0.25, 0.1], 'String', parametersText, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

