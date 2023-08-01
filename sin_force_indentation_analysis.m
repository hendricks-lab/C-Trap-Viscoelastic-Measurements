% Trying to reformate the code to deal with corrupted H5 files
% Going to be the exact same workflow as with the
% "force_extension_analysis_condensates.m" file
% 
% Curently being using to process data from July 27th experiments
%
% Make sure to change the save file names
%% Ben G. 2023
% Force extension analysis for single bead sinusoidal oscillations on a condensate.
% INPUTS: 
%   1. An experimental sinusoidal oscillation force file ("file")

%   2. A calibration file corresponding to the desired experimental
%   conditions ("desired_calibration")

% OUTPUTS:
%   1. An overlaid time trace of raw force signal and trap position during
%   droplet compression (figure 1).

%   2. A smoother and zeroed overlaid time trace of force signal and trap position during
%   droplet compression (figure 2).

%   3. A force indentation curve displaying droplet indentation (nm)
%   against the measured force signal (pN) - the positive force portion of
%   the plot is fitted using a linear regression model (figure 3).
%
% (Ben) July 28th -
%       Create a version that does not require loading  seperate calibration
%       files (calibration files are easily corrupter so loading them
%       individually may be risky).
%       Added an in_phase and out_phase to deal with difference in
%       experimental set up. You make the decision based on Figures 1 & 2.
%       set the desire model = 1 and proceed.
%

clear all

file = "5%OVP_50%TS_buffer_control_200nmAMP_1Hz.h5"; 
folder_path = "Figures";
force_1x = h5read(file, '/Force HF/Force 1x');

num_points = length(force_1x);
sampling_rate = h5readatt(file,'/Force HF/Force 2x','Sample rate (Hz)');
time_step = 1/sampling_rate;
t = [0:time_step:0+(num_points-1)*time_step];
xAOD = h5read(file, '/Trap position/1X'); 
yAOD = h5read(file, '/Trap position/1Y');
offset = h5readatt(file, '/Calibration/5%OVP_50%TS_2/Force 1x','Offset (pN)');

% to take initial guess on the zeroing force and desired time intervals
figure()
yyaxis left
plot(t, force_1x);
ylabel('Force (pN)');
yyaxis right
plot(t, xAOD);
ylabel('Trap position (um)');
xlabel("Time (s)")
title("Force and Zeroed Trap Position");

%%%%% DO YOU WANT TO SELECT A PORTION OF THE FORCE TRACE? 
%%%%% ENTER THE DESIRED TIME WINDOW BELOW

time_cut = 0; % set to 1 if you want to make a time cut

if time_cut == 1
    start = 0 ; %(s)
    stop = 9; %(s)
    trap_position = xAOD(((start*sampling_rate) +1):(stop*sampling_rate));
    calibrated_force_1x = force_1x(((start*sampling_rate) +1):(stop*sampling_rate));
    time = [0:time_step:0+(length(trap_position)-1)*time_step];
else
    trap_position = xAOD;
    calibrated_force_1x = force_1x;
    time = [0:time_step:0+(length(force_1x)-1)*time_step];
end

RF_1x = h5readatt(file, '/Calibration/5%OVP_50%TS_2/Force 1x','Response (pN/V)'); 
%RD_1x = h5readatt(file,'/Calibration/1/Force 1x','Rd (um/V)');
RF_1y = h5readatt(file, '/Calibration/5%OVP_50%TS_2/Force 1y','Response (pN/V)'); 
%RD_1y = h5readatt(file,'/Calibration/1/Force 1y','Rd (um/V)');
TK_calibrated = h5readatt(file, '/Calibration/5%OVP_50%TS_2/Force 1x','kappa (pN/nm)');

% because we use the correct calibration before hand - we do not need to
% use the recalibration portion of the code - we can simply use the raw
% force data

xAOD_min = min(trap_position);
xAOD_max = max(trap_position);
zeroing_position = (xAOD_max+xAOD_min)/2; 
zeroed_position = trap_position - zeroing_position;

zeroing_average = 8; % hardcoded zeroing guess from raw data from figure 1
zeroed_calibrated_force_1x = calibrated_force_1x - zeroing_average;

% Smooth the raw force data using a rolling window filter
window_size = 1000;
smoothed_data = movmean(zeroed_calibrated_force_1x, window_size);

figure()
yyaxis left
plot(time, smoothed_data);
ylabel('Force (pN)');
yyaxis right
plot(time, zeroed_position);
ylabel('Trap position (um)');
xlabel("Time (s)")
title("Smooth Zeroed Force and Trap Position")
%saveas(gcf, fullfile(folder_path, "Smoothed Data (5%50%200nm1Hz_2)"))

in_phase = 0;
if in_phase == 1
    center_disp = smoothed_data/TK_calibrated; % (nm)
    droplet_indentation = (zeroed_position*1000) - center_disp;
    % Filter data points to include only positive y-values
    positiveIndices = smoothed_data> 0;
    xPositive = droplet_indentation(positiveIndices);
    yPositive = smoothed_data(positiveIndices);
    
    % Perform linear regression using filtered data points
    coefficients = polyfit(xPositive, yPositive, 1);
    slope = coefficients(1);
    intercept = coefficients(2);
    
    % Generate values for the fitted line
    xFit = min(droplet_indentation):0.1:max(droplet_indentation);
    yFit = slope * xFit + intercept;
   
    % Figure 2. Plot the original data points and the fitted line
    figure()
    plot(droplet_indentation,smoothed_data); % Original data points in blue
    hold on;
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Fitted line in red
    hold off;
    % Add labels and title to the plot
    xlabel('Droplet Indentation (nm)');
    ylabel('Force (pN)');
    title('Droplet Force Indentation Curve');
    legend('Data Points', 'Fitted Line', 'Location', 'northwest');
    parametersText = sprintf('Slope: %.2f\nIntercept: %.2f', slope, intercept);
    annotation('textbox', [0.15, 0.7, 0.25, 0.1], 'String', parametersText, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
    %saveas(gcf, fullfile(folder_path, "Force Extension (5%50%200nm1Hz_2)"))
end 

out_phase = 1;
if out_phase == 1
    positive_trap_indices = zeroed_position > 0; % gets the indices where the trap positioin is greater than zero (the indentation) 
    positive_position = zeroed_position(positive_trap_indices); % only select the indentation portion of the displacement
    positive_force = smoothed_data(positive_trap_indices); % only select the force associated to the indentation portion of the displacement. 

    center_disp = positive_force/TK_calibrated;
    droplet_indentation = (positive_position*1000) - center_disp;
    
     % Perform linear regression using filtered data points
    coefficients = polyfit(droplet_indentation, positive_force, 1);
    slope = coefficients(1);
    intercept = coefficients(2);
    
    % Generate values for the fitted line
    xFit = min(droplet_indentation):0.1:max(droplet_indentation);
    yFit = slope * xFit + intercept;
  
    % Figure 2. Plot the original data points and the fitted line
    figure()
    plot(droplet_indentation,positive_force); % Original data points in blue
    hold on;
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Fitted line in red
    hold off;
    % Add labels and title to the plot
    xlabel('Droplet Indentation (nm)');
    ylabel('Force (pN)');
    title('Droplet Force Indentation Curve');
    legend('Data Points', 'Fitted Line', 'Location', 'northwest');
    parametersText = sprintf('Slope: %.2f\nIntercept: %.2f', slope, intercept);
    annotation('textbox', [0.15, 0.7, 0.25, 0.1], 'String', parametersText, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
    %saveas(gcf, fullfile(folder_path, "Force Extension (5%50%200nm1Hz_2)"))
end