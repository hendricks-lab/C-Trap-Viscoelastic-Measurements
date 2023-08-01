% BEN G 2023
% analysis code for linear and step wise trap displacements
% Multiple figure output to help the user situate themselves in the
% analysis. 

clear all

file = "contact_1.h5"; 
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

% trap zeroing will be slightly different - consider the trap positiona at
% start of experiement as zero

zeroing_position = xAOD(1);
zeroed_position = xAOD - zeroing_position;

% similar zeroing protocol to the sinusoisal waveforms 
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

% calculate displacement from center
center_disp = smoothed_data/TK_calibrated; % (nm)
droplet_indentation = (zeroed_position*1000) - center_disp;

figure()
yyaxis left
plot(time, droplet_indentation);
ylabel('Droplet Indentation (nm)');
yyaxis right
plot(time, zeroed_position);
ylabel('Trap position (um)');
xlabel("Time (s)")
title("Time Trace of Droplet Indentationa and Trap Position")
%saveas(gcf, fullfile(folder_path, "Smoothed Data (5%50%200nm1Hz_2)"))

figure()
yyaxis left
plot(time, droplet_indentation);
ylabel('Droplet Indentation (nm)');
yyaxis right
plot(time, smoothed_data);
ylabel('Force (pN)');
xlabel("Time (s)")
title("Time Trace of Droplet Indentationa and Force Signal")
%saveas(gcf, fullfile(folder_path, "Smoothed Data (5%50%200nm1Hz_2)"))

% try to automate finding the indices where a trap step occurs

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTINUE POSSIBLE STEPS WISE ANALYSIS. VERY SIMPLE TO DO OUTSIDE OF
% MATLAB BUT COULD BE COOL TO AUTOMATE. 
% set the threshold for detecting significant increases - i only want it to
% register jumps that are greater than the backgroun signal. Take the min
% and max values of initial background signal prior to force extention

%indent_min = min(droplet_indentation(1:10000)); 
%indent_max = max(droplet_indentation(1:10000));
%indent_threshold = indent_max - indent_min;

%force_min = min(smoothed_data(1:10000)); 
%force_max = max(smoothed_data(1:10000));
%force_threshold = force_max - force_min;

% Analysis of step wise jumps to determing the apparent elasticity of the
% droplets. Detect a steps wise jump in the force and indentation signal.

%step_indent = diff(droplet_indentation);
%step_force = diff(smoothed_data);

%indent_step_indices = find(abs(step_indent) > 0.5*indent_threshold);














