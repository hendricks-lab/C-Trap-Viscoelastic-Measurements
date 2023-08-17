% Ben G - 2023
% The following code is used to fit protein condensate force indentation
% curves. 
%
% INPUTS - A raw force indentation H5 files.
% USER INTERACTION/INPUTS - 
%   - enter the correct calibration file: '/Calibration/4/Force 1x','Offset
%   (pN)' --> the 4 corresponds to the calibration file position.
%   - enter the time interval that corresponds to your force trace
%   - enter the initial force value to zero the data
%   - manually enter the desired fitting ranges for the Hertz and Hookian
%     model
% OUTPUTS - 
%   - FIGURE 1: raw force and trap position data
%   - FIGURE 2: smoothed, zeroed, and selected force trace
%   - FIGURE 3: force indentation curve fitted with a hertz contact model
%   and linear hookian spring model.
%   - FIGURE 4: effective elastic modulus with respect to droplet
%   indentation (DONT USE)
%   - FIGURE 5: log scale elastic modulus (Pa) with respect to droplet
%   indentation. Fitted parameter appear in red as a horizontal line

clear all

file = "linear_displacement_1.5um.h5"; % enter your file here - make sure it is in yuor current matlab directory or a filepath had been added to the code
h5disp(file)
folder_path = "Figures";
force_1x = h5read(file, '/Force HF/Force 1x');

num_points = length(force_1x);
sampling_rate = h5readatt(file,'/Force HF/Force 2x','Sample rate (Hz)');
time_step = 1/sampling_rate;
t = [0:time_step:0+(num_points-1)*time_step];
xAOD = h5read(file, '/Trap position/1X'); 
yAOD = h5read(file, '/Trap position/1Y');
offset = h5readatt(file, '/Calibration/4/Force 1x','Offset (pN)');

% First plot is used to take initial guess on the zeroing force and desired time intervals
figure()
yyaxis left
plot(t, force_1x);
ylabel('Force (pN)');
yyaxis right
plot(t, xAOD);
ylabel('Trap position (um)');
xlabel("Time (s)")
title("Force and Zeroed Trap Position");

time_cut = 1; % set to 1 if you want to make a time cut

if time_cut == 1
    start = 0.2; %(s)
    stop = 1.25; %(s)
    trap_position = xAOD(((start*sampling_rate) +1):(stop*sampling_rate));
    calibrated_force_1x = force_1x(((start*sampling_rate) +1):(stop*sampling_rate));
    time = [0:time_step:0+(length(trap_position)-1)*time_step];
else
    trap_position = xAOD;
    calibrated_force_1x = force_1x;
    time = ([0:time_step:0+(length(force_1x)-1)*time_step]);
end

RF_1x = h5readatt(file, '/Calibration/4/Force 1x','Response (pN/V)'); 
%RD_1x = h5readatt(file,'/Calibration/1/Force 1x','Rd (um/V)');
RF_1y = h5readatt(file, '/Calibration/4/Force 1y','Response (pN/V)'); 
%RD_1y = h5readatt(file,'/Calibration/1/Force 1y','Rd (um/V)');
TK_calibrated = h5readatt(file, '/Calibration/4/Force 1x','kappa (pN/nm)');

% trap zeroing will be slightly different - consider the trap positiona at
% start of experiement as zero

zeroing_position = xAOD(1);
zeroed_position = trap_position - zeroing_position;

% similar zeroing protocol to the sinusoisal waveforms 
zeroing_average = 21.7 ; % hardcoded zeroing guess from raw data from figure 1
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
% (BEN G 2023 - August 3) Incorporate the option for droplet indentation or
% extension

indent = 1;
if indent == 1
    center_disp = smoothed_data/TK_calibrated; % (nm)
    droplet_indentation = (zeroed_position*1000) - center_disp;

% slightly change the geometry and frame of reference for droplet extension
else
    x_0 = zeroed_position(1);
    x_trap = x_0 - zeroed_position;
    center_disp = smoothed_data/TK_calibrated;
    droplet_indentation = x_trap - center_disp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN THE NON-LINEAR REGRESSION - make this into a callable function!!!
% NEXT: find a way to remove the peaks

bead_radius = 0.250;  %um
droplet_radius = 1.000; %um
eff_R = 1/(1/bead_radius + 1/droplet_radius);
x = droplet_indentation./1000; % put indentation in um
y = smoothed_data./1000; % nN
poisson_ratio = 0.3; %lit for soft matter spheres at ambient temperature

%determine the fit range
fit_range_hertz = x <= 0.500; % this gives us the indices that we want to keep
fit_range_hook = x >= 0.500;
x_fit_hertz = x(fit_range_hertz);
y_fit_hertz = y(fit_range_hertz);
x_fit_hook = x(fit_range_hook);
y_fit_hook = y(fit_range_hook);

% define the non-linear model - where params(1) is the effective elastic
% modulus
model_function_hertz = @(params,x) (4/3)*params(1)*eff_R^(1/2)*x.^(3/2); % apply on indentation range comparable to the sie of the probe radius
model_function_hook = @(params,x) params(1)*x + params(2); % apply a regular hookian spring model at larger indentation.
% initial parameter guess
initial_params_hertz = [5000];
initial_params_hook = [1, -100];

opt_params_hertz = lsqcurvefit(model_function_hertz, initial_params_hertz, x_fit_hertz, y_fit_hertz);
opt_params_hook =  lsqcurvefit(model_function_hook, initial_params_hook, x_fit_hook, y_fit_hook);
fitted_values_hertz = model_function_hertz(opt_params_hertz,x_fit_hertz);
fitted_values_hook = model_function_hook(opt_params_hook,x_fit_hook);
% Plot original data and fitted curve
figure();
%plot(droplet_indentation, smoothed_data);
plot(x, y, '-', x_fit_hertz, fitted_values_hertz, '-', x_fit_hook, fitted_values_hook, '-', 'MarkerSize', 1);
legend('Data', 'Hertzian Contact Model', 'Hookian Spring Model');
xlabel('Droplet Indentation (\mum)');
ylabel('Force (nN)');
title('Force vs Droplet Indentation');
legend('show', 'Location','northeast');
%saveas(gcf, fullfile(folder_path, "Fitted Force Indentation Curves_4"));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ben G - 2023 
% save the treated variable to a folder to be able to plot them aferwards

% make this into a user input at somepoint
want_to_save = 0;

if want_to_save == 1
    savedirectory = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\Data_Saves";
    savefilename = file + "_data.mat";
    savefilepath = fullfile(savedirectory, savefilename);
    save(savefilepath, 'x', 'y');
    disp("You have save the processed data.")
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEN G 2023 - Trying to extract the effective modulus as a function of
% droplet indentation and comapring it to the parameter obtained in the
% non-linear lsq curve fittting. 

E_data_hertz = (3/4)*y.*eff_R^(-1/2).*x.^(-3/2); % isolating the effective mod from the hertz contact equation
figure()
plot(x,E_data_hertz, "o", 'MarkerSize',3);
ylim([-20,20])
xlabel("Droplet Indentation (nm)");
ylabel("Effective Modulus - E*"); % nN/um^2 = 10^5 pascals
hold on;
best_fit_hertz = abs(opt_params_hertz);
line([min(x), max(x)], [best_fit_hertz, best_fit_hertz], 'Color', 'red', 'LineStyle', '-');
hold off;

Elastic_mod_hertz = E_data_hertz.*(1-poisson_ratio^2)*(10^5); % multiply by 100000 to put it in Pa
best_fit_hertz_pascal = best_fit_hertz*(1-poisson_ratio^2)*(10^5);

figure()
semilogy(x,Elastic_mod_hertz, "o", 'MarkerSize',2)
xlabel("droplet indentation (\mum)");
ylabel("Elastic Modulus (Pa)");
hold on;
line([min(x), max(x)], [best_fit_hertz_pascal, best_fit_hertz_pascal], 'Color', 'red', 'LineStyle', '-');
hold off;
title("Elastic Modulus with Respect to Droplet Indentation - Hertz Contact Model");
%saveas(gcf, fullfile(folder_path, "Elastic Modulus with Respect to Droplet Indentation - Hertz Contact Model_4"));





