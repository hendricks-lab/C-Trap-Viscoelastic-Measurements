%% Ben G. - 2024
% Trying to automate the indentation analysis for higher throughput

%clear all

directoryPath = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\feb_2_analysis";
keyword = "measurement";
dropletRadiusDirectory = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\feb_2_analysis\DropletRadius\February_2_droplet_radii.xlsx";
radiusMatrix = readmatrix(dropletRadiusDirectory);
files = dir(directoryPath);
files = files(~[files.isdir]); % removes directories from the file
probeRadius = 0.5; % um


%% Create the data structure that will hold all of parameter values associated with each droplet
numGroups = length(files);
numParameters = 6; % this can be changed based on the number of the indentation curves that we opt for
dropletResultDataStruct = struct();
% create some other 
trapStiffness = [];
dropletStifness = []; 
dropletPreStress_05 = [];

for i = 1:numGroups
    groupName = sprintf('Group%d', i);
    dropletResultDataStruct.(groupName) = struct();

    for j = 1:numParameters
        varName = sprintf('Var%d', j);
        dropletResultDataStruct.(groupName).(varName) = NaN;
    end
end    

for k = 1:length(files)
    
    %% Accesing the individual files of the analysis folder
    % Get the file name.
    file = files(k).name;
   
    if contains(file,keyword)
        
        % Construct full file path.
        fullFilePath = fullfile(directoryPath, file);
        h5disp(fullFilePath)
        infoStructure = h5info(fullFilePath);
        calibrationName = infoStructure.Groups(2).Groups.Name;
        group_name = sprintf('Group%d', k); % assign the group name for the data structure
    
        force_1x = h5read(fullFilePath, '/Force HF/Force 1x');
        num_points = length(force_1x);
        sampling_rate = h5readatt(fullFilePath,'/Force HF/Force 2x','Sample rate (Hz)');
        time_step = 1/sampling_rate;
        t = [0:time_step:0+(num_points-1)*time_step];
    
        xAOD = h5read(fullFilePath, '/Trap position/1X'); 
        yAOD = h5read(fullFilePath, '/Trap position/1Y');
        offset = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','Offset (pN)');
       
        RF_1x = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','Response (pN/V)');
        RF_1y = h5readatt(fullFilePath, string(calibrationName) + '/Force 1y','Response (pN/V)');
        TK_calibrated = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','kappa (pN/nm)');
        trapStiffness = [trapStiffness, TK_calibrated];

        if true == false
            figure(1)
            yyaxis left
            plot(t, movmean(force_1x,5000));
            ylabel('Force (pN)');
            yyaxis right
            plot(t, xAOD);
            ylabel('Trap position (um)');
            xlabel("Time (s)");
            title("Smooth Zeroed Force and Trap Position");
        end
    
        %% trying to automate the recognition of indentation intervals and the segmentation of time vectors 
        
        trapDisplacement = 0.35; % um
        displacementTime = 0.5; % seconds
        numIndices = round(displacementTime * sampling_rate);
    
        possibleStartTimes = [];
        possibleEndTimes = [];
        chuncksize = 5000; % can be fine tune. 
    
        % loop throught the data now
        for i = 1:chuncksize:(length(t) - numIndices)
            increase = xAOD(i + numIndices) - xAOD(i);
            
            if increase >= trapDisplacement
                possibleStartTimes = [possibleStartTimes, t(i)];
                possibleEndTimes = [possibleEndTimes, t(i + numIndices)];
            end
        end
    
        % select the actual start and end times
        trueStartTimes = [];
        trueEndTimes = [];
        for i = 2:3:length(possibleEndTimes)
            trueStartTimes = [trueStartTimes, possibleStartTimes(i)];
            trueEndTimes = [trueEndTimes, possibleEndTimes(i)];
        end 
    
        %% automate the data fitting
        
        dropletRadius = radiusMatrix(k); % based on the way the automated file renaiming code works - the first file that is read is the last measurement taken.
        
        
        for i = 1:length(trueEndTimes)
           
            var_name = sprintf('Var%d', i);

            window_size = 5000;
            start = trueStartTimes(i);
            ends = trueEndTimes(i);
            
            zeroing_position = xAOD(round(start*sampling_rate) - 100);
            zeroing_force_average = mean(force_1x(round(start*sampling_rate)-5000:round(start*sampling_rate)-500));
            zeroed_position = xAOD(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_position;
            zeroed_force = force_1x(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_force_average;
            time = ([0:(1/sampling_rate):(length(zeroed_force) -1)*(1/sampling_rate)]);
            smoothed_data = movmean(zeroed_force, window_size);
            
            % if you want to print the individual indentation curves
            if true == false
                figure()
                yyaxis left
                plot(time, smoothed_data);
                ylabel('Force (pN)');
                yyaxis right
                plot(time, zeroed_position);
                ylabel('Trap position (um)');
                xlabel("Time (s)");
                title("Smooth Zeroed Force and Trap Position");
                 
            end
            
            indent = 1;
            if indent == 1
                center_disp = smoothed_data/TK_calibrated; % (nm)
                droplet_indentation = (zeroed_position*1000) - center_disp;
            end
            
            %% Begin the hertzian fitting for the individual force indentation curves
            eff_R = 1/(1/probeRadius + 1/dropletRadius);
            x = droplet_indentation./1000; % um
            y = smoothed_data./1000; % nN
            poisson_ratio = 0.3;

            fit_range_hertz = x <= 0.4;
            x_fit_hertz = x(fit_range_hertz);
            y_fit_hertz = y(fit_range_hertz);

            model_function_hertz = @(params,x) (4/3)*params(1)*eff_R^(1/2)*x.^(3/2); % apply on indentation range comparable to the sie of the probe radius
            initial_params_hertz = [500];
            opt_params_hertz = lsqcurvefit(model_function_hertz, initial_params_hertz, x_fit_hertz, y_fit_hertz);
            fitted_values_hertz = model_function_hertz(opt_params_hertz,x_fit_hertz);

            %% Extracting the young's modulus
            E_data_hertz = (3/4)*y.*eff_R^(-1/2).*x.^(-3/2);
            best_fit_hertz = abs(opt_params_hertz);
            Elastic_mod_hertz = E_data_hertz.*(1-poisson_ratio^2)*(10^5); % multiply by 100000 to put it in Pa
            best_fit_hertz_pascal = best_fit_hertz*(1-poisson_ratio^2)*(10^5);

           
            figure(2);
            hold on;

            if true == true
                
                %plot(x,y,"r");
                p = polyfit(x,y,1);
                dropletStifness = [dropletStifness, p(1)];
                % xlabel("Droplet Indetation (\mum)");
                % ylabel("Force (nN)");
            end
            % xlabel("Droplet Indetation (\mum)");
            % ylabel("Force (nN)");
            % hold off;

            % saving this to the data structure
            dropletResultDataStruct.(group_name).(var_name) = best_fit_hertz_pascal;

        end
        xlabel("Droplet Indentation (\mum)");
        ylabel("Force (nN)");
        hold off;
        dropletResultDataStruct.(group_name).(sprintf('Var%d',(numParameters-1))) = dropletRadius; % assign the droplet radius as a variable in the droplet data strcuture. 
       
        %% determine the prestress on a droplet

        if true == true
            averageEndForce = mean(force_1x((length(force_1x)-5000):(length(force_1x))));
            dropletPreStress_05 = [dropletPreStress_05, (zeroing_force_average - averageEndForce)];
            %preStress = zeroing_force_average - averageEndForce;
            %dropletResultDataStruct.(group_name).(sprintf('Var%d',numParameters)) = preStress;
        end
    end

    if true == false
        saveFileName = "DropletResult.mat";
        saveFilePath = fullfile(directoryPath, saveFileName);
        save(saveFilePath,'dropletResultDataStruct');
    end
end

% find the range of trap stifnesses used for the experiments
trapStifnessMin = min(trapStiffness);
trapStiffnessMax = max(trapStiffness);

%% Plotting the histogram best fit stifnesses
if true == true
    figure; % Opens a new figure window
    initialBinEdges = linspace(min(dropletStifness), max(dropletStifness), 10);
    histogram(dropletStifness, initialBinEdges,'Orientation', 'horizontal','FaceColor', 'r');
    %title('Histogram of Data');
    xlabel('Counts');
    ylabel('Droplet Stiffness (pN/nm)');
    % Get the current axes
    ay = gca;
    
    % Change the format of the x-axis tick labels to show actual values instead of exponential notation
    ay.YAxis.Exponent = 0;
    ay.YAxis.TickLabelFormat = '%.4f';
end

%% plotting the histograms
if true == false
    figure;
    initialBinEdges_prestress = linspace(min(dropletPreStress_05), max(dropletPreStress_05), 7);
    histogram(dropletPreStress_05, initialBinEdges_prestress);
    hold on
    histogram(dropletPreStress,initialBinEdges_prestress);
    xlabel("Pre-Stress (pN)");
    ylabel("Counts");
    hold off
end










