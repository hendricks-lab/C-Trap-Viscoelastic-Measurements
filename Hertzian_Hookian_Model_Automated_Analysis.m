%% Ben G. - 2024
% Automated indentation analysis code.

% HOW IT WORKS: The data is segmented into N time point chuncks. It then
% iterates throught the entire time series using the N time point chunks.
% At each chunk, the code will check if there has been an increase in xAOD
% position. If there is an increase that matches the expected xAOD
% displacement - that time chunck is identified as being a potential
% indentation curve. 
%
% PARAMETERS TO CONTROL:
% trapDisplacement => i.e. if trap displacement is 400nm set to 350nm
% displacementTime => set to approximate time over which the displacement
% occurs
% chunkSize => if the trap displacement velocity is high - you maybe
% require a smaller chunk size - basically behaves as the sampling rate of
% the time series.
%
% Saving results data structure
% go to the end of the script and set the saving parameters 

clear all

directoryPath = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\mar_12_analysis";
keyword = "measurement";
dropletRadiusDirectory = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\mar_12_analysis\DropletRadius\March_12_droplet_radii.xlsx";
radiusMatrix = readmatrix(dropletRadiusDirectory);
files = dir(directoryPath);
files = files(~[files.isdir]); % removes directories from the file
probeDiameter = 0.5; % um

%% Create the data structure that will hold all of parameter values associated with each droplet
numGroups = length(files);
numParameters = 6; % this can be changed based on the number of the indentation curves that we opt for
dropletResultDataStruct = struct();
% create some other 
allDropletYoungModulusFit = [];
allDropletRadius = [];
trapStiffness = [];
dropletStifness = []; 
dropletPreStress = [];
maxDropletStiffness = [];
fittedslope1 = [];
fittedslope2 = [];

for i = 1:numGroups
    groupName = sprintf('Group%d', i);
    dropletResultDataStruct.(groupName) = struct();

    for j = 1:numParameters
        varName = sprintf('Var%d', j);
        dropletResultDataStruct.(groupName).(varName) = NaN;
    end
end    

for k = 5:5 %length(files)
    
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

        if TK_calibrated > 0.08
            disp("test");
            continue
        end
           
        if true == true 
            figure(1)
            % yyaxis left
            plot(t, movmean(force_1x,5000));
            ylabel('Force (pN)');
            % yyaxis right
            % plot(t, xAOD);
            % ylabel('Trap position (um)');
            xlabel("Time (s)");
            title("Smooth Zeroed Force and Trap Position");
        end
    
        %% trying to automate the recognition of indentation intervals and the segmentation of time vectors 
        
        trapDisplacement = 0.38; % um
        displacementTime = 0.48; % seconds
        numIndices = round(displacementTime * sampling_rate);
   
        possibleStartTimes = [];
        possibleEndTimes = [];
        chuncksize = 3500; % can be fine tune. 
    
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
        DropletYoungModulusFit = [];
        allSlopes = [];
        allCenterPoints = [];

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
            
            %% Ben G. 2024 - Try to extract the time derivative of Force vs time curve
            % if you want to print the individual indentation curves
            figure(2)
            hold on
            if true == false
                
                % Because of the extremely high sampling rate and
                % stochastic movement whitin the trap volume - a non
                % continuous differentiation approach will be used

                interval = 1000;
                diffForce = [];
                diffTime = [];
       
                for n = 1:(length(time)/interval)
                    diffForce = [diffForce, (smoothed_data(interval*n) - smoothed_data(interval*(n-1) +1))];
                    diffTime = [diffTime, ((n-1)*(interval/sampling_rate)) + (interval/2)/sampling_rate];
                    
                end

                forceDerivative = diffForce/(interval/sampling_rate);

                subplot(3,1,1);
                plot(time,smoothed_data, 'LineWidth', 2);
                hold on
                xlabel("Time (s)");
                ylabel("Force (pN)");
                title("Force vs. Time");

                subplot(3,1,2);
                plot(diffTime, forceDerivative, 'o','MarkerSize', 4);
                hold on
                xlabel("Time (s)");
                ylabel("Force Rate (pN/s)"); %  WHAT TO CALL THIS IDK??? 
                title("Force Rate vs. Time");

                subplot(3,1,3);
                plot(time, zeroed_position,'red', 'LineWidth', 2);
                hold on
                xlabel("Time (s)");
                ylabel("Trap Displacement (nm)");
                title("Trap Displacement vs. Time");
            end
            
            indent = 1;
            if indent == 1
                center_disp = smoothed_data/TK_calibrated; % (nm)
                droplet_indentation = (zeroed_position*1000) - center_disp;
            end
            
            [E_data_hertz, dropletYoungModulusFit, x, y, eff_R] = fitHertzModel(droplet_indentation, smoothed_data, probeDiameter, dropletRadius, 0.45);

            DropletYoungModulusFit = [DropletYoungModulusFit, dropletYoungModulusFit];
            

            %% Extracting the Apparent Droplet Stiffness with Respect to Indentation (pN/nm vs nm)
            
            if true == true
        
                [slopes, centerPoints] = calculateSlopes(droplet_indentation,smoothed_data,20, 2000);
                allSlopes = [allSlopes, slopes];
                allCenterPoints = [allCenterPoints, centerPoints];
                maxDropletStiffness = [maxDropletStiffness,max(slopes)];
                indentationDomain = linspace(0,250,25);
                
            end
      
            % plot of all indentation curves in the data set.
            figure(4);
            hold on;
            if true == false
                x_plot = x(1:1000:end);
                y_plot = y(1:1000:end);
                % hLine = plot(x,y, 'Color', [0, 0, 1, 0.25], LineWidth=1);
                % set(gcf, 'Renderer', 'OpenGL');
                % set(hLine, 'Color', [1, 0, 0, 0.5]); % Set the color with an alpha value
                % Scatter plot with semi-transparent markers
                % set(gcf, 'Renderer', 'OpenGL');
                % scatter(x_subset, y_subset, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1);

                % Use the patch function to plot lines with transparency
                for k = 1:length(x_plot)-1
                    patch('XData', x_plot(k:k+1), 'YData', y_plot(k:k+1), 'EdgeColor', [0.573, 0.643, 0.969], ...
                          'LineWidth', 2, 'EdgeAlpha', 0.1, 'FaceColor', 'none');
                    hold on; % Keep the plot active to draw the next segment
                end

                xlabel("Droplet Indentation (\mum)");
                ylabel("Force (nN)");
                hold on
            end
            % saving this to the data structure
            dropletResultDataStruct.(group_name).(var_name) = dropletYoungModulusFit;

        end
        averageModulus = mean(DropletYoungModulusFit)/(10^5); % divided by 10^5 to be in appropriate units to plot
        
         %% Simulated force indentation plot using average modulus values for a droplet. 
       
        if true == true
            simulatedIndentation = linspace(0,max(x), 20);
            simulatedForce = (4/3)*averageModulus*sqrt(eff_R)*(simulatedIndentation.^(3/2));
          
            %% Simulate force change with respect to indentation
            simulatedForceDerivative = 2*averageModulus*(eff_R)*sqrt(simulatedIndentation);
        end
            
        if true == false
            figure(4)
            hold on
            plot(simulatedIndentation,simulatedForce, 'ok', 'MarkerFaceColor','k','MarkerSize', 4, 'DisplayName','Simulated Force Identation Curve');
            dummyPlot = plot(NaN,NaN, 'ok', 'MarkerFaceColor','k','MarkerSize', 4);
            lgd = legend('show');
            legend(dummyPlot, 'Simulated Force Indentation Curve');
            lgd.Box = 'off';
            legend('Location', 'northwest');
        end

        xlabel("Droplet Indentation (\mum)");
        ylabel("Force (nN)");
        hold off;
        dropletResultDataStruct.(group_name).(sprintf('Var%d',(numParameters-1))) = dropletRadius; % assign the droplet radius as a variable in the droplet data strcuture. 
       
        %% determine the prestress on a droplet

        if true == true
            averageEndForce = mean(force_1x((length(force_1x)-5000):(length(force_1x))));
            dropletPreStress = [dropletPreStress, (zeroing_force_average - averageEndForce)];
            preStress = zeroing_force_average - averageEndForce;
            dropletResultDataStruct.(group_name).(sprintf('Var%d',numParameters)) = preStress;
        end
    end

    if true == false
        saveFileName = "DropletResult.mat";
        saveFilePath = fullfile(directoryPath, saveFileName);
        save(saveFilePath,'dropletResultDataStruct');
    end

    [slope1,slope2,intercept1,intercept2,xdomain] = logfitting(allSlopes,allCenterPoints);
    
    fittedslope1 = [fittedslope1, slope1];
    fittedslope2 = [fittedslope2, slope2];

    log_y_fit_1 = slope1 * log10(xdomain) + intercept1;
    % Convert the y values back from log scale to the original scale
    y_fit_1 = 10.^(log_y_fit_1);
    log_y_fit_2 = slope2 * log10(xdomain) + intercept2;
    % Convert the y values back from log scale to the original scale
    y_fit_2 = 10.^(log_y_fit_2);
    
    figure(3)
    % Overlay the linear fit line
    loglog(allCenterPoints,allSlopes,'o',Color=[0 0 0], LineWidth=2);
    hold on
    loglog((simulatedIndentation*1000), simulatedForceDerivative, "--",'LineWidth', 2, 'Color', 'r');
    hold on
    loglog(xdomain, y_fit_1, 'b-', 'LineWidth', 2);
    hold on
    loglog(xdomain, y_fit_2, 'g-', 'LineWidth', 2);
    hold on 
    % xlabel('Droplet Indentation (nm)');
    % ylabel('Effective Droplet Stifness (pN/nm)');
    hold off
end

if true == false

    simulatedIndentation = linspace(0,max(x), 20);
    R = 0.2;
    mod = 583/100000; % conversion factor for that the units work out. 
    simulatedForce = (4/3)*mod*sqrt(R)*(simulatedIndentation.^(3/2));
    figure(4)
    hold on
    plot(simulatedIndentation,simulatedForce, 'ok', 'MarkerFaceColor','blue','MarkerSize', 4, 'DisplayName','Simulated Force Identation Curve');
    dummyPlot = plot(NaN,NaN, 'ok', 'MarkerFaceColor','k','MarkerSize', 4);

end

%% Plotting the histogram best fit stifnesses - to be used for linear fits. 
if true == false 
    figure; % Opens a new figure window
    initialBinEdges = linspace(min(dropletStifness), max(dropletStifness), 10);
    histogram(dropletStifness, initialBinEdges,'Orientation', 'horizontal', 'FaceColor', 'r');
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
    initialBinEdges_prestress = linspace(min(dropletPreStress), max(dropletPreStress), 10);
    histogram(dropletPreStress,initialBinEdges_prestress);
    xlabel("Pre-Stress (pN)");
    ylabel("Counts");
end


%% saving the the Droplet results data structure to a folder in current directory

if true == false
    savingDir = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\DropletResults";
    experimentInfo = "april_18_newprotein_01FCRNA";

    newFileName = sprintf('dropletResultDataStruct_%s.mat',experimentInfo);
    savingFilePath = fullfile('DropletResults', newFileName);
    save(savingFilePath,"dropletResultDataStruct");
end

if true == false
    savingDir = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\DropletResults";
    experimentInfo = "mar_5_slope1_newprotein_01FRNA_MAXSTIFNESS";

    newFileName = sprintf('loflogfit_%s.mat',experimentInfo);
    savingFilePath = fullfile('DropletResults', newFileName);
    save(savingFilePath,"fittedslope1");
end


% droplet stifness distribution bootstrapping
if true == false
    [bootstrapMean,ci,bootstrapMeans] = boostrap(maxDropletStiffness,1000);
    % Display results
    fprintf('Bootstrap Mean: %f\n', bootstrapMean);
    fprintf('95%% Confidence Interval: [%f, %f]\n', ci(1), ci(2));
    
    % Fit a normal distribution to the bootstrap means
    pd = fitdist(bootstrapMeans, 'Normal');
    
    % Evaluate the probability density function (PDF) for the fitted distribution
    x_values = linspace(min(bootstrapMeans), max(bootstrapMeans), 100);
    pdf_values = pdf(pd, x_values);
    
    % Plot the histogram of the bootstrap means
    figure(5);
    histogram(bootstrapMeans,'Normalization','pdf',FaceColor='yellow');
    hold on;
    
    % Overlay the PDF of the fitted distribution
    plot(x_values, pdf_values, 'r-', 'LineWidth', 2);
    % xlabel('Mean Value');
    % ylabel('Probability Density');
    % title('Bootstrap Sample Means with Normal Fit');
    % legend('Bootstrap Means', 'Normal Fit');
    % hold off;

    xlabel("Effective Maximum Droplet Stiffeness (pN/nm)");
    ylabel('Probability Density');
    title('Bootstrap Sample Means with Normal Fit');
    %legend('Bootstrap Means', 'Normal Fit');
    hold off;
end








