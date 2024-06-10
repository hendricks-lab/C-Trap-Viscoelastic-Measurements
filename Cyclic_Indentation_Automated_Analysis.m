%% BEN G. - 2024 
% Code for analysis of cyclic indentation & retraction traces.

clear all

directoryPath = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\feb_20\droplet_4";
keyword = "yield";
%dropletRadiusDirectory = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\BEN_DATA_ANALYSIS\jan_30_analysis\DropletRadius\January_30_droplet_radii.xlsx";
%radiusMatrix = readmatrix(dropletRadiusDirectory);
files = dir(directoryPath);
files = files(~[files.isdir]); % removes directories from the file
probeRadius = 0.5; % um


%% Create the data structure that will hold all of parameter values associated with each droplet
% numGroups = length(files);
% numParameters = 6; % this can be changed based on the number of the indentation curves that we opt for
% dropletResultDataStruct = struct();
% % create some other 
trapStiffness = [];
% dropletStifness = []; 
% dropletPreStress_05 = [];

%% NO NEED FOR NOW
% for i = 1:numGroups
%     groupName = sprintf('Group%d', i);
%     dropletResultDataStruct.(groupName) = struct();
% 
%     for j = 1:numParameters
%         varName = sprintf('Var%d', j);
%         dropletResultDataStruct.(groupName).(varName) = NaN;
%     end
% end    
colors = ["b","g","r"];

for k = 1:1 %length(files)
    
    %% Accesing the individual files of the analysis folder
    % Get the file name.
    file = files(k).name;
   
    if contains(file,keyword)
        
        % Construct full file path.
        fullFilePath = fullfile(directoryPath, file);
        %h5disp(fullFilePath)
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
        dropletStiffness = [];
        dropletFitOffset = [];

        if true == true
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

        trapDisplacement = 1.79; % um
        displacementTime = 3.80; % seconds
        numIndices = round(displacementTime * sampling_rate);
        cyclingSize = 50000;
        
    
        possibleStartTimes = [];
        possibleEndTimes = [];
        chuncksize = 5000; % can be fine tune. 
    
        % loop throught the data now
        for m = 1:chuncksize:(length(t) - numIndices)
            increase = xAOD(m + numIndices) - xAOD(m);
            cyclingSize = 500000; % this can be optimized

            if increase >= trapDisplacement
                if length(possibleStartTimes) == 0
                    possibleStartTimes = [possibleStartTimes, t(m)];
                    possibleEndTimes = [possibleEndTimes, t(m + cyclingSize)];
                end
                condition = true;

                if length(possibleStartTimes) > 0 && t(m) > (possibleStartTimes(end) + 2)
                    possibleStartTimes = [possibleStartTimes, t(m)];
                    
                    while condition == true
                        if m + cyclingSize <= length(t)
                            possibleEndTimes = [possibleEndTimes, t(m + cyclingSize)];
                            condition = false;
                        else
                            cyclingSize = cyclingSize - 50000;
                        end
                    end
                end
            end
        end

        for i = 1:length(possibleStartTimes)

            window_size = 2000;
            start = possibleStartTimes(i);
            ends = possibleEndTimes(i);
            
            if round(start*sampling_rate) - 100 <= 0
                zeroing_position = xAOD(1);
                zeroing_force_average = mean(force_1x(1:round(start*sampling_rate) + 5000));
                zeroed_position = xAOD(1:round((ends*sampling_rate) + 5000)) - zeroing_position;
                zeroed_force = force_1x(1:round((ends*sampling_rate) + 5000)) - zeroing_force_average;

            else
                zeroing_position = xAOD(round(start*sampling_rate) - 100);
                zeroing_force_average = mean(force_1x(round(start*sampling_rate)-5000:round(start*sampling_rate)-500));
                zeroed_position = xAOD(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_position;
                zeroed_force = force_1x(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_force_average;
            end

            %zeroing_position = xAOD(round(start*sampling_rate) - 100);
            %zeroing_force_average = mean(force_1x(round(start*sampling_rate)-5000:round(start*sampling_rate)-500));
            %zeroed_position = xAOD(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_position;
            %zeroed_force = force_1x(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_force_average;
            time = ([0:(1/sampling_rate):(length(zeroed_force) -1)*(1/sampling_rate)]);
            smoothed_data = movmean(zeroed_force, window_size);

            % plots the force vs time & trap position vs time for all of
            % the cyclic indentations. 

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

            % Plot the indididual force indentation curves. 
            % condition_1 = true;
            % index = 0;
            % while condition_1 == true
            %     indentationCycleTime = ends - start;
            %     indentationCycleNumPoints = (round(indentationCycleTime*sampling_rate) - (1000*index))/2;
            % 
            %     if indentationCycleNumPoints < length(droplet_indentation)
            %         conditions = false;
            %     else 
            %         index = index + 1;
            %     end
            % end

            indentationCycleTime = ends - start;
            indentationCycleNumPoints = (round(indentationCycleTime*sampling_rate))/1.65;

                            
            % splitting up the traces into indentation and retraction
            % segments

            x_indentation = droplet_indentation(1:indentationCycleNumPoints - 1);
            x_retraction = droplet_indentation(indentationCycleNumPoints:end);
            y_indentation = smoothed_data(1:indentationCycleNumPoints - 1);
            y_retraction = smoothed_data(indentationCycleNumPoints:end);
            


            figure(2);
            hold on;
            if true == true

                % trying to plot the line of best fit onto the figures.
                %p = polyfit(x_indentation,y_indentation, 1);
                %dropletStiffness = [dropletStiffness, p(1)];
                %dropletFitOffset = [dropletFitOffset, p(2)];

                plot(x_indentation, y_indentation, 'blue');
                plot(x_retraction, y_retraction, 'red');
                %plot(zeroed_position,smoothed_data, colors(k));
                %xlabel("Droplet Indentation (nm)");
                xlabel("Trap Displacement (\mum)")
                ylabel("Force (pN)");
            end
        end

        xlabel("Droplet Indentation (nm)");
        ylabel("Force (pN)");
        hold off;

        % Area calculation between the curves - set to true is you want to
        % conduct this analysis
        if true == true
            % calculate the area between the indentation and retraction curves.
            % inverse the retraction array.
    
            x_retraction_inversed = flip(x_retraction);
            y_retraction_inversed = flip(y_indentation);
    
            % Set the domain where the graphs do not cross over
            domain_start = 1200 ;
            domain_end = 1650;
    
            % must interpolate the indentation and retraction domains. 
            common_start = max([min(x_indentation), min(x_retraction_inversed), domain_start]);
            common_end = min([max(x_indentation), max(x_retraction_inversed), domain_end]);
            x_common = linspace(common_start, common_end, max(numel(x_indentation), numel(x_retraction_inversed)));
            
            % Remove duplicate x-values and their corresponding y-values
            [x_indentation_unique, ia, ~] = unique(x_indentation);
            y_indentation_unique = y_indentation(ia);

            [x_retraction_unique, ia, ~] = unique(x_retraction_inversed);
            y_retraction_unique = y_retraction_inversed(ia);

            % Sort the unique values
            [x_indentation_sorted, sortIndex] = sort(x_indentation_unique);
            y_indentation_sorted = y_indentation_unique(sortIndex);

            [x_retraction_sorted, sortIndex] = sort(x_retraction_unique);
            y_retraction_sorted = y_retraction_unique(sortIndex);

            y_indentation_interp = interp1(x_indentation_sorted,y_indentation_sorted,x_common, 'linear', 'extrap');
            y_retraction_interp = interp1(x_retraction_sorted,y_retraction_sorted,x_common, 'linear', 'extrap');
    
            % calculate the difference between the curves
            difference = abs(y_indentation_interp-y_retraction_interp);
    
            % calculate the area between the curves
            area = trapz(x_common,difference);
        end
    end
end

