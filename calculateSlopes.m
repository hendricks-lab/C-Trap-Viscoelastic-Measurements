function [slopes,centerPoints] = calculateSlopes(indentation,force, numPoints, rangeN)
    
    if nargin < 4
        error('Not enough input arguments. Expected indentation, force, numPoints, and rangeN.');
    end
    
    % Determine the step size over witch the slop will be calculated
    stepSize = floor((length(force) - rangeN) / (numPoints - 1));
    
    % Initialize the slope and center vectors
    slopes = zeros(1, numPoints);
    centerPoints = zeros(1, numPoints);

    for i = 1:numPoints
        % Calculate the start and end indices of the range of points
        startIdx = (i - 1) * stepSize + 1;
        endIdx = startIdx + rangeN - 1;
        
        % Fit a line to the data within the current range
        p = polyfit(indentation(startIdx:endIdx), force(startIdx:endIdx), 1);
        
        % Store the slope and the center point of the range
        slopes(i) = p(1);
        centerPoints(i) = indentation(startIdx + floor(rangeN / 2));
    end
    
    % Now you have the slopes at the 10 points and the center points of the ranges
    % If you want to plot these slopes as a separate figure:
    % figure (3);
    % plot(centerPoints, slopes, 'o', 'LineWidth', 2);
    % hold on
    % xlabel('Indentation (nm)');
    % ylabel('Effective Droplet Stiffness (pN/nm)');
    % title('Effective Droplet Stiffness with Respect to Indentation');
end

