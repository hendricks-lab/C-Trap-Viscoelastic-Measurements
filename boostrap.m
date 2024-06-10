function [bootstrapMean,ci,bootstrapMeans] = boostrap(valueVector,numBootstrap)
    
     % initialize the zeros array
    bootstrapMeans = zeros(numBootstrap,1);

    % remove the NaN values.
    valuesSansNaN = valueVector(~isnan(valueVector));

    % boostrap from this new vector.
    numObservations = numel(valuesSansNaN);

    % Perform bootstrapping
    for i = 1:numBootstrap
        % Generate a bootstrap sample: resample the clean data with replacement
        bootstrapSample = valuesSansNaN(randsample(numObservations, numObservations, true));
        
        % Compute the mean of the bootstrap sample
        bootstrapMeans(i) = mean(bootstrapSample);
    end

     % Calculate the mean and 95% confidence interval of the bootstrap means
    bootstrapMean = mean(bootstrapMeans);
    ci = prctile(bootstrapMeans, [2.5, 97.5]);
end

