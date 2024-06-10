function [E_data_hertz,dropletYoungModulusFit, x, y, eff_R] = fitHertzModel(indentation,force,probeDiameter,dropletRadius,poissonRatio)

    eff_R = 1 / (1/(probeDiameter/2) + 1/(dropletRadius/2));

     % Convert indentation and force data to the appropriate units (if necessary)
    x = indentation / 1000; % um to nm
    y = force / 1000; % pN to nN

    % Define the fitting range based on indentation values
    fit_range_hertz = x <= 0.30;
    x_fit_hertz = x(fit_range_hertz);
    y_fit_hertz = y(fit_range_hertz);

    % Define the Hertz model
    model_function_hertz = @(params, x) (4/3)*params(1)*eff_R^(1/2)*x.^(3/2);

    % Initial parameters for the Hertz model fitting
    initial_params_hertz = [500]; % Example initial modulus value

    % Perform the least squares curve fitting
    [opt_params_hertz, ~] = lsqcurvefit(model_function_hertz, initial_params_hertz, x_fit_hertz, y_fit_hertz);

    % Calculate the absolute value of the fitted modulus parameter
    best_fit_hertz = abs(opt_params_hertz);
    best_fit_hertz_pascal = best_fit_hertz*(1-poissonRatio^2)*(10^5);

    % Calculate the Young's modulus from the Hertz model
    E_data_hertz = (3/4)*best_fit_hertz.*eff_R^(-1/2).*x_fit_hertz.^(-3/2);

    % Convert to Pascal
    Elastic_mod_hertz = E_data_hertz.*(1-poissonRatio^2)*(10^5); % Conversion factor included
   
    % Collect the fitted modulus values
    dropletYoungModulusFit = best_fit_hertz_pascal;
end

