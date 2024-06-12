% Ben G. 2024
% Membrane bending model

clear all

% Part 1. Membrane bending energy as a function of membrane indentation.
% Interested in the GUV bending energies over a range of indendation depths
% and internal pressures.

simulatedDropletRadius = 1000; % nm
simulatedIndentation = linspace(0,150,50); % nm
simulatedIndentationRadius = sqrt(simulatedDropletRadius.^2 - (simulatedDropletRadius-simulatedIndentation).^2); % nm
simulatedIndentationVolume = (1/3)*(pi)*(simulatedIndentation.^2).*((3.*simulatedDropletRadius)-simulatedIndentation);

[U_MBP1, U_MBP2, U_MBP3, U_MBP4] = membrane_bending_energy(simulatedIndentationVolume);

figure(1)
plot(simulatedIndentation, U_MBP1, '--',LineWidth=3);
hold on
plot(simulatedIndentation, U_MBP2, '--', LineWidth=3);
hold on
plot(simulatedIndentation, U_MBP3, '--', LineWidth=3);
hold on
plot(simulatedIndentation, U_MBP4, '--', LineWidth=3);
set(gca, 'YScale', 'log')
xlabel("Membrane Indentation (nm)");
ylabel("Membrane Bending Energy (J)");

% Part 2. Using hookes law, establish the total mechanical stress felt by
% the membrane during bending

E_modulus = 1000; % Pa
epsion_strain = linspace(0,(pi - 1),25);

[sigma_stress] = membrane_bending_stress(epsion_strain,E_modulus);

figure(2)
plot(epsion_strain,sigma_stress, '--',LineWidth=3);
xlabel("Strain")
ylabel("Stress (Pa)")

dropletRadiusArray = linspace(simulatedDropletRadius/2,simulatedDropletRadius, 10);

figure(3)
red = [1,0,0];
blue = [0,0,1];
hold on
for n = 1:length(dropletRadiusArray)
    % Interpolate the color for the current line
    fraction = (n-1)/(length(dropletRadiusArray)-1);
    current_color = (1-fraction) * blue + fraction * red;

    Y_dc = (sigma_stress.*(dropletRadiusArray(n))*(10^(-9)))/2;
    plot(sigma_stress,Y_dc, 'Color',current_color,'LineStyle','-',LineWidth=3);
    
end    

xlabel("Stress (Pa)");
ylabel("Droplet - Medium Interfacial Tensions (N/m)")

% Part 3. JKR theory for the elastic energy penalties of deforming the
% membrane
E = [100,1000,10000];
R = [500,750,1000];
area_under_curve_below_x = [];

for i = 1:length(E)

    [U_penalties,U_hertz,U_pressure, U_adhesive, U_total] = energy_summation(1000,R(i), 5000, 1000,simulatedIndentationVolume,(10^(-3)),simulatedIndentation);
    
    figure(4)
    %plot(simulatedIndentation,U_penalties, '--',LineWidth=2,Color=red);
    %hold on
    %plot(simulatedIndentation,U_hertz, '--',LineWidth=3);
    %hold on
    %plot(simulatedIndentation, U_pressure, '--',LineWidth=3);
    %hold on
    %plot(simulatedIndentation, -U_adhesive, '--',LineWidth=2,Color=blue);
    %hold on
    plot(simulatedIndentation, U_total, '-',LineWidth=2);
    hold on
    yline(0);
    hold on

    negative_U_total = U_total; % Make a copy of U_total
    negative_U_total(negative_U_total > 0) = 0;  % Set positive values to zero
    area_under_curve_below_x = [area_under_curve_below_x, abs(trapz(simulatedIndentation, negative_U_total))];

end
% Determine the regions of favorable and unfavorable free energy

% for n = 2:length(U_total)
%     if U_total(n) >= 0
%         index = n;
%         break
%     end
% end
% 
% fill([simulatedIndentation(1:index), fliplr(simulatedIndentation(1:index))], [U_total(1:index), zeros(size(U_total(1:index)))],[0.3010 0.7450 0.9330], 'LineStyle', 'none');
% hold on
% fill([simulatedIndentation(index:end), fliplr(simulatedIndentation(index:end))], [U_total(index:end), zeros(size(U_total(index:end)))],[0.9 0.2 0.3], 'LineStyle', 'none');

xlabel("Membrane Indentation (nm)");
ylabel("Free Energy (J)");

figure(5)
plot(R, area_under_curve_below_x, 'o', LineWidth=2,Color=[0 0 0]);
hold on





















