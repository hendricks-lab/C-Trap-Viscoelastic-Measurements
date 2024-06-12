function [U_MBP1, U_MBP2, U_MBP3, U_MBP4] = membrane_bending_energy(V)
    %km = 12.5*(1.38*10^(-23)); % J
    P_0 = [(100),(10^3),(10^4),(10^5)]; % Turgor Pressure of the vesicle (Pa); Harmandaris & Deserno 2006
    %U_MBE = (11) * (sqrt((pi^3)/(32)))*(h_0./R_0)*km;
    %U_MBE(1) = 0;
    U_MBP1 = ((P_0(1))*(V.*(10^-9)^3));
    U_MBP2 = ((P_0(2))*(V.*(10^-9)^3));
    U_MBP3 = ((P_0(3))*(V.*(10^-9)^3));
    U_MBP4 = ((P_0(4))*(V.*(10^-9)^3));
    
    %U_MBT = U_MBE + U_MBP;
end

