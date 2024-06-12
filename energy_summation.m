function [U_penalties,U_hertz,U_pressure, U_adhesive, U_total] = energy_summation(E_mod, R_d, R_v,P, V, y_dc, d_indent)


    R_eff = 1/((1/R_d) +(1/R_v));
    E_eff = 1/(((((1+0.45)^2)/E_mod)+(((1+0.45)^2)/1000)));
    c = (8/15)*(sqrt(3));
    U_hertz = (c)*(E_eff)*(sqrt(R_eff*(10^(-9)))).*((d_indent*(10^(-9))).^(5/2)); % factor of -9 to putt nm into m
    U_pressure = P.*(V.*(10^-9)^(3));
    U_disipative = (P*(pi)*((R_d*(10^(-9))))*(d_indent.*(10^(-9))).^2)/2;
    U_ycm = pi*(y_dc).*(d_indent.*(10^-9)).^2;
    U_penalties = U_hertz + U_pressure + U_ycm + U_disipative;
    U_adhesive = (R_eff*(10^(-9)))*y_dc.*(d_indent*10^(-9));
    U_total = U_penalties-U_adhesive;

end

