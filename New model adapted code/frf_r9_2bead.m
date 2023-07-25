function [mag,phase]=frf_r9_2bead(u,f,Rbead,Zbead,kT,falias)%,f0);

% r8:
% added abs(Phase) as the output  (Loïc May 2019) but why?
% why??!

% r9 (Loïc) August 2019
% added gamma_water at 37 degrees (3.2577e-06) [pN*s/nm] 

%r9_2bead (Ben) July 2023
% modifying the TF to fit the new model.


%u=[Minf, fp, Dfit, phi]
j=sqrt(-1);
Minf=1/u(1);
gammar=u(2); % is the drag coefficient as a multiple of that of water.
alpha=u(3);
ktrap=u(4);
kcyt0=u(5);
kcyt2 = u(6); % (Ben 2023) aded new term. 
kcyt1=u(7);
%keq0=u(6);
%w0=keq0/(gammar*9.42e-6);
%kcyt1=(u(6)-ktrap-kcyt0)/(w0^alpha);
% wm0=(gammar*9.42e-6)/m;
% kcyt1=(keq0-ktrap-kcyt0)/(wm0^alpha);

m=u(8)*1e-21; % (BEN 2023) reindexed this portion
nu=u(9)*1e12;

gamma_frf=gammar*3.2577e-6.*ones(size(f));%gamma_r5(u,f,Rbead,Zbead,kT);

wfrf=2*pi.*f;

% kcyt=kcyt0+kcyt1.*wfrf.^alpha;
% keq=ktrap+kcyt;
% num_real=-m.*wfrf.^2-imag(gamma_frf).*wfrf + kcyt;
% num_imag=real(gamma_frf).*wfrf;
% den_real=keq-m.*wfrf.^2-imag(gamma_frf).*wfrf;
% den_imag=real(gamma_frf).*wfrf;
% Mth=Minf.*sqrt(num_real.^2+num_imag.^2)./...
%     sqrt(den_real.^2+den_imag.^2);
% PHth=atan2(num_imag,num_real)-...
%     atan2(den_imag,den_real);
% Mpd=(2*pi.*falias)./sqrt(wfrf.^2+(2*pi.*falias).^2);

% (Ben - July 2023) first itteration of the new transfer function. 
% Assuming no mass & no difference in close and long range mechanics.
% (BEN 2023) added kcyt2 term to the denominator of the tranfer function.
Num = kcyt0 + gamma_frf.*wfrf.*j;
Den = (gamma_frf - 0.5.*gamma_frf).*wfrf.*j + (kcyt0 + kcyt2 + ktrap);

%Num=kcyt0-m.*wfrf.^2+gamma_frf.*wfrf.*j+(kcyt1.*(j.*wfrf).^(alpha))./gamma(alpha);
%Den=ktrap+kcyt0-m.*wfrf.^2+gamma_frf.*wfrf.*j+(kcyt1.*(j.*wfrf).^(alpha))./gamma(alpha);

Mth=Minf.*sqrt(real(Num).^2+imag(Num).^2)./...
     sqrt(real(Den).^2+imag(Den).^2);
 PHth=atan2(imag(Num),real(Num))-...
     atan2(imag(Den),real(Den));
Mpd=(2*pi.*falias)./sqrt(wfrf.^2+(2*pi.*falias).^2);

PHpd=-atan2(wfrf,2*pi*falias);
mag=Mth.*Mpd;
phase=(180/pi).*(PHth+PHpd);
phase = abs(phase); % why?

end