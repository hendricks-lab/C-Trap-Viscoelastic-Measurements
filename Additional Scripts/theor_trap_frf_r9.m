function Mthr=theor_trap_frf_r9(u,fdata,ydata1,phdata1,W_M,W_PH,ydata2,fsamp,Rbead,Zbead,kT,falias)%,f0);

% r7 - added a different weight vector for the phase (Loïc April 2019)
% r8 - added frf_r8 (Loïc May 2019)
i=sqrt(-1);

n1=numel(ydata1); % lenght of the mag data
n2=numel(phdata1); % lenght of the phase data
ffrf=fdata(1:n1); % frequency of the mag data
%fph=fdata([n1+1:n2]); % frequency of phase data (same as frequency of mag data)
fps=fdata((n1+n2+1):end); % frequency of Power Spectrum data
W_mag = W_M.*numel(fps)/(n1+n2); % Weight of the mag, scaled by the the number of points used for PS, and divided by the number of points used for the FRF (mag + phase)
W_phase = W_PH.*numel(fps)/(n1+n2); % Weight of the phase, scaled by the the number of points used for PS, and divided by the number of points used for the FRF (mag + phase)
%Wf=Wf.*exp(-ffrf./100);%Wf-ffrf./Wf;%
%Wf=Wf.*exp(-ffrf./50);

[Mth,PHth]=frf_r9(u,ffrf,Rbead,Zbead,kT,falias);
Mthr1=W_mag.*(Mth-ydata1)./ydata1;

if isempty(phdata1)==1 %zph=0
    PHthr1=[];
else
    PHthr1=W_phase.*(PHth-phdata1)./phdata1; 
end
     
Pth_t=ps_r6(u,fps,fsamp,Rbead,Zbead,kT,falias);%,
IPth=(cumtrapz(fps,Pth_t)-ydata2)./ydata2;
IPth(1)=0;

Mthr=[Mthr1,PHthr1,IPth];
end