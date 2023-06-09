%To fit theoretical frf, ps of a viscous solution
%Fits only beta, gamma and ktrap
%Kcyt0, kcyt1 and alpha are set to zero
%
% r7 - added a different weight vector for the phase (Loïc April 2019)
% r8 - added frf_r8 (Loïc May 2019)
%
%23May2021 - Ora 
%            Changed to frf_r9
%      
%24May2021 - Ora
%            uu=[u(1:2),0,u(3),0,0,u(4:5)]; %old
%            uu=[u(1:2),0,u(4),0,0,u(7:8)]; %new
%
%
function Mthr=theor_trap_frf_gly_r9(u,fdata,ydata1,phdata1,W_M,W_PH,ydata2,fsamp,Rbead,Zbead,kT,falias)%,f0);
                                 
i=sqrt(-1);
% u =  [Beta, GammaR, Alpha, k_trap, k_cyt0, k_cyt1, mass, nu];
uu=[u(1:2),0,u(4),0,0,u(7:8)]; %fz=0;

n1=numel(ydata1);
n2=numel(phdata1);
ffrf=fdata(1:n1);
%fph=fdata([n1+1:n2]);
fps=fdata((n1+n2+1):end);
W_mag = W_M.*numel(fps)/(n1+n2); % Weight of the mag, scaled by the the number of points used for PS, and divided by the number of points used for the FRF (mag + phase)
W_phase = W_PH.*numel(fps)/(n1+n2); % Weight of the phase, scaled by the the number of points used for PS, and divided by the number of points used for the FRF (mag + phase)
%Wf=Wf.*exp(-ffrf./100);%Wf-ffrf./Wf;%
%Wf=Wf.*exp(-ffrf./50);

[Mth,PHth]=frf_r9(uu,ffrf,Rbead,Zbead,kT,falias);
Mthr1=W_mag.*(Mth-ydata1)./ydata1;

if isempty(phdata1)==1 %zph=0
    PHthr1=[];
else
    PHthr1=W_phase.*(PHth-phdata1)./phdata1; 
end

Pth_t=ps_r6(uu,fps,fsamp,Rbead,Zbead,kT,falias);
IPth=(cumtrapz(fps,Pth_t)-ydata2)./ydata2;
IPth(1)=0;

Mthr=[Mthr1,PHthr1,IPth];

end