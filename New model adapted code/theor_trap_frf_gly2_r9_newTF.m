%To fit theoretical frf, ps of a mostly viscous solution
%Fits beta, gamma, ktrap, kcyt1 and alpha
%Kcyt0 is set to zero
%
%23May2021 - Ora 
%            Changed to frf_r9
%      
%24May2021 - Ora
%            uu=[u(1:2),0,u(3),0,0,u(4:5)]; %old
%            uu=[u(1:4),0,u(6:8)]; %new

%24July2023 - Ben
%        modyfying the input to include the second cytoplasmic elasticity
%        term - kcyt2
%
%
%
%
%

function Mthr=theor_trap_frf_gly2_r9(u,fdata,ydata1,phdata1,W_M,W_PH,ydata2,fsamp,Rbead,Zbead,kT,falias)%,f0);
                                 
i=sqrt(-1);

% u =  [Beta, GammaR, Alpha, k_trap, k_cyt0, k_cyt2 k_cyt1, mass, nu];%
% (BEN) reindexed to have 9 terms in the uu array
uu=[u(1:4),0,u(6:9)]; %fz=0;

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