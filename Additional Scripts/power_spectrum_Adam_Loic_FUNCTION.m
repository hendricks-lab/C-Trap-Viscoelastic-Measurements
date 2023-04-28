% Umich                          14.08.2007              adam g. hendricks
% Edited + extra comments        07.03.2019              Loïc Chaubet
%
% This code calculates the power spectrum (PS) of a signal by using the 
% Pwelch PSD estimate and by integrating over a finite frequency range. 
% It also substracts any general trend observed in the signal (before PSD). 
%
% NOTE: using a 10x smaller Nfft yields 10x more averaging, which is very 
% similar to what this code does if Nint = 10. One would need to test the 
% two to see what's best. However, this power spectrum is used strickly for 
% the thermal response (thicker baseline band), which may not affect the
% downstream analysis in any significant way. It would be better to verify
% yourself and with Adam.
% 
% PS = [frequency; power spectrum]
% data = [timevector,signal]
% Nwindow = number of points per window
% Noverlap = number of points to overlap (usually 0.5*Nwindow)
% Nfft = number of points to calculate FFT (frequency resolution)
% Nint = number of frequencies to integrate for power spectrum (Nint>=2)
% Fs = sampling frequency
% Fps = frequency of PS
% Mps = magnitude of PS (i.e. density, in Power/Hz)

function [Fps,Mps]=power_spectrum_Adam_Loic_FUNCTION(data,Nwindow,Noverlap,Nfft,Nint,Fs)

%remove trends by fitting a 2nd order polynomial to the whole signal, and
%substracting the fit from the data. This does not seem to have any
%significant impact unless the data has a strong trend.
pfit=polyfit(data(1,:),data(2,:),2); 
data(2,:)= data(2,:) - polyval(pfit,data(1,:));

% Using pwelch method to obtain the PSD estimate
[Pxx,F]=pwelch(data(2,:),Nwindow,Noverlap,Nfft,Fs);

%integrate the PSD to find power spectrum
Pxx_sum=cumtrapz(F,Pxx); % integrated
ji=1;
jf=ji+Nint;
k=1;
while jf<=length(Pxx)
    Fps(k)=0.5*(F(ji)+F(jf)); % frequency bin
    Mps(k)=(Pxx_sum(jf)-Pxx_sum(ji))/(F(jf)-F(ji)); % density of power over the frequency bin
    ji=ji+Nint;
    jf=ji+Nint;
    k=k+1;
end
Fps = Fps';
Mps = Mps';
%PS=[Fps',Mps']; % 