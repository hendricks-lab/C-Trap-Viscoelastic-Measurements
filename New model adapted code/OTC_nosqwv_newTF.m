% For calibrating optical trap in situ
% Fit magnitude, phase response and power spectrum to get Beta, GammaR, Alpha, k_trap, k_cyt0, k_cyt1

% 17.06.16 - agh load existing transfer function data, fit using simple VE model
% version to fit Loic's multiharmonic data
%
% Jan. 2019 - Modified/commented/cleaned - Loïc
% r7: April 2019 (Loïc)
%       Added phase fitting with low weight and removed coherence-based weighing
%       Added strict bounds on alpha across the two fits
%       Will add interpolation for extra points??!
% r8: May 2019 (Loïc)
%       Changed the automatic initial conditions set-up
%       Using r8 verion of fit funciton and frf trap function
%       Added MSE for more user guidance in terms of fitness of fits (only
%       stored in a few of the WT (BB) data sets..I manually compared the MSE values when fitness decision was difficult.
%
% r9: August 2019 (Loïc)
%       Further added comments
%       Set the gamma of water (drag coefficient) to the right one (not that it changes the fitting anyway, but if we want a quick comparison between the drag coeff. of water and the cytoplasm, we need the right one)
%       Added semi-automated initial guess finding based on previous final fitted values and empirical & analytical equations
%       Further cleaned up the code (Sept. 2019)
%
% r10: Sept. 2019 (Loïc)
%       Close to the final version of the code..!
%       Organized the code to make it easier to read (I hope)
%       Further added comments
%       Made the interactive parts more user-friendly
%
% r10.1: Oct. 2020 (Ora)
%       Super similar to fitting_adam_method_r10, minor bug fixes/user-friendly changes
%       Added an option to use fixnfit (fit parameters iteratively in two subsets)
%       Fixnfit works best in cases where the hardcoded guess is used
%       (often when manual selection gives errors)
%
%r10.2: Feb. 2021 (Ora)
%       Changed min frequency on power spectrum fitting from hardcoded
%       value to selected value
%       Added fixnfit grouping option and fit order option
%       Added possibility to use multistart (doesn't work very well)
%       Added lots of fun information to the file name
%
%r10.3: Mar. 2021 (Ora)
%       Added option to fit yPS instead of xPS
%
%r10.4: 14 Oct 2021 (Ora, but technically Adam)
%      Changed power spectrum scaling    
%
%r11: July 12 2023 (Ben)
%     changing the frf_r9 function to frf_r9_2bead which contains the new
%     transfer function.
%
%r12: July 24 2023 (Ben)
%     adding new parameters that fit the 2 bead free body diagram and the
%     new transfer function. Start by introducing new medium stifness and
%     keep the same drag coefficient. Changing the called functions to
%     those this added parameters. Calling the new transfer function

clc, clear all, %close all %, warning off
addpath("C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\loic_codes-master\loic_codes-master\adams-method") %path where some functions are stored

% file path, file name, initial guess condition, fitting weights (TF weight, relative to power spectrum), min freq. fit for power spectrum, fitting phase or not and a few constants which are usually leftuntouched.

%set directory to multiharmonic measurement
datdir = "C:\Users\Ben\Dropbox\My PC (LAPTOP-L13ADGT7)\Documents\MATLAB\DATA\13July2023\"; % defining the path to where the transferfunction.mat files are
addpath(datdir)

%set filename to multiharmonic measurement
flist{1} = {"21.07.23-13.54-TF-voltage data table_5amp_3xsep_morefrequencies_newcal_minus.mat"}; %name of file to load
filename =  fullfile(datdir,flist{1}{1}); % selecting the file, using the file path and file name

%fitting parameters to set to your liking
% Ben 2023 - semi_auto_ig should be set to zero for non-vivo measurements.
semi_auto_ig = 0; % ig = Initial Guess. note: this only works for IN CELL calibration, as this was based on the data I collected in the cells. In any other viscoelastic solution (or just viscous), this won't work.
load_file_ig = 0; % load the initial guess conditions from a previous iteration that you will have saved
user_ig = 1; % this uses a fully hardcoded initial guess from the user and, similar to the load_file_ig, completely skips the interactive steps for finding a proper initial guess 
multi = 0; %if set to 1, a multistart algorithm will be used to sample initial guesses (requires supporting functions)

usefixnfit = 1; % 0 - fit all parameters together, anything else - use fix'n'fit (-Ora)
fitOrder = 1; %0 varies ktrap 2nd, 1 varies ktrap first

groupem = 1; %ktrap, kcyt0 & kcyt1 vs alpha beta gamma
%groupem = 2; %beta, ktrap, vs kcyt0, kcyt1, alpha, gamma

alphaLB = 0.75; %lower bound on the alpha value (-Ora)
alphaUB = 1; %upper bound on the alpha value (-Ora)

% set initial guess
beta_guess = 20000; % manual input (nm/V) - Ben 2023 the Ctrap provides all of it's distance parameter in um, make sure to convert back to nanometers. 
gammar_guess = 10; %1; %1.68; % manual input, this means the drag coeff. of the cytoplasm is "ig.GammaR" times the one of water
alpha_guess = 0.75; %0.75; % manual input % the higher the alpha - WILL BE SET TO 0 IN COMPLETELY VISCOUS REGIMES
ktrap_guess = 0.98; %0.05; %0.1; % manual input % Ben 2023 - find an estimate for this guess using an active calibration
kcyt0_guess = 0.001; % manual input
kcyt2_guess = 0.01; % BEN 2023 - long range elastic component
kcyt1_guess = 0.0000; %0.002; % manual input % Ben 2023 - set this to zero as there is not currently a frequency dependence term built into the TF

cohthreshh = 0.75; %set coherence threshold for points used in the fit
minDatPoints = 0.5; %min percentage (if you multiplied by 100%) of selected data points from magnitude response required for fitting (-Ora)
 
Wmag = 0.075; %5e-2; % 0.5e-3; weight applied to the magnitude. This is key in balancing the quality of the fit between the FRF and the PS. A low weight here will fit the power spectrum better, but the TF worse. In general, it's better to prioritize the power spectrum as it is more reliable because of the number of points.
Wphase = 0.05.*Wmag;%5e-3; % weight applied to the phase (fitting the phase with low weight usually helps, and sometime up to 0.2*Wmag can be used. But sometime, it's best to leave this low, say 0.01*Wmag or less)
Wlast_points = 0; % Add more weight to the last point (or points) in the TF? This may help fitting the higher frequency points.
zph = 1; % zph = 1 --> fits the phase, zph = 0 --> the phase data is ignored

%some important constants (make sure to check temp and Rbead!)
i = sqrt(-1); % imaginary "i"
Temp = 27; % [deg C]
kT = 1.381e-2*(Temp+273.15); % [K]
Rbead = 500; % Radius of bead [nm]
Zbead = 2000; % Height of bead [nm] (not currently being used (August 2019) but it's something that could be incorporated.

numPS = input('How many power spectra are in this file (1 - just x, 2 - x & y)? ');

if numPS == 1
    load(filename,'Fexc','H','FP','PYm','Pyy','Pyy2','fsamp','C');
    FPx = FP;
    PYm_x = PYm;
    disp('Loading yPS will cause calibration to be done using yPS.')
    getPSy = input('Would you like to load yPS from a separate file? (1 - yes, 2 - no)? ');
    if getPSy == 1
        yPSpath = input('What is the path for the yPS file? (include single quotes) ');
        addpath(yPSpath);
        yPSfile = input('What is the filename for yPS? (include single quotes and file type (.mat)) ');
        load(yPSfile);
        PSfit = '-yPS';
        FPy = FP;
        PYm_y = PYm;
    else
        PSfit = '-xPS';
    end
elseif numPS == 2
    load(filename,'Fexc','H','FP_y','PYm_y','FP_x','PYm_x','Pyy','Pyy2','fsamp','C');
    usePS = input('Which PS do you want to use for calibration? (1 - xPS, 2 - yPS) ');
        if usePS == 1
            FP = FP_x;
            PYm = PYm_x;
            PSfit = '-xPS';
        elseif usePS == 2
            FP = FP_y;
            PYm = PYm_y;
            PSfit = '-yPS';
        end
end

if semi_auto_ig + load_file_ig + user_ig > 1 % this is in case more than one ig condition is toggled to "1"
    disp('Choose only one initial guess (ig) condition..!' )
    return
elseif semi_auto_ig + load_file_ig + user_ig == 0 % this is in the case none of the ig conditions are selected
    disp('Choose one initial guess (ig) condition..!' )
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hall = H; %record original FRF
Fexcall = Fexc; %record original excitation frequencies 
Call = C; % record original magnitude squared coherence
kph = 1; % constant multiplier factor for the photodiode (leave it at 1 unless otherwise specified)
H = H.*exp(i*kph*pi).*exp(i*0.000035*2*pi.*Fexc); %correct for phase lag between acquisition channels 
% NOTE: phase shift measured through applying high frequency oscilattions to a bead stuck on coverslip and finding the correcting factor that would yield a zero phase at those frequencies
% The idea is that because the bead is stuck to the coverslip, any applied force (through the offset of the lazer with respect to the bead) will
% yield a direct measurment, without any viscous contribution of the surrounding medium. i.e., it would be a purely elastic response, i.e. phase should be zero.
%PYm = 2*PYm; % one sided vs. two sided-power spectrum???

%Adam's new thingy
if exist('PYm','var')==1
    PYm = 2.*PYm*mean(diff(FP)); % one sided vs. two sided-power spectrum???
else
  PYm = 2*mean(diff(FP_y)).*PYm_y; disp('try scaling PYm by frequency bin size')
end

% HANDLING THE TRANSFER FUNCTION (TF) POINTS:
%   - Add points to help the fit on higher frequency points? Which points to keep for the fitting? Note: Remember that changing the weight of individual points using W_vector can achieved all that.
%   - Plots the TF points, and ask the user to select the range of TF points to fit
%   - Coherence thresholding
%   - Zoom in on the range of points selected
addpoint = 0; %input('Add copies of a point at the end? 0 = no, 1 = yes ');
if addpoint == 1
point_to_copy = input('Which point is to be copied? (index)');
Npoints = input('Copied how many times?');
end 
if length(H) == length(Fexc) && length(H) == length(C) % checking that freqs, magnitude and coherence are all the same length (-Ora)
    jkeep = 1:length(H); % keep all points. Can change this if you wish to discard the last few points
else
    fprintf('Houston, we''ve got a problem - length(H) = %d, length(Fexc) = %d, and length(C) = %d\n', length(H), length(Fexc), length(C))
end

H = H(jkeep);
Fexc = Fexc(jkeep);
C = C(jkeep);
Mfrf = 20*log10(abs(H));
PHfrf = unwrap(angle(H)); 
PHfrf = abs(PHfrf)*180/pi; % why take the absolute value only?
if addpoint == 1 % if we want to add artificial points at the end of the magnitude plot to help the fit to flatten out after the pole
    act_length = length(Fexc); % keeping track of actual real length of Fexc
    for jj = 1:Npoints
        H(act_length+jj) = H(point_to_copy); % copying H
        C(act_length+jj) = C(point_to_copy); % copying C
        PHfrf(act_length+jj) = PHfrf(point_to_copy); % copying PH
        Pyy(act_length+jj)=Pyy(point_to_copy); % copying Pyy
        Pyy2(act_length+jj)=Pyy2(point_to_copy); % copying Pyy2
        Mfrf(act_length+jj) = 0.999*Mfrf(point_to_copy); % shifting the newly added magnitudes up, just a little bit
        Fexc(act_length+jj) = Fexc(point_to_copy)*(1.5*jj+0.5); % shifting the newly added frequencies to the right
    end
end

figure, subplot(311), semilogx(Fexc,Mfrf,'o','linewidth',2)
ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3],'fontsize',14,'fontweight','bold');
subplot(312), semilogx(Fexc,PHfrf,'o','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3],'fontsize',14,'fontweight','bold');
hfrf=gcf;
subplot(313), semilogx(Fexc,C,'o','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz'), set(gca,'fontsize',14,'fontweight','bold')
disp('Click on minimum acceptable frequency')
txt = 'Click on MIN freq';
subplot(311)
t1 = text(0.012,max(Mfrf)+0.5*(min(Mfrf)-max(Mfrf))/2,0,txt,'FontSize',12);
[freq_low, chert]=ginput(1);
delete(t1)
disp('Click on maximum acceptable frequency')
txt = 'Click on MAX freq';
t1 = text(0.012,max(Mfrf)+0.5*(min(Mfrf)-max(Mfrf))/2,0,txt,'FontSize',12);
[freq_high,chert]=ginput(1);
delete(t1)
% MALINA WAS HERE (somewhere around Nov. 2018?)
% Keeping the data points selected by the user, and that is above the coherence threshold. One could add a coherence threshold check or something. However, low coherence points are usually 
% at lower frequencies which are not used in this calibration. So, right now, this threshold is set to very low and is basically not being used.  

freq_ind = find(Fexc>freq_low&Fexc<freq_high);
%Fexc=Fexc(freq_ind);
%H=H(freq_ind);
%PHfrf = PHfrf(freq_ind);
%C=C(freq_ind);
find_coh=find(C>cohthreshh&C<0.99999999&Fexc>freq_low&Fexc<freq_high);

% if coherence threshold is set too high (or data has very low coherence),
% too many points might get thrown out to get a good fit, so this loop lets
% you reset the num points required or the coherence threshold without
% having to start over. On the other hand, if it gets annoying, you can 
% comment out the whole 'while' loop, but the code will crash if find_coh 
% is empty or has a length of 1. (-Ora)

while (length(find_coh) < floor(length(freq_ind).*minDatPoints)) 
    fprintf('Not enough data points (%d%%, require %0.2f%%) above coherence threshold (%0.4f). Max coherence = %0.4f, median coherence = %0.4f\n', length(find_coh)/length(freq_ind), minDatPoints, cohthreshh, max(C(freq_ind)), median(C(freq_ind))); %find second-to-largest coherence max(C(C<max(C))) %max(C(find(Fexc>freq_low&Fexc<freq_high))))
    what_to_change = input('Would you like to change num points required (input 1) or the coherence threshold (input 2)? ');
    if what_to_change == 1
        minDatPoints = input('What minimum percentage of data points would you like to use? ');
    else
        cohthreshh = input('What coherence threshold would you like to use? ');
    end
    find_coh=find(C>cohthreshh&C<0.99999999&Fexc>freq_low&Fexc<freq_high);
end

%filter for data within appropriate range. All the data is already stored 
%in Fexcall, Hall, Call. (-Ora)
Fexc=Fexc(find_coh);
H=H(find_coh);
PHfrf = PHfrf(find_coh);
C=C(find_coh);
Mfrf=Mfrf(find_coh);
Pyy=Pyy(find_coh);
Pyy2=Pyy2(find_coh);
% Setting x and y plotting ranges so that we get a zoomed in view on the
% points selected (to better see the shape of the curve). I think this can
% be done with "tight" function on MATLAB...but oh well (Sept 9th 2019)
x_upper = freq_high;
x_lower = freq_low;
y_upper = max(Mfrf);
y_lower = min(Mfrf);

%GLOBAL FOR MULTI
if multi == 1
    global absHglob PHfrfglob Wglob Wphaseglob ydat10glob Rbeadglob Zbeadglob faliasglob fsampglob kTglob
    absHglob = abs(H(:))';
    PHfrfglob = PHfrf;
    Rbeadglob = Rbead;
    Zbeadglob = Zbead;
    fsampglob = fsamp;
    kTglob = kT;
end

if addpoint == 1 % i.e. if points were added, they are marked and mentionned in the title
    figure(hfrf), hold on, subplot(311), semilogx(Fexc(end-Npoints+1:end),Mfrf(end-Npoints+1:end),'o','markerfacecolor','g','linewidth',2)
    ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
    title(['FRF with ' num2str(Npoints) ' artificially added points at the end '])
else
    figure(hfrf), hold on, subplot(311), semilogx(Fexc,Mfrf,'o','markerfacecolor','b','linewidth',2)
    ylabel('Mag, dB'), set(gca,'Xlim',[0.01 5e3]);
end
set(gca,'Xlim',[x_lower x_upper]), set(gca,'Ylim',[y_lower y_upper]); % this is where the zoom in happens
subplot(312), hold on, semilogx(Fexc,PHfrf,'o','markerfacecolor','b','linewidth',2)
ylabel('Phase, deg'), xlabel('f, Hz'), set(gca,'Xlim',[0.01 5e3]);
subplot(313), hold on, semilogx(Fexc,C,'o','markerfacecolor','b','linewidth',2)
ylabel('Coherence'), xlabel('f, Hz')

if semi_auto_ig == 1 % in that case, the user has to click values multiple times

    % Manual selection of the zero, pole and DC gain on the TF
    subplot(311)
    disp('Click on zero corner frequeny (Hz)') % asking for the zero
    txt = 'Click on zero corner frequeny (Hz)';
    t1 = text(x_lower,y_upper+0.5*(y_lower-y_upper)/2,0,txt,'FontSize',12);
    [ig.fz_frf,ig.Pz]=ginput(1); % [frequency of the zero, magnitude of the zero (dB)]
    delete(t1)
    disp('Click on pole corner frequency (Hz)') % asking for the pole
    txt = 'Click on pole corner frequency (Hz)';
    t1 = text(x_lower,y_upper+0.5*(y_lower-y_upper)/2,0,txt,'FontSize',12);
    [ig.fp_frf,ig.MinfdB]=ginput(1); % [frequency of the pole, magnitude of the pole (dB)]
    delete(t1)
    disp('Click on DC gain') % asking for DC gain
    txt = 'Click on MIN freq. range for DC gain';
    t1 = text(x_lower,y_upper+0.5*(y_lower-y_upper)/2,0,txt,'FontSize',12);
    [DC_low_freq,DC_low_mag] = ginput(1); % [frequency of the left-most selected point, it's magnitude]
    delete(t1)
    txt = 'Click on MAX freq. range for DC gain';
    t1 = text(x_lower,y_upper+0.5*(y_lower-y_upper)/2,0,txt,'FontSize',12);
    [DC_high_freq,DC_high_mag] = ginput(1); % [frequency of the right-most selected point, it's magnitude]
    delete(t1)
    jDC = find(Fexc>=DC_low_freq & Fexc<=DC_high_freq); % finds the points in the TF that are in the range of the frequencies manually selected
    DC_gain = mean(Mfrf(jDC)); % Taking the mean to get the mean magnitudes of the points. We know that the DC gain is always:   DCgain = (1/Beta)*(kcyt0/(kcyt0 + kcyt2 + ktrap)), one can derive this.

    % Plotting a few extra things to guide the user (pole and zero locations and such)
    vline_y_z_frf = [min(Mfrf) max(Mfrf)];
    vline_x_z_frf = [ig.fz_frf ig.fz_frf];
    vline_y_p_frf = [min(Mfrf) max(Mfrf)];
    vline_x_p_frf = [ig.fp_frf ig.fp_frf];
    plot(vline_x_z_frf,vline_y_z_frf,':r','linewidth',4), plot(vline_x_p_frf,vline_y_p_frf,':r','linewidth',4) % plotting vertical lines to indicate where the zero and pole were selected
    figure, loglog(FP,PYm,'b','linewidth',2),hold on, hps=gcf;
    set(gca,'Ylim',[1e-14 1]);
    xlabel('f (Hz)'), ylabel('Power Spectrum (V^2*s)')
    vline_y_p_ps = logspace(-10,-6,2); 
    vline_x_p_ps = [ig.fp_frf ig.fp_frf];
    plot(vline_x_p_ps,vline_y_p_ps,':r','linewidth',4); % plotting a vertical line in the power spectrum
    freqplot = logspace(-2,4,350); % this makes 350 log-evenly spaced points, from 10^-2 to 10^4, and is used for plotting the functions. 
    y = 0.0001*freqplot.^-2;
    loglog(freqplot,y,'m','linewidth',1) % plotting a slope of 2 in the power spectrum. Slope smaller than 2 means a subdiffusive movement of the bead.

    % Finds the slope of the power spectrum past the pole corner frequency, and also defines (based on user input click) the filtering corner frequency (falias) 
    disp('Click on the range of frequency to calculate slope. The max freq. will also be used as the maximum freq. limit of fitting')
    disp('The min freq. for the power spectrum fitting is selected later')
    txt = 'Click on the range of frequency to calculate slope. The max freq. will also be used as the maximum freq. limit of fitting';
    t1 = text(0.02,10*max(PYm),0,txt,'FontSize',9); 
    txt = 'The min freq. for the power spectrum fitting is predetermined with the variable fmin';
    t1_p2 = text(0.02,5*max(PYm),0,txt,'FontSize',9);
    [slope_psmin,ymin] = ginput(1);
    [slope_psmax,ymax] = ginput(1);
    delete(t1),delete(t1_p2)
    jps = find(FP>=slope_psmin & FP<=slope_psmax); % finds the index that are in the range of frequencies selected to calculate the slope.
    pslope = polyfit(log(2*pi.*FP(jps)),log(PYm(jps)),1); % finds the slope using a first order polynomial (y = mx+b). Multiply by 2*pi to convert from Hz to rad/s. 
    % NOTE: polyfit gives the coefficients of the polynomial of degree N. In this case, the degree is 1, i.e. mx + b. With m = pslope(1), b = pslope(2)
    loglog(FP,exp(polyval(pslope,log(2*pi.*FP))),'r') % plotting the slope on the power spectrum graph
    % NOTE: polyval evaluates a polynomial at the given x-value. In this case, it finds the power spectral density value at each frequency following the slope and initial value.
    fmax = slope_psmax;
     disp('Select min freq for power spectrum fitting')
    [fmin,ymin] = ginput(1); 
    if abs(pslope(1)) > 2
        disp('---------------------     Warning!     ---------------------------')
        disp('Filtering may be occuring in the range where the slope was selected as the magnitude of the slope is greater than 2')
        disp('------------------------------------------------------------------') 
        txt = 'NOTE:'; 
        t1 = text(0.02,1*max(PYm),0,txt,'FontSize',16,'Color','r'); 
        txt = 'Filtering may be occuring in the range where the slope was selected';
        t1_p2 = text(0.02,0.33*max(PYm),0,txt,'FontSize',11,'Color','r');
        txt = 'as the magnitude of the slope is greater than 2';
        t1_p3 = text(0.02,0.13*max(PYm),0,txt,'FontSize',11,'Color','r');   
    end
    
else % in that case, there is a minimal amount of values that the user has to click
    
    figure, loglog(FP,PYm,'b','linewidth',2),hold on, hps=gcf;
    xlabel('f (Hz)'), ylabel('Power Spectrum (V^2*s)')
    txt = 'Click on the MIN then MAX frequency to define the frequency range to fit';
    t1 = text(0.02,10*max(PYm),0,txt,'FontSize',9); 
    txt = 'The min freq. for the power spectrum fitting is predetermined with the variable fmin';
    t1_p2 = text(0.02,5*max(PYm),0,txt,'FontSize',9); 
    disp('Click on the MIN then MAX frequency to define the frequency range to fit')
    %disp('The min freq. for the power spectrum fitting is predetermined with the variable fmin')
    [fmin,ymin] = ginput(1); 
    [fmax,ymax] = ginput(1); 
    delete(t1), delete(t1_p2)
    freqplot = logspace(-2,4,350); % this makes 350 log-evenly spaced points, from 10^-2 to 10^4, and is used for plotting the functions. 
end

fminsies = ['-fmin',num2str(fmin)];

disp('Click on corner frequency of QPD response (i.e. filtering frequency due to the time response of the QPD) ' )
disp('If there is no filtering visible, click as far as possible to the right (outside of plot box)' )
txt = 'Click on corner frequency of QPD response (i.e. filtering frequency due to the time response of the QPD) ';
t1 = text(0.02,10*max(PYm),0,txt,'FontSize',9); 
txt = 'If there is no filtering visible, click as far as possible to the right (outside of plot box)';
t1_p2 = text(0.02,5*max(PYm),0,txt,'FontSize',9); 
[falias,yfalias]=ginput(1);%6e3;%20e3; %falias:  QPD has filtering effect. first order low pass filter
delete(t1),delete(t1_p2)

if multi == 1
    faliasglob = falias;
end

% Spikes in the PS to be removed for thermal calibration (both random parts noise peaks and also the active response from the input oscillations. 
% It woud be nice to have this automated, or semi-automated. Or, one could use the true passive spectrum by taking a ~20s measurement without laser excitation.
% However, in theory, it should not matter as in the frequency domain, the baseline IS the passive spectrum.
%Fexc_PS_reject=[Fexcall,235:250,650:730,260:280,905:915,330:395,660:720,760:800,845,686,382,923.5,6172,7445,7490,8000,8800,2950,810:830];%600:760];%1700:2200,2920,2815,5250,6500:6580,7600:7650];%,random_noise,random_noise2,random_noise4,random_noise9,2815,3853,random_noise10,random_noise11];%random_noise2,random_noise8,random_noise5];


% Ben 2023 check that the CTRAP frequencies are not getting rejected. 
Fexc_PS_reject=[Fexcall,55:57,70:71,85:180,235:250,650:730,260:280,905:915,330:395,660:720,760:800,845,686,382,923.5,6172,7445,7490,8000,8800,2950,810:830];%600:760];%1700:2200,2920,2815,5250,6500:6580,7600:7650];%,random_noise,random_noise2,random_noise4,random_noise9,2815,3853,random_noise10,random_noise11];%random_noise2,random_noise8,random_noise5];

jfit1 = find(FP>fmin & FP<fmax & isnan(PYm)~=1); % selecting the frequency range INDEX in POWER SPECTRUM to be used for fitting 
for ke = 1:numel(Fexc_PS_reject)
   je = find(abs(FP(jfit1)-Fexc_PS_reject(ke))>0.01*Fexc_PS_reject(ke)); % cutting out the frequency peaks that we don't want to fit
   jfit1=jfit1(je);
end
fdata1 = FP(jfit1); % selecting the frequency points in the PS 
ydata1 = cumtrapz(FP(jfit1),PYm(jfit1));%interp1(fx,IPSx,fdata,'pchip') % integrating the power spectrum
ydat10 = ydata1; % I think this is to save a copy of the selected data.
figure, semilogx(fdata1,ydata1,'b','linewidth',2), hips=gcf; % this plots the integrated power spectrum, in the frequency range being fitted
xlabel('f (Hz)'), ylabel('$\int P_{th}$ ($V^2$)'), axis tight
figure(hps), semilogx(fdata1,PYm(jfit1),'y') % this plots the PS points that are used for the power spectrum fit, in yellow 

if multi == 1
    ydat10glob = ydat10';
end

%%%%%%%%%%%%%%%% INITIAL GUESS VALUES
% (BEN 2023) YOU ARE ONLING USING THIS FIRST HARDCODED SECTION
if semi_auto_ig == 0 && load_file_ig == 0 && user_ig == 1 % this will use the initial guess manually entered below
    
    % - Ben G. perform calibration with the pylake software and provide
    % accurate hardcoded guesses.
    paramset = 'hardcoded-guess'; %add guess type into file name (-Ora)

    ig.k_trap = ktrap_guess; %0.1; %0.1; % manual input
    ig.Beta = beta_guess; %550; %1000; %500; % manual input
    ig.Alpha = alpha_guess; %0.75; % manual input
    ig.GammaR = gammar_guess; %1.68; % manual input, this means the drag coeff. of the cytoplasm is "ig.GammaR" times the one of water
    ig.k_cyt0 = kcyt0_guess; %0.02; % manual input
    ig.k_cyt2 = kcyt2_guess; % manual input - (BEN july 2023)
    ig.k_cyt1 = kcyt1_guess; %0.002; % manual input
    
elseif semi_auto_ig == 0 && load_file_ig == 1 && user_ig == 0 % this will use the initial guess loaded from a previous file
    
    paramset = 'previous-fit'; %add guess type into file name (-Ora)
    
    cd(datdir) % sets the file directory to the current one, to make sure the .mat to load is in the same file.
    load('u2_OLD_saved.mat')
    u0_loaded = u2_OLD; % saving the loaded variables into a new variable. 
    ig.Beta = u0_loaded(1);
    ig.GammaR = u0_loaded(2);
    ig.Alpha = u0_loaded(3);
    ig.k_trap = u0_loaded(4);
    ig.k_cyt0 = u0_loaded(5);
    ig.k_cyt1 = u0_loaded(6);
    
elseif semi_auto_ig == 1 && load_file_ig == 0 && user_ig == 0  % this will use a semi-automated initial guess setting (to be more finely tuned, August 2019, Loïc). 
    
    paramset = 'manual-selection'; %add guess type into file name (-Ora)
    
    ig.Beta = ((ig.Pz)*(-198) - 11378 + (ig.MinfdB)*(-264) - 15595)/2; % phenomenological expression
    ig.k_cyt1 = ((ig.Pz)*(0.000261) + 0.0196 + (ig.MinfdB)*(0.000241)  + 0.0178)/2; % phenomenological expression
    ig.k_cyt0 = ((ig.Pz)*(0.0040) + 0.3057 + (ig.MinfdB)*(0.0033) + 0.2459)/2; % phenomenological expression
    ig.k_trap = (ig.k_cyt0)/(10^(DC_gain/20)*ig.Beta)-(ig.k_cyt0); % theoretical expression
    ig.GammaR = 0.5; %1.68; % manual input, this means the drag coeff. of the cytoplasm is "ig.GammaR" times %ORAORAORA
    % note about alpha: it is not very accurate to pick alpha based on the slope of the power spectrum after the corner frequency. I think it's because the power spectrum is often not clean enough to manually pick a proper range for the slope. 
    if abs(pslope(1)) < alphaUB+1 && abs(pslope(1)) > alphaLB+1 % this means alpha is between 0.5 and 1.0, which is what it's supposed to be if there is no filtering (roll down of power spectrum)
        ig.Alpha = abs(pslope(1))-1; % alpha + 1 = slope of the power spectrum after the corner frequency. i.e., alpha = slope - 1. This alpha is the same as in MSD
    elseif abs(pslope(1)) > alphaUB+1  % this means alpha is greater than 1.0, which suggests filtering and the alpha is set to the mean alpha from previous experiments, and "falias" (which is in stricly speaking, f_filter) will take care of adding a filtering effect 
        %ig.Alpha = 0.5639;
        ig.Alpha = 0.8*alphaUB; % this should be set to the upper limit, i.e. 1.0, but I find that 1.0 is too much.
    elseif abs(pslope(1)) < alphaLB+1 % this means alpha is smaller than 0.5, which suggests aliasing (flattening of power spectrum). The alpha is set to the minimum physical value it can take
        ig.Alpha = alphaLB;
    end
end

%  [Beta      GammaR           Alpha            k_trap          k_cyt0           k_cyt2             k_cyt1               mass    nu]; 
u0=[ig.Beta, ig.GammaR,       ig.Alpha,         ig.k_trap,      ig.k_cyt0,       ig.k_cyt2,         ig.k_cyt1,           10,     1];   % (BEN july 2023) added kcyt_2 % setting up tight upper and lower limit for this first fit
lb=[500      0.1*ig.GammaR    0.95*ig.Alpha     0.1*ig.k_trap   0.1*ig.k_cyt0    0.1*ig.k_cyt2      0.9*ig.k_cyt1        1e-1    0]; % (BEN) using the same bounds as with the k_cyt0
ub=[2500     10*ig.GammaR     1.05*ig.Alpha     3*ig.k_trap     10.*ig.k_cyt0    10.*ig.k_cyt2      1.1*ig.k_cyt1        1e1     1.01];
%test_ben_2023 = length(u0);
%disp(u0(9))
% This is to make sure alpha is never below 0.5 or above 1.0
if u0(3) < alphaLB
    u0(3) = alphaLB;
elseif u0(3) > alphaUB
    u0(3) = alphaUB;
end

% This is to make sure the bounds on alpha is never below 0.5 or above 1.0
if lb(3) < alphaLB
    lb(3) = alphaLB;
end
if ub(3) > alphaUB
    ub(3) = alphaUB;
end

% this plots the initial fit in black
[mf0,phf0] = frf_r9_2bead(u0,Fexc,Rbead,Zbead,kT,falias); 
figure(hfrf), subplot(311)  %set(gca,'Xlim',[0.01 5e3]), set(gca,'Ylim',[-80 -70]),hold on,
semilogx(Fexc,20.*log10(mf0),'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
subplot(312), hold on
semilogx(Fexc,phf0,'k','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
figure(hps), hold on
psf0 = ps_r6(u0,freqplot,fsamp,Rbead,Zbead,kT,falias);
loglog(freqplot,psf0,'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')

    %setting up some fit stuff
    if zph == 0  % setting the proper matrix size for the case where the phase is fitted or not (one can go around this by using weights to control the fitting, instead of truncating/resizing/et...!)
        PHfrf = [];
        xdata = [Fexc,fdata1'];
        ydata = [zeros(size(H)),zeros(size(ydat10'))];
    else
        xdata = [Fexc,Fexc,fdata1'];
        ydata = [zeros(size(H)),zeros(size(H)),zeros(size(ydat10'))];
    end

    W_vector = ones(1,length(Fexc)); % weight multiplier applied to each FRF (mag and phase) points. Can change this if you want different weight applied to different points.
    if Wlast_points == 1
        W_vector(end) = 2;
        %W_vector(end-1) = 2;
    end
    W = Wmag.*W_vector;
    Wphase = Wphase.*W_vector;
    
    if multi == 1
        Wglob = W;
        Wphaseglob = Wphase;
    end
    
if multi == 1
    
    %set fitting options
    options = optimoptions('lsqcurvefit', 'TolFun', 1e-15, 'TolX', 1e-15, 'StepTolerance', 1e-15, 'MaxFunEvals',1e4, 'MaxIter', 1e3);%, 'DiffMinChange', 1e-5);

    %set lower and upper bound scaling
    lbmult = 0.01;
    ubmult = 10; 

    %set lower and upper bounds
    lb = [u0(1:2).*lbmult, 0.5, u0(4:7).*lbmult]; %(BEN) - changed index 6 to 7 to account for added parameter
    ub = [u0(1:2).*ubmult, 1, u0(4:7).*ubmult];

    problem = createOptimProblem('lsqcurvefit','objective',@theor_trap_frf_multi,'x0',u0,'xdata',xdata,'ydata',ydata,'lb',lb,'ub',ub, 'options', options);

    ms = MultiStart('PlotFcns',@gsplotbestf);
    [paramsList(2,:),resNormList(2,:)] = run(ms,problem,50);
    u0 = paramsList(2,:);
end

if usefixnfit == 0 % will fit all parameters simultaneously (-Ora)
    
    fitmethod = 'fit-r10-'; %add fit method to file name (-Ora)
    %%%%%%%%%%%%%%%%%%%%%%%% FIRST FIT %%%%%%%%%%%%%%%%%%%%%%%

    options_fit1 = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-3);
    [u2i,resnorm,residual]=lsqcurvefit(@theor_trap_frf_r9_newTF,u0,xdata,ydata,lb,ub,options_fit1,abs(H(:))',PHfrf,W,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias); 

    % This plots the intermediate fit in green, from the parameters found from the first fit above.
    [mf1,phf1]=frf_r9_2bead(u2i,Fexc,Rbead,Zbead,kT,falias);
    figure(hfrf), subplot(311), hold on
    semilogx(Fexc,20.*log10(mf1),'g','linewidth',2)
    set(gca,'fontsize',14,'fontweight','bold')
    subplot(312), hold on
    semilogx(Fexc,phf1,'g','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
    figure(hps), hold on
    psf1=ps_r6(u2i,freqplot,fsamp,Rbead,Zbead,kT,falias);
    loglog(freqplot,psf1,'g','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
    
    %k_trap=2*pi*(fc_frf-fc_zero)*Gamma;
    mi=u2i(end-1); %*10^-21 g
    nui=u2i(end); %*10^12 nm^2/s
    betai=u2i(1);
    gammari=u2i(2);
    alphai=u2i(3);
    k_trapi=u2i(4);
    k_cyt0i=u2i(5);
    %k_eq0=u2(6);
    %w_0=k_eq0/(gammar*9.42e-6);
    k_cyt1i=u2i(7);%(k_eq0-k_trap-k_cyt0)/(w_0^alpha);
    k_cyt2i = u2i(6); % (BEN 2023) added to accomodate new term

    disp(['INTERMEDIATE FIT --------------------'])
    disp(['Beta = ', num2str(betai), ' nm/V'])
    disp(['Gamma = ',num2str(gammari),' *gamma0'])
    disp(['Alpha = ',num2str(alphai)])
    disp(['k_trap = ',num2str(k_trapi),' pN/nm'])
    disp(['k_cyt,0 = ', num2str(k_cyt0i), ' pN/nm'])
    disp(['k_cyt,1 = ',num2str(k_cyt1i),' pN/nm'])
    disp(['k_cyt2 = ',num2str(k_cyt2i),' pN/nm']) % (BEN 2023) added to accomodate new term
    disp(['Beta*k_trap = ',num2str(betai*k_trapi),' pN/V'])

    %%%%%%%%%%%%%%%%%%%%%%%% SECOND FIT %%%%%%%%%%%%%%%%%%%%%%%
    %%%%% [Beta     GammaR        Alpha          k_trap        k_cyt0        k_cyt2       k_cyt1          mass       nu]  
    lb2=  [200      0.1*u2i(2)    0.90*u2i(3)    0.1*u2i(4)    0.1*u2i(5)    0.1*u2i(6)   0.5*u2i(7)      1e-14      0.01]; % (BEN 2023) added kcyt2 terms - following similar bound to kcyt0   
    ub2=  [5000     10*u2i(2)     1.1*u2i(3)     10*u2i(4)     10*u2i(5)     10*u2i(6)    2*u2i(7)        1e2        1.01];
    W2=W;%.*C(:)';

    % This is to make sure alpha is never below 0.5 or above 1.0
    if lb2(3) < alphaLB
        lb2(3) = alphaLB;
    end
    if ub2(3) > alphaUB
        ub2(3) = alphaUB;
    end

     % This fits using the parameters found from the first fit
    options_fit2 = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-6);
    [u2,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(@theor_trap_frf_gly3_r9_newTF,u2i,xdata,ydata,lb2,ub2,options_fit2,abs(H(:))',PHfrf,W2,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias);
    u2_ci=nlparci(u2,residual,'jacobian',full(jacobian)); %nlparci to calculate confidence intervals on the fit

else % fits parameters iteratively in two subsets (-Ora)
    
    fitmethod = ['fixnfit-r10','-fitOrder',num2str(fitOrder),'-fitGroups',num2str(groupem),'-']; %add fit method to file name (-Ora)
    
    %for 2 holds
    %lbscale and ubscale need to be the same length, with each number repeated at least once per hold
    lbscale = [0.1 0.1 0.1 0.1 0.25 0.25 0.5 0.5 0.75 0.75 0.9 0.9];
    ubscale = [10 10 10 10 4 4 2 2 1.25 1.25 1.1 1.1];

    %store fit parameters
    paramsList(1,:) = [ig.Beta, ig.GammaR, ig.Alpha, ig.k_trap, ig.k_cyt0, ig.k_cyt2, ig.k_cyt1, 10, 1];
    if multi == 1
        paramsList(2,:) = u0;
    end
    
    %set up initial guess
    u0 = paramsList(1+multi,:);
    disp(length(u0));
    % this is where the cod is messing up

    % this plots the initial fit in black
    [mf0,phf0] = frf_r9_2bead(u0,Fexc,Rbead,Zbead,kT,falias); 
    figure(hfrf), subplot(311)  %set(gca,'Xlim',[0.01 5e3]), set(gca,'Ylim',[-80 -70]),hold on,
    semilogx(Fexc,20.*log10(mf0),'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')
    subplot(312), hold on
    semilogx(Fexc,phf0,'k','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
    figure(hps), hold on
    psf0 = ps_r6(u0,freqplot,fsamp,Rbead,Zbead,kT,falias);
    loglog(freqplot,psf0,'k','linewidth',2),set(gca,'fontsize',14,'fontweight','bold')


    % The fitting
    for fitNum = 1:length(lbscale)

        %set the initial guess
        u0 = paramsList(fitNum+multi,:); %use last fit's parameters as initial guess

        lbmult = lbscale(fitNum); %set this iteration's lower bound
        ubmult = ubscale(fitNum); %set this iteration's upper bound
        lbfix = 1;
        ubfix = 1;
        
        if groupem == 1
             if mod(fitNum,2) == fitOrder
            %vary ktrap, kcyt0, kcyt2, kcyt1 
            lb = [u0(1:3).*lbfix, u0(4:7).*lbmult, 0, 0];   %  (BEN 2023) changed the index to account for added ter,
            ub = [u0(1:3).*ubfix, u0(4:7).*ubmult, 100, 100];
            else
            %vary beta, gamma, alpha 
            lb = [u0(1:3).*lbmult, u0(4:7).*lbfix, 0, 0];
            ub = [u0(1:3).*ubmult, u0(4:7).*ubfix, 100, 100];
            end
        elseif groupem == 2
            if mod(fitNum,2) == fitOrder
            %vary ktrap, beta
            lb = [u0(1).*lbmult, u0(2:3).*lbfix, u0(4).*lbmult, u0(5:7).*lbfix, 0, 0]; % (BEN 2023) - reindexded to account for new parameter
            ub = [u0(1).*ubmult, u0(2:3).*ubfix, u0(4).*ubmult, u0(5:7).*ubfix, 100, 100];
            else
            %vary kcyt0, kcyt2, kcyt1, gamma, alpha 
            lb = [u0(1).*lbfix, u0(2:3).*lbmult, u0(4).*lbfix, u0(5:7).*lbmult, 0, 0];
            ub = [u0(1).*ubfix, u0(2:3).*ubmult, u0(4).*ubfix, u0(5:7).*ubmult, 100, 100];
            end
        else
            disp('You''ve got problems, my friend...')
        end

        % This is to make sure alpha is never below 0.5 or above 1.0
        if u0(3) < alphaLB
            u0(3) = alphaLB;
        elseif u0(3) > alphaUB
            u0(3) = alphaUB;
        end

        % This is to make sure the bounds on alpha is never below 0.5 or above 1.0
        if lb(3) < alphaLB
            lb(3) = alphaLB;
        end
        if ub(3) > alphaUB
            ub(3) = alphaUB;
        end

        options_fit1 = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-3);
        [paramsList(fitNum+1+multi,:), resNormList(fitNum+1+multi), resList(fitNum+1+multi,:)]=lsqcurvefit(@theor_trap_frf_gly3_r9_newTF,u0,xdata,ydata,lb,ub,options_fit1,abs(H(:))',PHfrf,W,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias); 

    end

    %find fit with minumum residual
    bestFit = find(resNormList == min(resNormList(2:end))); %get the best fit out of all iterations
    resFit = resNormList(bestFit); %residual of best fit

    bestFit = bestFit(1) %in case multiple fits have the same residual, select the first
    resFit = resFit(1)

    W2 = W;
    %%%%%%%%%%%%%%%%%%%%%%%% BEST FIT %%%%%%%%%%%%%%%%%%%%%%%

    % This fits the final parameters to get jacobian etc
    lb2 = paramsList(bestFit,:);
    ub2 = paramsList(bestFit,:);
    options_fit2 = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',1e5,'MaxIter',2e3);%,'DiffMinChange',1e-6);
    [u2,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(@theor_trap_frf_gly3_r9_newTF,paramsList(bestFit,:),xdata,ydata,lb2,ub2,options_fit2,abs(H(:))',PHfrf,W2,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias);
    u2_ci=nlparci(u2,residual,'jacobian',full(jacobian)); %nlparci to calculate confidence intervals on the fit
end 

errfit = theor_trap_frf_gly3_r9_newTF(u2,xdata,abs(H(:))',PHfrf,W2,Wphase,ydat10',fsamp,Rbead,Zbead,kT,falias); % this function gives the weighted residual
jfrf = 1:numel(Fexc);
if isempty(PHfrf) == 1
    jph=[];
    %PHfrf_fit=zeros(size(Fexc));
    jps=(jfrf(end)+1):numel(errfit);
else
    jph=jfrf+numel(Fexc);
    %PHfrf_fit=((errfit(jph).*10./W)+1).*PHfrf;
    jps=(jph(end)+1):numel(errfit);
end

[Mfrf_fit,PHfrf_fit]=frf_r9_2bead(u2,Fexc,Rbead,Zbead,kT,falias);

%PHfrf_fit = abs(PHfrf_fit);% this makes the angles all positive... why do we do this?


%Mfrf_fit=((errfit(jfrf)./W)+1).*abs(H(:))';
IPS_fit=(errfit(jps)+1).*ydat10';  % this takes the integrated power spectrum weighted residual and transforms it into non-weighted non-residual 
PS_fit=ps_r6(u2,freqplot,fsamp,Rbead,Zbead,kT,falias);
    
% unweighted residual = data - fit
Res_mag = abs(H) - Mfrf_fit;
Res_phase = PHfrf - PHfrf_fit;
Res_IPS = ydat10' - IPS_fit;

% Meas squared error (for the user to compare fitness, when it's very difficult to decide between two fits)
MSE_mag = (1/length(Res_mag))*sum(Res_mag.^2);
MSE_phase = (1/length(Res_phase))*sum(Res_phase.^2);
MSE_IPS = (1/length(Res_IPS))*sum(Res_IPS.^2);
   

% This plots the final fit in cyan, from the parameters found from the second (no fixnfit) or best (fixnfit) fit above.
figure(hfrf), subplot(311), hold on,
set(gca,'fontsize',14,'fontweight','bold')
semilogx(Fexc,20.*log10(Mfrf_fit),'c','linewidth',2) %%%%%%%
set(gca,'Xlim',[0.9*x_lower 1.11*x_upper]), set(gca,'Ylim',[1.01*y_lower 0.99*y_upper]); % this is where the zoom in happens
%set(gca,'Xlim',[10e-3 10e3]), set(gca,'Ylim',[-77 -68]);
title(['FRF with mag MSE ' num2str(MSE_mag) ' and phase MSE ' num2str(MSE_phase) ' '])
subplot(312), hold on,
semilogx(Fexc,PHfrf_fit,'c','linewidth',2)%+(180/pi).*unwrap(angle(exp(-j*1e-6.*Hw)))
set(gca,'fontsize',14,'fontweight','bold')

figure(hps), hold on
loglog(freqplot,PS_fit,'c','linewidth',2)
set(gca,'Xlim',[1e-1 1e4]), set(gca,'Ylim',[1e-14 1e-4]); % this is where the zoom in happens
title(['PS with integrated PS MSE ' num2str(MSE_IPS) ''])


% IPth1=cumtrapz(fdata1,Pth_t0);%excited_power_spectrum(u,fx(jfit),falias,Fexc,Adrive1,ydat0,f0);
% figure(hips), hold on, semilogx(fdata1,IPth1,'c','linewidth',2)
figure(hips), hold on, semilogx(fdata1,IPS_fit,'c','linewidth',2)

P12=20.*log10(Pyy)-20.*log10(Pyy2);
figure, semilogx(Fexc,P12,'bo')
xlabel('Frequency (Hz)'), ylabel('Pyy(Fexc) - Pyy(2*Fexc)')
set(gca,'fontsize',14,'fontweight','bold')


%k_trap=2*pi*(fc_frf-fc_zero)*Gamma;
m=u2(end-1); %*10^-21 g
nu=u2(end); %*10^12 nm^2/s
beta=u2(1);
gammar=u2(2);
alpha=u2(3);
k_trap=u2(4);
k_cyt0=u2(5);
%k_eq0=u2(6);
%w_0=k_eq0/(gammar*9.42e-6);
k_cyt1=u2(7);%(k_eq0-k_trap-k_cyt0)/(w_0^alpha);
k_cyt2=u2(6); %(BEN 2023) added

disp(['FINAL FIT -------------------------------'])
disp(['Beta = ', num2str(beta), ' nm/V  +/- ' num2str(u2(1)-u2_ci(1,1)) ''])
disp(['Gamma = ',num2str(gammar),' *gamma0  +/- ' num2str(u2(2)-u2_ci(2,1)) ''])
disp(['Alpha = ',num2str(alpha), '  +/- ' num2str(u2(3)-u2_ci(3,1)) ''])
disp(['k_trap = ',num2str(k_trap),' pN/nm  +/- ' num2str(u2(4)-u2_ci(4,1)) ''])
disp(['k_cyt,0 = ', num2str(k_cyt0), ' pN/nm  +/- ' num2str(u2(5)-u2_ci(5,1)) ''])
disp(['k_cyt,1 = ',num2str(k_cyt1),' pN/nm  +/- ' num2str(u2(7)-u2_ci(7,1)) ''])
disp(['k_cyt,2 = ',num2str(k_cyt2),' pN/nm  +/- ' num2str(u2(6)-u2_ci(6,1)) ''])
disp(['Beta*k_trap = ',num2str(beta*k_trap),' pN/V'])

forexcel = [beta,u2(1)-u2_ci(1,1),gammar,u2(2)-u2_ci(2,1),alpha,u2(3)-u2_ci(3,1),k_trap,u2(4)-u2_ci(4,1),k_cyt0,u2(5)-u2_ci(5,1),k_cyt1,u2(6)-u2_ci(6,1)];

k_cyt_tot=k_cyt0+k_cyt1.*(2*pi.*FP).^alpha; %Question - why isn't this alpha-1?
Gp=k_cyt_tot./(6*pi*Rbead);
figure, loglog(FP,Gp,'linewidth',2)
xlabel('f (Hz)'), ylabel('Gp ($pN/nm^2$)')

ww=2*pi.*FP;
ig.Gamma0 = 3.2577e-06; % manual input, drag coefficient of water at 37 degrees [pN*s/nm]

g=k_cyt0 + k_cyt1.*((j.*ww).^alpha)./gamma(alpha)...
    + gammar*ig.Gamma0.*j.*ww;
Gp=g./(6*pi.*Rbead);
v=imag(g)./ww; % this is the drag coefficient - not the dynamic viscosity - the line below fixes that.
v = v.*(1000000/(3000*3.14159265)); 
figure, hold on
loglog(FP,real(Gp)*1e6,'color','b','linewidth',2)
ylabel('Gp - Storage Modulus')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')
%hold on % see if we can get them to plot on different graphs.
loglog(FP,imag(Gp)*1e6,'color','b','linewidth',2)
ylabel('Gpp - Loss Modulus')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')
figure, hold on
loglog(FP,v,'color','b','linewidth',2)
ylabel('viscosity')
set(gca,'Xscale','log'), set(gca,'Yscale','log')
set(gca,'fontsize',14,'fontweight','bold')


% this saves the viscoelastic fit
simplesave = input('Save simple fit and exit now? ');
if simplesave == 1

%matfilename = strrep(filename,'.mat','-'); % this replaces the ".mat" of the original file name, to "-"
%savename=string([matfilename,'OTC-nosqwv-','.mat']);
%save(savename)
save(filename);

return
end
