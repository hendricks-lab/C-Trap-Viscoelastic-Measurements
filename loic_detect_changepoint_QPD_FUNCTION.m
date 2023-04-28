%%% Loïc Chaubet
%%% March 2019
%%%
%%%
%%% Cropping QPD signal that comes from manual drift correction during OT
%%% experiments, or any other sudden and persistant movement of the bead
%%% with respect to the trap. This code takes a moving average of the 
%%% standard deviation of xQUAD and xAOD, and looks at the ratio of 
%%% squared. This way, any abornmal signal seem to show up quite easily. 
%%% Then, to select which peaks are abnormal, a threshold value is compared
%%% to the standard dev. of the ratio. The peak heights and positions are obtained, 
%%% along with the width of the peaks. The width are somewhat dependant on 
%%% the window used for detection (which is set to be equal to the width of
%%% the moving window (mov_wind). Then, an integer amount of excitation 
%%% periods (1/fexc) is selected to cover the width of the peak. The periods 
%%% are also centered so that the period that contains the position
%%% of the peak is exactly centered with respect to the latter, while the 
%%% other periods are distributed equally on the left and right side of the
%%% peak. 
%%%
%%% NOTE: the locations of the peaks are indepedant of fexc, but the peaks
%%% that will be cut out of the data are dependant on fexc. 


function [indices_to_cut,all_peaks_loc,Npeaks_cut,peaks_cut] = loic_detect_changepoint_QPD_FUNCTION(xQUAD,xAOD,fexc,fs,mov_wind,std_th,cut_th,plot_cut)



xQUAD_mstd = movstd(xQUAD,mov_wind);
xAOD_mstd = movstd(xAOD,mov_wind);
ratio = xQUAD_mstd./xAOD_mstd;
n = 1:length(xQUAD);
dratio2 = (ratio-mean(ratio)).^2;
std_threshold = std_th*std(dratio2);
std_threshold_plot = [];
std_threshold_plot(1:length(ratio)) = std_threshold;





[pks,all_peaks_loc,w] = findpeaks(dratio2,'MinPeakHeight',std_threshold,'MinPeakDistance',mov_wind); % gives the peak value, the peak locations and the peak widths in INDICES

if isempty(pks) == 1 % if empty --> 1
    indices_to_cut = [];
    all_peaks_loc = [];
    Npeaks_cut = [];
    peaks_cut = [];
    disp('No peaks detected above threshold')
    return
elseif isempty(fexc) ==1 % this is for the quick check
    indices_to_cut = [];
    Npeaks_cut = [];
    peaks_cut = [];
    return
end
    
Texc = (1/fexc)*fs; % period of excitation [# samples]

if plot_cut == 1
    figure,
    plot(n,xQUAD), hold on
    for k = 1:length(all_peaks_loc)

    range_location_line_plot = 10*mean(xQUAD);
    location_line_Yplot = -range_location_line_plot:0.1*range_location_line_plot:range_location_line_plot;
    location_line_Xplot(1:length(location_line_Yplot)) = all_peaks_loc(k);
    plot(location_line_Xplot,location_line_Yplot,'k')

    end
    xQUAD_signal=gcf;
end

if plot_cut == 1
    figure, 
    plot(n,dratio2,'b',n,std_threshold_plot,'k'), hold on
    plot(all_peaks_loc,pks,'mv','MarkerSize',12,'MarkerFaceColor','m')
    std_ratio=gcf;
end

%highlights the peaks in red
if plot_cut == 1
    for k = 1:length(all_peaks_loc)
    w_plot = dratio2(all_peaks_loc(k)-w(k):all_peaks_loc(k)+w(k));
    X = all_peaks_loc(k)-w(k):all_peaks_loc(k)+w(k);
    plot(X,w_plot,'r')
    end

    figure(xQUAD_signal)
    for k = 1:length(all_peaks_loc)
    w_plot = xQUAD(all_peaks_loc(k)-w(k):all_peaks_loc(k)+w(k));
    X = all_peaks_loc(k)-w(k):all_peaks_loc(k)+w(k);
    plot(X,w_plot,'r')
    end
end

if plot_cut == 1
    figure(std_ratio),hold on
    title(['Fexc = ' num2str(fexc) '']) 
end
% centers the excitation period with the peak, and calculates the largest
% amount of periods that can fit into the width of the peak. Highlights the
% parts to cut in green
for j = 1:length(all_peaks_loc)
    n_cycles = 2*w(j)/Texc; % how many times can a full period [# samples] can be cut in twice of the width [# samples] of the peak. Twice because that way it seems to cover the full width of the peak better.
    if n_cycles < 1 && n_cycles > cut_th % even though the width of the peak is smaller than the width of one period, the cutting still happens as the ratio of width/period is still not too small
       w_allowed = Texc; % actual width allowed to be cut (smaller than the full width of the peak)
       index_centered_cut = all_peaks_loc(j)-floor(w_allowed/2):all_peaks_loc(j)+floor(w_allowed/2); % this makes sure that the indices are centered about the peak, and equally distributed left and right
       if max(index_centered_cut) > length(xQUAD) % this means the peak is close to the end of the data
           index_centered_cut = min(index_centered_cut):length(xQUAD); % this sets the index to be cut to be from the left of the peak until the end of the data set
       elseif min(index_centered_cut) < 0 % this means the peak is close to the beginning of the data
           index_centered_cut = 1:max(index_centered_cut); % this sets the index to be cut to be from the beginning of the data until the right of the peak
       end
       w_centered = dratio2(index_centered_cut); 
       xx{j} = index_centered_cut; % x indices of the data that must be cut
       
       if plot_cut == 1
           plot(xx{j},w_centered,'g')
       end
    elseif n_cycles < 1 && n_cycles < cut_th % now the ratio of width/period is way too small, meaning that if a cut is done at the minimum lenght of one period, too much data will be removed
       xx{j} = []; % this makes sure that nothing is cut if the conditions above are met
    else % in that case, more than one period fits in 2x width of the peak
       n_cycles = round(n_cycles);
       w_allowed = n_cycles*Texc; % actual width allowed to be cut (smaller than the full width of the peak)
       index_centered_cut = all_peaks_loc(j)-floor(w_allowed/2):all_peaks_loc(j)+floor(w_allowed/2);
       if max(index_centered_cut) > length(xQUAD) % this means the peak is close to the end of the data
           index_centered_cut = min(index_centered_cut):length(xQUAD); % this sets the index to be cut to be from the left of the peak until the end of the data set
       elseif min(index_centered_cut) < 0 % this means the peak is close to the beginning of the data
           index_centered_cut = 1:max(index_centered_cut); % this sets the index to be cut to be from the beginning of the data until the right of the peak    
       end
       w_centered = dratio2(index_centered_cut); % 
       xx{j} = index_centered_cut; % x indices of the data that must be cut, centered to have half of the allowed width on one side, the other half on the other side
       if plot_cut == 1
           plot(xx{j},w_centered,'g')
       end
    end
    peaks_cut(j) = ~isempty(xx{j}); % if it's not empty, i.e. if the peak was cut, this outputs "1" 
    Npeaks_cut = nnz(peaks_cut);
end

%title(['Moving average standard deviation ratio with window size = ' num2str(mov_wind) ' samples and standard devation threshold = ' num2str(std_th) ''])
%legend('Ratio squared','Std threshold','Peaks')
indices_to_cut = [];

% places all the indices in one line vector
for i = 1:length(xx)
    indices_to_cut= [indices_to_cut,xx{i}];
end

if plot_cut == 1
    figure(xQUAD_signal)
    title(['Fexc = ' num2str(fexc) '']) 
end

for j = 1:length(all_peaks_loc)
    n_cycles = 2*w(j)/Texc; % how many times can a full period be cut in the total width (2x w) of the peak
    if n_cycles < 1 && n_cycles > cut_th
       w_allowed = Texc; % actual width allowed to be cut (smaller than the full width of the peak)
       index_centered_cut = all_peaks_loc(j)-floor(w_allowed/2):all_peaks_loc(j)+floor(w_allowed/2);
       if max(index_centered_cut) > length(xQUAD) % this means the peak is close to the end of the sample
           index_centered_cut = min(index_centered_cut):length(xQUAD); % this sets the index to be cut to be from the left of the peak until the end of the data set
       end
       w_centered = xQUAD(index_centered_cut); % 
       xx{j} = index_centered_cut; % x indices of the data that must be cut, centered to have half of the allowed width on one side, the other half on the other side
       if plot_cut == 1
           plot(xx{j},w_centered,'g')
       end
    elseif n_cycles < 1 && n_cycles < cut_th % now the ratio of width/period is way too small, meaning that if a cut is done, too much of the good data will also be removed
       xx{j} = [];
    else
       n_cycles = round(n_cycles);
       w_allowed = n_cycles*Texc; % actual width allowed to be cut (smaller than the full width of the peak)
      index_centered_cut = all_peaks_loc(j)-floor(w_allowed/2):all_peaks_loc(j)+floor(w_allowed/2);
       if max(index_centered_cut) > length(xQUAD) % this means the peak is close to the end of the sample
           index_centered_cut = min(index_centered_cut):length(xQUAD); % this sets the index to be cut to be from the left of the peak until the end of the data set
       end
       w_centered = xQUAD(index_centered_cut); % 
       xx{j} = index_centered_cut; % x indices of the data that must be cut, centered to have half of the allowed width on one side, the other half on the other side
       if plot_cut == 1
           plot(xx{j},w_centered,'g')
       end
    end
    peaks_cut(j) = ~isempty(xx{j}); % if it's not empty, i.e. if the peak was cut, this outputs "1" 
    Npeaks_cut = nnz(peaks_cut); % finds the number of non-zero elements (i.e. number of cuts)
end

%Npeaks = length(pks);
%figure(xQUAD_signal)
%plot(n(indices_to_cut),xQUAD(indices_to_cut),'c')
%L = length(n)
%n(indices_to_cut) = [];
%L_cut = length(n)
%R = L_cut/L

%end


