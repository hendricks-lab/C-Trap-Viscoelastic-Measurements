%%% Loïc Chaubet
%%% March 2019
%%% This function calculates the power spectrum, FRF, and coherence,
%%% as well as the added functionality to cut out any part of the photodiode 
%%% signal that shows abrupt and persistant change (likely from sudden stage
%%% movement when using the joystick during drift correction, for
%%% example).

%%% NOTICE HERE THAT the mean sumQUAD signal in volts output is an
%%% APPROXIMATION because the sum signal channel length does not get
%%% updated as compared to the xQUAD and xAOD channel. Similarily, the
%%% power spectrum graph is also directly from the raw signal, and thus is
%%% from a slighlty larger sample than the FRF (when cutting is enabled)
%%% This should not matter as the corner frequency should not be affected.

function [H,Fexc,C,FP,PYm,Pyy,Pyy2,fsamp,length_used_out,sumQUAD_mean,Npeaks_out,all_loc_out,K_int_out,frequency_range_tolerance,NFFT,number_windows,data_lost,Noverlap]= Step1_TF_Loic_2018b_FUNCTION(N_channels,kpd,AODx_input,QUADx_input,nfftcoef,hannw,plot_range,cut_peaks_data,mean_threshold,cut_threshold,mov_wind_s,plot_cut,Fexc_in,fsamp_in,cx,fileID,manualcut,manualcut_data_end)
    fsamp = fsamp_in;
   
    if N_channels == 6
        %%%% Reading from text file (needs to be in the current folder)
        rawdata = textscan(fileID, '%f %f %f %f %f %f', 'HeaderLines', 0);
    elseif N_channels == 4
        rawdata = textscan(fileID, '%f %f %f %f', 'HeaderLines', 0);
    end
    
    if manualcut == 1 % this changes the raw data so that only a truncated part of the data is used
        rawdata{QUADx_input} = rawdata{QUADx_input}(1:end-manualcut_data_end*fsamp);
        rawdata{3} = rawdata{3}(1:end-manualcut_data_end*fsamp); 
        rawdata{AODx_input} = rawdata{AODx_input}(1:end-manualcut_data_end*fsamp);
    end
    
    xQUAD = rawdata{QUADx_input}.*kpd; % [V]
    sumQUAD_mean = mean(rawdata{3}); % [V] % just to give an idea of the sum.
    xAOD = (rawdata{AODx_input}-75).*cx; % [nm]
    
   

    Ldata = length(xAOD);
    Nfft_pre = 2^(nextpow2(Ldata)-1); %computationally optimized to get the 
    %fastest fourier transform. Nfft should never exceed the length of the
    %data, otherwise zero padding happens.       

    Nfft = floor(Nfft_pre/fsamp_in)*fsamp_in; % 
    % Nfft = number of Fast Fourier Transform points for the whole signal.

    % large Nfft = high frequency resolution (i.e. small bin size), less
    % averaging (power spectrum very coarse), but response could leak to adjacent 
    % bins causing lower magnitudes of peaks. 

    % small Nfft = low frequency resolution (i.e. bing bin size) more averaging
    % (smooth PS), but coherence could be low because of the averaging now
    % happening over a very large frequency bin.

    % Nfft also needs to be an integer multiple of the excitation period when
    % used for the transfer function. Here, for thermal response, it does not
    % matter and it's better to use the definition above (or is it?).

    Nwindow = Nfft; % Here we arbitrarily take the window to be the same as Nfft for simplicity. Note, if we take Nwindow = Nfft/10, then we get very similar (not quite exactly the same) response as if we did Nfft = Nfft/10, and Nwindow = Nfft.            
    Noverlap = floor(Nwindow/2); % floor(Nfft_pre/2); % very standard to have 50% overalp. This could be fine tuned to maximize response, especially at lower frequencies where there are very little full cycles (i.e. Nfft is LARGE)
    Nint = 10; % number of points to integrate
    n = 1:Ldata;
    ctime=n'./fsamp_in; % constructing the time vector and flipping the row-colum to match the data

    [FP,PYm] = power_spectrum_Adam_Loic_FUNCTION([ctime,xQUAD]',Nwindow,Noverlap,Nfft,Nint,fsamp_in);
    
     
    % this first call to the function is to check if there are peaks at all.
     if cut_peaks_data == 1
         mov_wind = mov_wind_s*fsamp_in;
         [a1,a2,a3,a4] = loic_detect_changepoint_QPD_FUNCTION(xQUAD,xAOD,[],[],mov_wind,mean_threshold,[],0);
     end
     
     
     for k = 1:length(Fexc_in) % looping through the frequencies of excitation
        k
        xQUAD = rawdata{QUADx_input}.*kpd; % [V]
        xAOD = (rawdata{AODx_input}-75).*cx; % [nm]          
        
        if cut_peaks_data == 1
            if isempty(a2) == 1 % if there are no peaks, then the function stops and the code stops trying to cut
                cut_peaks_data = 0;
            else % if there are peaks, the code finds the peaks and cuts accordingly 
                [indices_to_cut,all_loc,Npeaks(k),loc_peaks(k,:)] = loic_detect_changepoint_QPD_FUNCTION(xQUAD,xAOD,Fexc_in(k),fsamp_in,mov_wind,mean_threshold,cut_threshold,plot_cut); % Note here that "peaks" should be the same for all frequencies as it is independant of centering. 
                xQUAD(indices_to_cut) = [];
                xAOD(indices_to_cut) = [];
            end  
        end
        
        L = length(xAOD); % re-storing the proper data length
        data_time = L/fsamp
        [K_int(k),NFFT(k)] = NFFT_generator_loic_v4_constant_windows(Fexc_in(k),fsamp_in,data_time,nfftcoef);
       
        
        if data_time > 230
            if Fexc_in(k) < 0.060 % for f = 0.023, and f = 0.053
                K_overlap(k) = 0.85;
            elseif Fexc_in(k) < 0.11 && Fexc_in(k) > 0.08 % for f = 0.097
                K_overlap(k) = 0.80; % 0.70
            elseif Fexc_in(k) < 0.21 && Fexc_in(k) > 0.14 % for f = 0.17
                K_overlap(k) = 0.75; % 0.60
            elseif Fexc_in(k) < 0.40 && Fexc_in(k) > 0.22 % for f = 0.27
                K_overlap(k) = 0.70; % 0.50 
            elseif Fexc_in(k) < 0.60 && Fexc_in(k) > 0.40 % for f = 0.43
                K_overlap(k) = 0.65; % 0.50
            else
                K_overlap(k) = 0.50;
            end
            
        else
            
            if Fexc_in(k) < 0.060 % for f = 0.023, and f = 0.053
                K_overlap(k) = 0.90;
            elseif Fexc_in(k) < 0.11 && Fexc_in(k) > 0.08 % for f = 0.097
                K_overlap(k) = 0.80; % 0.70
            elseif Fexc_in(k) < 0.21 && Fexc_in(k) > 0.14 % for f = 0.17
                K_overlap(k) = 0.75; % 0.60
            elseif Fexc_in(k) < 0.40 && Fexc_in(k) > 0.22 % for f = 0.27
                K_overlap(k) = 0.70; % 0.50 
            elseif Fexc_in(k) < 0.60 && Fexc_in(k) > 0.40 % for f = 0.43
                K_overlap(k) = 0.65; % 0.50
            else
                K_overlap(k) = 0.50;
            end
        end
        % Changing the % overlap at the lowest frequencies as it seems to help get a better/more realistic response (i.e. no wierd flat coherence at 1.0) 
%         if Fexc_in(k) < 0.060 % for f = 0.023, and f = 0.053
%             K_overlap = 0.85;
%         elseif Fexc_in(k) < 0.11 && Fexc_in(k) > 0.08 % for f = 0.097
%             K_overlap = 0.85; % 0.70
%         elseif Fexc_in(k) < 0.21 && Fexc_in(k) > 0.14 % for f = 0.17
%             K_overlap = 0.80; % 0.60
%         elseif Fexc_in(k) < 0.40 && Fexc_in(k) > 0.22 % for f = 0.27
%             K_overlap = 0.75; % 0.50 
%         elseif Fexc_in(k) < 0.60 && Fexc_in(k) > 0.40 % for f = 0.43
%             K_overlap = 0.65; % 0.50
%         else
%             K_overlap = 0.50;
%         end
% 
%         if Fexc_in(k) < 0.060 % for f = 0.023, and f = 0.053
%             K_overlap = 0.85;
%         elseif Fexc_in(k) < 0.11 && Fexc_in(k) > 0.08 % for f = 0.097
%             K_overlap = 0.80; % 0.70
%         elseif Fexc_in(k) < 0.21 && Fexc_in(k) > 0.14 % for f = 0.17
%             K_overlap = 0.75; % 0.60
%         elseif Fexc_in(k) < 0.40 && Fexc_in(k) > 0.22 % for f = 0.27
%             K_overlap = 0.70; % 0.50 
%         elseif Fexc_in(k) < 0.60 && Fexc_in(k) > 0.40 % for f = 0.43
%             K_overlap = 0.65; % 0.50
%         else
%             K_overlap = 0.50;
%         end

 
        Noverlap(k) = floor(K_overlap(k)*NFFT(k)); % overlap now in [samples]
        
        if NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 2;
            data_lost(k) = (L - NFFT(k) - Noverlap(k) + NFFT(k))/L; % loss in percent
        else
            number_windows(k) = 1;
            data_lost(k) = (L - NFFT(k))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 3;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 4;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 5;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 6;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 7;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 8;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;
        end
        
        if NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) <= L
            number_windows(k) = 9;
            data_lost(k) = (L - (NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k) - Noverlap(k) + NFFT(k)))/L;  
        end
        
        
        
        if hannw == 1
            wind = hann(NFFT(k));
        else
            wind = hamming(NFFT(k));
        end
        % Using Hanning window is standard and should work most of the time.
        % Hann window is better at distinguishing Fexc from a broad range of 
        % frequencies around it (little side lobes), while Hamming is better 
        % at distinguishing Fexc from a narrow range of frequencies. Since 
        % we know our input Fexc, we can purposely look near that frequency
        % using the Pxx (PSD of input) and easily find where the maximum 
        % response is. For the same resolution bandiwdth (i.e. minimum 
        % frequency separation we can detect), the Hann window requires a 
        % slightly longer (10%) signal compared to Hamming. In brief, the 
        % Hamming function retains more information at each window, due to 
        % the boundaries not forced to zero, at the cost of more effects 
        % from discountinuities. However, it using the Hamming window seems 
        % to yield a better response through more averaging, and better
        % coherence.
        %
        % see: https://www.tek.com/blog/window-functions-spectrum-analyzers
        % Hamming seems BETTER!
        
         
        % At this point, xAOD and xQUAD have either been left unmodified or
        % have been slightly truncated at parts due to the peaks being cut
        % out. Note that the data being cut out should be an integer multiple 
        % of the period of excitation of interest. Now to further minimize
        % effects from discontinuities, the data is truncated so that an
        % integer multiple of the Fexc period fits exactly in the data.
        % That functionality is already implemented in the Welch's method
        % from MATLAB (see reference for pwelch)
        %
        % "     If you cannot divide the length of x exactly into 
        %       an integer number of segments with 50% overlap, x 
        %       is truncated accordingly.                           "
       
        length_used(k) = length(xQUAD);
      

        [Txy,FTF] = tfestimate(xAOD,xQUAD,wind,Noverlap(k),NFFT(k),fsamp_in); %Txy is the quotient of the Cross Power Spectral Density of X and Y (Pxy) and the Power Spectral Density of X (Pxx), i.e. Pxy/Pxx.
        cxy = mscohere(xAOD,xQUAD,wind,Noverlap(k),NFFT(k),fsamp_in); % Estimating coherence using magnitude-squared corelation between input X and output Y
        Pyy_p = pwelch(xQUAD,wind,Noverlap(k),NFFT(k),fsamp_in); % PSD of the output (xQUAD)
        Pxx = pwelch(xAOD,wind,Noverlap(k),NFFT(k),fsamp_in); % PSD of the input (xAOD)

        fresolution = (FTF(end)-FTF(1))/length(FTF);
        fresolution_multiplier = floor((K_int(k).^(1/3))./3); % this is found empirically, sets the range to detect the signal to very narrow at lower frequencies (which makes sense as the signal/noise is low, because of the low amount of cycles in one data set), and very wide at higher frequencies (also where we have high signal/noise)
        
        if fresolution_multiplier == 0
            fresolution_multiplier = 1;
        end
        
        fnear_index = find(FTF <= Fexc_in(k)+1.05*fresolution_multiplier*fresolution & FTF >= Fexc_in(k)-fresolution_multiplier*fresolution); % index of the frequencies near the frequency of excitation. 1.05 additional multiplier is to make sure that the upper bound is not missed

        % looking directly at the excitation signal (xAOD)
        Pxx_near = Pxx(fnear_index); % Pxx near the frequency of excitation
        Pxx_max_index_max = fnear_index(find(Pxx_near == max(Pxx_near))); % this is the global INDEX value, indicating where Pxx [nm^2/Hz] is max
        Pxx_max_freq_max = FTF(Pxx_max_index_max); % this is the frequency at which the Pxx is max
         

        fnear_broad_index = find(FTF <= Fexc_in(k)+10*fresolution_multiplier*fresolution & FTF >= Fexc_in(k)-10*fresolution_multiplier*fresolution); % index of the frequencies close to the freq. excitation, over a wider range, for plotting purposes
        
        % this is to check if the coherence is a flat line at 1.0 (which
        % means something is off, either NFFT is too high, or the %overlap is not high enough, or both)
        plot_range_once = 0;
        if (1 - mean(cxy(fnear_broad_index))) < 0.001
            plot_range_once = 1;
            disp(' --------------- WARNING --------------- ')
            disp([' Coherence is a flat line at 1.0, for Fexc = ' num2str(Fexc_in(k)) '. Need to adjust NFFT or %overlap'])         
        end
        
        if plot_range_once == 1 || plot_range == 1
            figure, plot(FTF(fnear_broad_index),cxy(fnear_broad_index),'b'), hold on
            plot(FTF(fnear_index),cxy(fnear_index),'g'), plot(FTF(Pxx_max_index_max),cxy(Pxx_max_index_max),'ko')
            C_range_plot=gcf;
            
            figure, plot(FTF(fnear_broad_index),20.*log10(abs(Txy(fnear_broad_index))),'b'), hold on
            plot(FTF(fnear_index),20.*log10(abs(Txy(fnear_index))),'g'), plot(FTF(Pxx_max_index_max),20.*log10(abs(Txy(Pxx_max_index_max))),'ko')
            M_range_plot=gcf;
            
            figure, semilogy(FTF(fnear_broad_index),Pxx(fnear_broad_index),'b'), hold on
            plot(FTF(fnear_index),Pxx_near,'g'), plot(FTF(Pxx_max_index_max),Pxx(Pxx_max_index_max),'ko')
            Pxx_range_plot=gcf;
        end
        
                
        selection_based = 'based on Pxx max magnitude';
        
        if plot_range_once == 1 || plot_range == 1
            
            figure(C_range_plot),% plot(Cmax_freq_max,cxy(Cmax_index_max),'kx'),
            legend('Broad range (10x)','Actual range investigated','Selected value'),title('Coherence')
            title(['Coherence at ' num2str(Fexc_in(k)) ' , maximized at ' num2str(Pxx_max_freq_max) ' , ' num2str(selection_based) ' , used data length = ' num2str(length_used(k)) ' over ' num2str(Ldata) ''])
            
            figure(M_range_plot) %plot(Mmax_freq_max,20.*log10(abs(Txy(Mmax_index_max))),'kx')
            legend('Broad range (10x)','Actual range investigated','Selected value'),title('Magnitude')
            title(['Magnitude at ' num2str(Fexc_in(k)) ' , maximized at ' num2str(Pxx_max_freq_max) ' , ' num2str(selection_based) ' , used data length = ' num2str(length_used(k)) ' over ' num2str(Ldata) ''])
            
            figure(Pxx_range_plot), %plot(Mmax_freq_max,abs(20.*log10(abs(Txy(Mmax_index_max)))),'kx')
            legend('Broad range (10x)','Actual range investigated','Local max detected'),title('Pxx (PSD of AOD signal)')
            title(['Input PSD at ' num2str(Fexc_in(k)) ' maximized at ' num2str(Pxx_max_freq_max) ' , used data length = ' num2str(length_used(k)) ' over ' num2str(Ldata) ''])      
        end
        
        
        
         if length(Pxx_max_index_max) ~= 1 % this happens in the case where the maxium happens at more than one index (length > 1), or that there are none (length = 0)
             disp('There are two or more maximum coherence/magnitude index, or none. Change the range tolerance to find a proper maximum, ')
             disp('make sure the frequency of excitation vector is the right one')
             disp('and make sure data channels chosen correspond to the actual data channels')
             return
         end
        
        Fexc(k) = FTF(Pxx_max_index_max);
        H(k) = Txy(Pxx_max_index_max);
        Mfrf(k) = 20.*log10(abs(H(k)));
        C(k) = cxy(Pxx_max_index_max);
        PHfrf(k) = unwrap(angle(Txy(Pxx_max_index_max)));
        PHfrf(k) = abs(PHfrf(k) - round(PHfrf(k)./pi).*pi)*180/pi;% why??! 
        Pyy(k) = Pyy_p(Pxx_max_index_max);
        Fexc2x_index = find(abs(FTF-2*Fexc(k)) == min(abs(FTF-2*Fexc(k)))); % getting the index of the power at 2xFexc, why???!!
        Pyy2(k)=Pyy_p(Fexc2x_index); % power at 2xFexc, why??
    end
    
    
    if cut_peaks_data == 1
        Npeaks_out = Npeaks;
        all_loc_out = all_loc;
    else
        Npeaks_out = [];
        all_loc_out = [];
    end
    K_int_out = K_int;
    frequency_range_tolerance = floor((K_int.^(1/3)));
    length_used_out = length_used;
    
end


        
