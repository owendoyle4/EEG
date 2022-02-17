function [results] = getFFTpwelchFooof(params ,FFTBoolean, pwelchBoolean, fooofBoolean)
    %%%%%%%%
    % INPUTS
    % params:
    %   IDs: A list of subject IDs who will get incorporated into the
    %   subsequent calculations.
    %   condition: Frequency of AM, 10 or 20.
    %   freqRangeName: A string of the frequency range (Theta, Alpha, or Beta).
    % 
    % FFTBoolean: If true, return the spectrum and x_axis (frequencies) for 
    % the FFT as calculated by plotFFT.
    % 
    % pwelchBoolean: If true, return the pwelch spectra for individual
    % blocks and the average of the pre blocks, along with the
    % corresponding x_axis (frequencies).
    % 
    % fooofBoolean: If true, return the original, fooofed, aperiodic, and
    % aperiodic-corrected power spectra along with the corresponding x_axis
    % (frequencies).
    % 
    % powerBoolean: If true, return the spectral power over the bandwidth associated with
    % freqRangeName. The power is the integral under the fooofed, aperiodic-corrected pwelch spectrum.
    %
    %
    % OUTPUTS
    %
    % results: Includes a variety of spectra and frequency vectors based on
    % the input booleans.
    %%%%%%%%
    
    addpath('/Users/owen/Library/CloudStorage/Box-Box/kTMP_EEG/code/analysis/BrainCap/individualFreqAnalysisFuncs');
    addpath('/Users/owen/Library/CloudStorage/Box-Box/kTMP_EEG/code/');
    addpath('/Users/owen/Library/CloudStorage/Box-Box/kTMP_EEG/code/fooof/');
    
    IDs = params.IDs;
    condition = params.condition;
    blocks = params.blocks; 
    
    numIDs = length(IDs);
    numBlocks = length(blocks);
    numChans = 31; % exclude nose reference
    
    %% set necessary variables depending on booleans
    
    % load sample block data that is used for initilizing dimensions
    if FFTBoolean || pwelchBoolean
        cd(strcat('/Users/owen/Library/CloudStorage/Box-Box/kTMP_EEG/SDATA/06RO'));
        fileName = '06RO_10hz_SDATA.mat';
        load(fileName);
        
        SDATA.data = bsxfun(@minus,SDATA.data(:,1:end),SDATA.data(:,strcmp(SDATA.info.channel_labels, 'ECG'))); % nose ref
        SDATA.events.event_times = trimData(SDATA, 4, 'first');
        blockTimes = SDATA.events.event_times(1,:);
        sr = SDATA.info.sampling_rate;
        blockData = SDATA.data(blockTimes(1):blockTimes(2),:);
    end
        
    
    if FFTBoolean
        [FFTspectrum, ~] = plotFFT(blockData,sr);
        FFTspectra = zeros(numChans,numBlocks,numIDs,size(FFTspectrum,1));
    end
        
    
    if pwelchBoolean
        % Resolution Required 
        resolutonReqd = 0.1;
        % Calculate number of FFT points (NFFT) required for the required resolution fs/NFFT = resolutionReqd
        NFFT = sr / resolutonReqd;
        numAverages = 30;
        aveBlockLen = sr*60*4;
        window = aveBlockLen/numAverages; % This makes number of FFT averages to X samples/ window = numAverages
        overlap = 0;
        
        [s,x_axis] = pwelch(blockData,window,overlap,NFFT,sr);
        
        PSD = zeros(numChans,numBlocks,numIDs,size(s,1));
        x_axes_temp = zeros(numBlocks,numIDs,size(s,1));
        blockAvePSDPreBlocks = zeros(numChans,numIDs,size(s,1));

        
        if fooofBoolean
            settings = struct();
            settings.peak_width_limits = [5 20];
            settings.max_n_peaks = 4;
            settings.min_peak_height = 0.01;
            settings.peak_threshold = 0.1;
            settings.aperiodic_mode = 'fixed';
            settings.verbose = true;

            f_range = [2 40];
            return_model = true;

            fooof_results = fooof(x_axis, s(:,1), f_range, settings, return_model);

            ap_corrected = zeros(numIDs,numChans,numBlocks,length(fooof_results.power_spectrum));
            fooof_power_spectrum = zeros(numIDs,numChans,numBlocks,length(fooof_results.power_spectrum));
            fooofed_spectrum = zeros(numIDs,numChans,numBlocks,length(fooof_results.power_spectrum));
            ap_fit = zeros(numIDs,numChans,numBlocks,length(fooof_results.power_spectrum));

            ap_corrected_avePreBlocks = zeros(numIDs,numChans,length(fooof_results.power_spectrum));
            x_axis_fooof = fooof_results.freqs;
        end
    end
    
%     clear FFTspectrum s x_axis fooof_results ap_corrected;
    
	%% load in SDATA
    for IDidx = 1:length(IDs)
        tic;
        ID = IDs(IDidx);
        % cd to folder
        cd(strcat('/Users/owen/Library/CloudStorage/Box-Box/kTMP_EEG/SDATA/',ID));

        % load SDATA
        if exist(strcat(ID,'_',string(condition),'hz_SDATA_interp_grant.mat'), 'file')
            fileName = strcat(ID,'_',string(condition),'hz_SDATA_interp_grant.mat');
        elseif exist(strcat(ID,'_',string(condition),'hz_SDATA_interp.mat'), 'file')
            fileName = strcat(ID,'_',string(condition),'hz_SDATA_interp.mat');
        else
            fileName = strcat(ID,'_',string(condition),'hz_SDATA.mat');
        end

        load(fileName);
        disp(fileName);

        % nose reference
        SDATA.data = bsxfun(@minus,SDATA.data(:,1:end),SDATA.data(:,strcmp(SDATA.info.channel_labels, 'ECG')));
        
        % trim data to be a uniform 4 minutes, does not remove artifacts
%         SDATA.events.event_times = trimData(SDATA, 4, 'last');
        
        chans = ~strcmp(SDATA.info.channel_labels, 'ECG');
        sr = SDATA.info.sampling_rate;
        
        for blockIdx = 1:length(blocks)
            block = blocks(blockIdx);
            disp(block);   
           
           % trim block's data to remove artifacts
           blockData = trimBlock(SDATA, block, chans, 4, 'last', true);

            if FFTBoolean
                [spectrumFFT,x_axisFFT] = plotFFT(blockData, sr);
                FFTspectra(:,blockIdx,IDidx,:) = spectrumFFT';
            end

            if pwelchBoolean
                [spectrumPwelch,x_axisPwelch] = pwelch(blockData,window,overlap,NFFT,sr);
                PSD(:,blockIdx,IDidx,:) = spectrumPwelch';
                x_axes_temp = x_axisPwelch * (sr/(2*pi)); 
            end


            if fooofBoolean
                for chan = 1:numChans
                    power_spectrum = spectrumPwelch(:,chan); 
                    freqs = x_axisPwelch;
                    fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
                    fooof_power_spectrum(IDidx,chan,blockIdx,:) = fooof_results.power_spectrum;
                    fooofed_spectrum(IDidx,chan,blockIdx,:) = fooof_results.fooofed_spectrum;
                    ap_fit(IDidx,chan,blockIdx,:) = fooof_results.ap_fit;
                    ap_corrected(IDidx,chan,blockIdx,:) = fooof_results.fooofed_spectrum-fooof_results.ap_fit;
                end
            end
            
        end
        
        % average the preblocks together
        if pwelchBoolean
            patPSDPreBlocks = squeeze(PSD(:,1:2,IDidx,:));
            blockAvePSDPreBlocks(:,IDidx,:) = squeeze(mean(patPSDPreBlocks,2));
        end
        
        if fooofBoolean
            pat_ap_correctedPreBlocks = squeeze(ap_corrected(IDidx,:,1:2,:));
            ap_corrected_avePreBlocks(IDidx,:,:) = squeeze(mean(pat_ap_correctedPreBlocks,2));
        end
        toc;
    end
    
    
    %% store everything in structure
    results = struct();
    results.chans = SDATA.info.channel_labels(chans);
   
    if FFTBoolean
        results.fft_spectra = FFTspectra;
        results.fft_x_axis = x_axisFFT;
    end
    
    if pwelchBoolean
        results.pwelch_spectra = PSD;
        results.pwelchs_PreBlocksAve = blockAvePSDPreBlocks;
        results.pwelch_x_axis = x_axisPwelch;
    end
    
    if fooofBoolean
        results.fooof_power_spectrum = fooof_power_spectrum;
        results.fooofed_spectrum = fooofed_spectrum;
        results.ap_fit = ap_fit;
        results.ap_corrected = ap_corrected;
        results.ap_corrected_PreBlocksAve = ap_corrected_avePreBlocks;
        results.fooof_x_axis = x_axis_fooof;
        results.fooof_settings = settings;
    end
end
