%--%
close all; clear; clc;

% apply Wiener filter to spectrogram for enhancing bioacoustic recordings.
% audioPath = 'D:\Project_NoiseReduction\data\noisy_data\Bush Stone-curlew\se_399_188294_20101015_022030_4_852999999999156.wav';

filename = {'.\data\noisy_data\Bush Stone-curlew\se_399_188294_20101015_022030_4_852999999999156.wav';
    '.\data\noisy_data\Grey Fantail\se_399_188293_20101014_104704_3_940999999998894.wav';
    '.\data\noisy_data\Southern Boobook\sc_site_1_577_230339_20111011_041359_22_951999999999316.wav';
    '.\data\noisy_data\Tawny Grassbird\s4_413_190155_20100816_072939_5_017000000007101.wav'};

nFile =length(filename);
for iFile = 1:nFile
    
    filePath = filename{iFile};
    disp(filePath)
    
    [initialSignal, fs] = audioread(filePath);
    
    %--noisy signal is manually selected
    IS = 0.25; % seconds of Initial silence
    
    denoiseMethods ={ 'wiener1D',   'SSBoll',  'SSBerouti', 'SSTowsey', 'MMSESTSA' };
    
    nMethods = length(denoiseMethods);
    for iMethod  = 1:nMethods
        denoiseMethod = denoiseMethods{iMethod};
        
        %--parameter setting
        param.winLen = 512;
        param.SP = 0.4;
        
        switch denoiseMethod
            
            case 'MMSESTSA'
                enhancedSignal = MMSESTSA84(initialSignal, fs, IS, param);
                
            case 'SSBerouti'
                enhancedSignal = SSBerouti79(initialSignal, fs, IS, param);
                
            case 'SSBoll'
                enhancedSignal = SSBoll79(initialSignal,fs,IS,param);
                
            case 'SSTowsey'
                enhancedSignal = SSTowsey(initialSignal,fs,IS,param);
                
            case 'wiener1D'
                % enhancedSignal = wiener2D(initialSignal,fs,IS,param);
                enhancedSignal = WienerScalart96(initialSignal,fs,IS,param);
                
        end
        
        %--evalution metrics
        minLen = min(length(initialSignal), length(enhancedSignal));
        initialSignal = initialSignal(1:minLen);
        enhancedSignal = enhancedSignal(1:minLen);
        
        initialNoise = initialSignal(1:floor(0.25*fs));
        enhancedNoise = enhancedSignal(1:floor(0.25*fs));
        
        if iFile == 1
            powerInitialSignal = sum(initialSignal(8820:24696).^2);
            powerEnhancedSignal = sum(enhancedSignal(8820:24696).^2);
        elseif iFile == 2
            powerInitialSignal = sum(initialSignal(108486:123039).^2);
            powerEnhancedSignal = sum(enhancedSignal(108486:123039).^2);
        elseif iFile == 3
            powerInitialSignal = sum(initialSignal(72352:106253).^2);
            powerEnhancedSignal = sum(enhancedSignal(72352:106253).^2);
        else
            powerInitialSignal = sum(initialSignal(72168:215049).^2);
            powerEnhancedSignal = sum(enhancedSignal(72168:215049).^2);
        end
        
        powerInitialNoise =  sum(initialNoise.^2);
        PowerEnhancedNoise =  sum(enhancedNoise.^2);
        
        %--evaluation metrics
        SnNR_noisy(iFile, iMethod) = 10 * log10(powerInitialSignal / powerInitialNoise);
        SnNR_denoise(iFile, iMethod) = 10 * log10(powerEnhancedSignal / PowerEnhancedNoise);
        success_ratio(iFile, iMethod) = log10( var(initialNoise) / var(enhancedNoise));
        PSNR(iFile, iMethod) = 20* log10(max(powerInitialSignal) / sqrt(calMSE(initialSignal, enhancedSignal)));
        
    end
end

