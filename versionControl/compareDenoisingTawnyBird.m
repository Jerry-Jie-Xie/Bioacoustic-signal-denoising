%--%
close all; clear; clc;

% apply Wiener filter to spectrogram for enhancing bioacoustic recordings.
audioPath = '.\data\noisy_data\Tawny Grassbird\s4_413_190155_20100816_072939_5_017000000007101.wav';

[initialSignal, fs] = audioread(audioPath);

%--noisy signal is manually selected
IS = 0.15; % seconds of Initial silence

denoiseMethods ={ 'wiener2D',   'SSBoll',  'SSBerouti', 'SSTowsey', 'MMSESTSA' };
nMethods = length(denoiseMethods);
for iMethod  = 1:nMethods
    denoiseMethod = denoiseMethods{iMethod};
    
    %--parameter setting
    param.winLen = 1024;
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
            
        case 'wiener2D'
            enhancedSignal = wiener2D(initialSignal,fs,IS,param);
            
        case 'waveletDenoisingFrame'
            
            Y = calculateSpectrogram(initialSignal, param);
            figure; imagesc(Y); axis xy; colorbar;
            title('Initial Spectrogram');
            enhancedSignal = waveletDenoisingFrame(initialSignal,fs,IS,param);
            X = calculateSpectrogram(enhancedSignal, param);
            figure; imagesc(X); axis xy; colorbar;
            title('Denoised Spectrogram');
            
        case 'waveletDenoising'
            Y = calculateSpectrogram(initialSignal, param);
            figure; imagesc(Y); axis xy; colorbar;
            title('Initial Spectrogram');
            wavelet_level = 8;
            wavelet_type = 'db2';
            enhancedSignal = wden(initialSignal,'sqtwolog','s','mln', wavelet_level, wavelet_type);
            %enhancedSignal = wden(initialSignal,'modwtsqtwolog','s','mln', wavelet_level, wavelet_type);
            X = calculateSpectrogram(enhancedSignal, param);
            figure; imagesc(X); axis xy; colorbar;
            title('Denoised Spectrogram');
            
        otherwise
            disp('WRONG method!')
            
    end
    
    %--evalution metrics
    minLen = min(length(initialSignal), length(enhancedSignal));
    initialSignal = initialSignal(1:minLen);
    enhancedSignal = enhancedSignal(1:minLen);
    
    initialNoise = initialSignal(1:floor(0.15*fs));
    enhancedNoise = enhancedSignal(1:floor(0.15*fs));
    
    powerInitialSignal = sum(initialSignal(72168:215049).^2);
    powerEnhancedSignal = sum(enhancedSignal(72168:215049).^2);
    powerInitialNoise =  sum(initialNoise.^2);
    PowerEnhancedNoise =  sum(enhancedNoise.^2);
    
    if size(enhancedSignal, 2) > 1
        enhancedSignal = enhancedSignal';
    end
    
    %--evaluation metrics
    SnNR_noisy(iMethod) = 10 * log10(powerInitialSignal / powerInitialNoise);
    SnNR_denoise(iMethod) = 10 * log10(powerEnhancedSignal / PowerEnhancedNoise);
    success_ratio(iMethod) = log10( var(initialNoise) / var(enhancedNoise));
    PSNR(iMethod) = 20* log10(max(powerInitialSignal) / sqrt(calMSE(initialSignal, enhancedSignal)));
    
end
