clc; close all; clear;
%Compile with MATLAB Wavelet Toolbox 3.0
%
%De-Noising Speech Signal Using Discrete Wavelet Transform
%
%The de-noising procedure proceeds in three steps:
%[1]
%Decomposition. Choose a wavelet, and choose a level N. Compute the wavelet
%decomposition of the signal s at level N.
%[2]
%Detail coefficients thresholding. For each level from 1 to N, select
%a threshold and apply soft thresholding to the detail coefficients.
%[3]
%Reconstruction. Compute wavelet reconstruction based on the original
%approximation coefficients of level N and the modified detail coefficients
%of levels from 1 to N.
%
%Programmed by Student of Lambung Mangkurat University
%Nadya Amalia

%load example of speech sound
cleanPath = 'C:\Jie\work\data\NIPS2014\Noise Reduction clean\nips4b_birds_trainfile001.wav';
[truesignal,Fs] = audioread(cleanPath);
noisedPath = 'C:\Jie\work\data\NIPS2014\Noise Reduction noisy\nips4b_birds_trainfile001.wav';
[truesignalN,Fs] = audioread(noisedPath);

% cleanPath = 'C:\Jie\work\data\NIPS2014\Noise Reduction noisy\nips4b_birds_trainfile001.wav';
% [truesignal,Fs] = audioread(cleanPath);

% % snrValue = snrValueArray(iSnr);
% snrValue = -10;
% % add white noise
% wavin3 = load('.\noise\white.mat'); % load babble noise
% wavin3 = wavin3.white;
% [truesignalN,~] = addnoise(truesignal, wavin3(1:length(truesignal)), snrValue);

k = 1;
%---------------------------%
%           DWT :           %
%   wavelet decomposition   %
%---------------------------%
levelArray = [5:15];
parameter = zeros(3,1000);
% levelArray = 5;
for iLevel = 1:length(levelArray)
    
    level = levelArray(iLevel);
    
    wnameCell = {'db1','db2','db3','db4','db5','db6','db7','db8','db9','db10'};
    for iName = 1:length(wnameCell)
        
        wt = wnameCell{iName};
        
        %---------------------------%
        %        thresholding       %
        %---------------------------%
        %TPTR = 'rigrsure', adaptive threshold selection using principle of Stein's
        %Unbiased Risk Estimate
        %TPTR = 'heursure', heuristic variant of the first option
        %TPTR = 'sqtwolog', threshold is sqrt(2*log(length(X)))
        %TPTR = 'minimaxi', minimax thresholding
        
        thresholdingName = {'heursure','rigrsure','minimaxi','sqtwolog','modwtsqtwolog'};
        for iThresh = 1:length(thresholdingName)
            
            tptr = thresholdingName{iThresh};
            
            %scal = 'one'; % Use model assuming standard Gaussian white noise.
            %scal = 'sln'; % Basic model with unscaled noise
            scal = 'mln'; % Use a level-dependent estimation of the level noise
            denoised = wden(truesignalN,tptr,'s',scal,level,wt);
            
            if size(denoised,1)<size(denoised,2)
                denoised = denoised';
            end
            
            %---------------------------%
            err(k) = max(abs(truesignalN-denoised));
            
            %---------------------------%
            %       compute SNR         %
            %---------------------------%
            %SNR - Signal to Noise Ratio
            % SNR = snr(truesignal,truesignalN);
            NoisySNR = 20*log10(norm(truesignal(:)) / norm (truesignal(:)-truesignalN(:)));
            % SNR = snr(truesignal,denoised);
            %DenoisedSNR(k) = 20*log10(norm(truesignal(:)) / norm (truesignal(:)-denoised(:)));
            DenoisedSNR(k) = GetSNR(truesignal, truesignal - denoised);
            
            parameter(1,k) = iLevel;
            parameter(2,k) = iName;
            parameter(3,k) = iThresh;
            
            k = k + 1;
        end
    end
end

% figure(3)
% subplot(3,1,1); plot(truesignal); title('True Speech Signal');
% xlabel('Samples'); ylabel('Amplitude');
% subplot(3,1,2); plot(truesignalN); title('Noisy Speech Signal');
% xlabel('Samples'); ylabel('Amplitude');
% subplot(3,1,3); plot(denoised); title('De-noised Speech Signal');
% xlabel('Samples'); ylabel('Amplitude');
% figure(4)
% subplot(1,3,1); specgram(truesignal,512,Fs); title('True Speech Signal');
% subplot(1,3,2); specgram(truesignalN,512,Fs); title('Noisy Speech Signal');
% subplot(1,3,3); specgram(denoised,512,Fs); title('De-noised Speech Signal');




