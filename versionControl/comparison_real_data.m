% Compare multiple noise reduction methods
close all; clear; clc;
%======================%
noiseReductionMethods = {'Boll', 'Berouti',...
    'MultiBandKamath', 'MMSESTSA', 'Wiener',...
    'WaveletAllFrames', 'WaveletStudentT'};

noiseReductionMethods = {'WaveletStudentT'};
%=========Whistle=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Eastern Whipbird\01.wav';

%=========Click=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Eastern Yellow Robin\02.wav';

%=========Slurs=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Little Friarbird\04.wav';

%=========Warble=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Indian Peafowl\04.wav';

%=========Block=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Australian Masked Owl\03.wav';

%=========Stacked Harmonic=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Australian Owlet_Nightjar\01.wav';

%=========Oscillation=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Rufous Whistler\07.wav';

% audioPath = 'C:\Jie\work\data\NIPS2014\NIPS4B_BIRD_CHALLENGE_TRAIN_TEST_WAV\train\nips4b_birds_trainfile001.wav';
% [cleanSignal, fs] = audioread(audioPath);

% fig = figure('rend','painters','pos',[300 300 300 450]);
% [I1,F1,T1] = wav_to_spec(cleanSignal(:,1), fs, 512, 0.5);
% imagesc(I1); axis xy;
% title('Eastern Yellow Robin', 'FontSize',14);
% xlabel('Time index', 'FontSize',12);
% ylabel('Frequency bin index', 'FontSize',12);
% print(fig,'Click','-dpng')

cleanPath = 'C:\Jie\work\data\NIPS2014\Noise Reduction clean\nips4b_birds_trainfile001.wav';
[cleanSignal,fs] = audioread(cleanPath);

noisedPath = 'C:\Jie\work\data\NIPS2014\Noise Reduction noisy\nips4b_birds_trainfile001.wav';
[noisedSignal,fs] = audioread(noisedPath);

%======================%
nMethods = length(noiseReductionMethods);
outSNR = cell(1, nMethods);
for iMethod = 1:nMethods
    
    disp(iMethod);
    %======================%
    % [cleanSignal, fs] = audioread(audioPath);
    % cleanSignal = cleanSignal(:,1); % mono to stero
    % S1 = mean(cleanSignal.^2); %// actual signal power
    
    % normSignal = cleanSignal / sum(cleanSignal);
    % S2 = mean(normSignal.^2); %// actual signal power
    
    % ADD noise of specified SNR
    
    % add white noise
    %wavin3 = load('.\noise\white.mat'); % load babble noise
    %wavin3 = wavin3.white;
    %[noisedSignal,~] = addnoise(cleanSignal, wavin3(1:length(cleanSignal)), snrValue);
    
    % add white noise
    %noisedSignal_white = awgn(cleanSignal, snrValue, 'measured');
    
    awgnNoise = cleanSignal - noisedSignal;
    
    snr_value_1 = GetSNR(cleanSignal, awgnNoise);
    snr_value_2 = snr(cleanSignal, awgnNoise);
    
    %======================%
    %noisedSignal = noisedSignal(:,1); % mono to stero
    
    % noiseReductionMethod = 'MultiBandKamath';
    noiseReductionMethod = noiseReductionMethods{iMethod};
    
    IS = 0.1; % seconds of Initial silence
    
    plot_type = 'no';
    %======================%
    winLenArray = 2.^[8:10];
    %winLenArray = floor(0.1*fs);
    
    nWinLen = length(winLenArray);
    outSNR_array = zeros(nWinLen,1);
    orgSNR_array = zeros(nWinLen,1);
    for iWin = 1:nWinLen
        
        winLen = winLenArray(iWin);
        
        SP = 0.4; % shift percentage
        
        switch noiseReductionMethod
            
            % case 'Bandpass'
            
            case 'Berouti'
                noiseReducedSignal = SSBerouti79(noisedSignal,fs,IS,winLen,SP);
                
            case 'Boll'
                noiseReducedSignal = SSBoll79(noisedSignal,fs,IS,winLen,SP);
                
            case 'Towsey1'
                W = winLen; % samples
                nfft=W;
                wnd=hamming(W);
                
                NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments
                Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
                
                y=segment(noisedSignal,W,SP,wnd);
                Y=fft(y,nfft);
                YPhase=angle(Y(1:floor(end/2)+1,:)); %Noisy Speech Phase
                Y=abs(Y(1:floor(end/2)+1,:)).^Gamma; %Specrogram
                %Y1 = (Y - min(Y(:)))/(max(Y(:))-min(Y(:)));                 
                Y1 = Y;
                
                % do ostu
                [level,EM] = graythresh(Y1);
                
                Y1(Y1 < level) = 0;
                               
                Y2 = wiener2(Y1,[5 5]);

                numberOfFrames=size(Y2,2);
                T1 = length(noisedSignal)/fs;
                rT = max(T1)/numberOfFrames;

                X = TowseyQUT(Y2,rT,IS);
                                   
                % figure;
                % subplot(211); imagesc(X); axis xy;
                % subplot(212); imagesc(Y2); axis xy;
                
                noiseReducedSignal=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(SP*W));
                
            case 'Towsey2'
                W = winLen; % samples
                nfft=W;
                wnd=hamming(W);
                
                NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments
                Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
                
                y=segment(noisedSignal,W,SP,wnd);
                Y=fft(y,nfft);
                YPhase=angle(Y(1:floor(end/2)+1,:)); %Noisy Speech Phase
                Y1=abs(Y(1:floor(end/2)+1,:)).^Gamma; %Specrogram
                
                numberOfFrames=size(Y,2);
                FreqResol=size(Y,1);
                
                Y = 255*(Y - min(Y(:)))/(max(Y(:))-min(Y(:))); 
                
                X = withoutSubbandModeIntensities(Y);
                                
                noiseReducedSignal=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(SP*W));
                
            case 'MultiBandKamath'
                noiseReducedSignal = SSMultibandKamath02(noisedSignal,fs,IS,winLen,SP);
                
            case 'MMSESTSA'
                noiseReducedSignal = MMSESTSA84(noisedSignal,fs,IS,winLen,SP);
                
            case 'Wiener'
                noiseReducedSignal = WienerScalart96(noisedSignal,fs,IS,winLen,SP);
                
            case 'WaveletAllFrames'
                
                % https://ww2.mathworks.cn/help/wavelet/ref/wden.html
                lev = 13;
                %wname = 'sym8';
                wname = 'db10';
                %tempSignal = wden(noisedSignal,'minimaxi','h','mln',lev,wname);
                %tempSignal = wden(noisedSignal,'sqtwolog','h','mln',lev,wname);
                %tempSignal = wden(noisedSignal,'modwtsqtwolog','s','mln',lev,wname);
                tempSignal = wden(noisedSignal,'rigrsure','s','mln',lev,wname);
                
                if size(tempSignal,1) < size(tempSignal,2)
                    noiseReducedSignal = tempSignal';
                else
                    noiseReducedSignal = tempSignal;
                end
                
            case 'WaveletSingleFrame_2'
                
                noiseReducedSignal = WaveletSingleFrame_2(noisedSignal,fs,IS,winLen,SP);
                
            case 'WaveletStudentT'
                
                noiseReducedSignal = WaveletStudentT(noisedSignal,fs,IS,winLen,SP);
                
        end
        
        
        switch plot_type
            
            case 'yes'
                
                % Display images
                %======================%
                % figure;
                % subplot(211);
                % orgSNR = snr(orgSignal, fs);
                % plot(orgSignal);
                % title(['SNR of Original Signal: ', num2str(orgSNR)])
                %
                % subplot(212);
                % % bollSNR = snr(noiseReducedSignal);
                % bollSNR = snr(noiseReducedSignal, fs);
                % plot(noiseReducedSignal);
                % title(['SNR after Boll: ', num2str(bollSNR)])
                %======================%
                figure;
                subplot(311);
                % orgSNR = snr(orgSignal, fs);
                [I1,F1,T1] = wav_to_spec(cleanSignal, fs, 512, 0.5);
                imagesc(I1); axis xy;
                % title(['Clean Signal: ', num2str(orgSNR)])
                title('Clean Signal')
                
                subplot(312);
                orgSNR = GetSNR(cleanSignal, awgnNoise);
                % orgSNR = snr(cleanSignal, awgnNoise);
                [I2,F2,T2] = wav_to_spec(noisedSignal, fs, 512, 0.5);
                imagesc(I2); axis xy;
                title(['SNR of Noisy Signal: ', num2str(orgSNR)])
                
                subplot(313);
                if length(cleanSignal)~= length(noiseReducedSignal) % do zero padding
                    noiseReducedSignal = [noiseReducedSignal;...
                        zeros(length(cleanSignal) - length(noiseReducedSignal), 1)];
                end
                
                % enhancedSNR = snr(cleanSignal, cleanSignal - noiseReducedSignal);
                enhancedSNR = GetSNR(cleanSignal, cleanSignal - noiseReducedSignal);
                
                [I3,F3,T3] = wav_to_spec(noiseReducedSignal, fs, 512, 0.5);
                imagesc(I3); axis xy;
                title(['SNR after Noise Reduction: ', num2str(enhancedSNR)])
                
            case 'no'
                
                orgSNR_array(iWin) = GetSNR(cleanSignal, awgnNoise);
                
                if length(cleanSignal)~= length(noiseReducedSignal) % do zero padding
                    noiseReducedSignal = [noiseReducedSignal;...
                        zeros(length(cleanSignal) - length(noiseReducedSignal), 1)];
                end
                
                outSNR_array(iWin) = GetSNR(cleanSignal, cleanSignal - noiseReducedSignal);
        end
    end
    %======================%    
    % figure;
    % plot(diffSNR_array);
    % ylabel('SNR improvement (dB)');
    % title(noiseReductionMethod);
    
    outSNR{iMethod} = outSNR_array;
    
end

outSNRmat = cell2mat(outSNR);

plotPattern = {'-o', '-^', '-x', '-s', '-h', '-v', '->'};

if length(outSNRmat) == 1
    figure;
    plot(outSNRmat, 'LineWidth',1);
    xticklabels(noiseReductionMethods);
    xtickangle(45);
else
    fig = figure;
    [M,N] = size(outSNRmat);
    for i = 1:N
        plot(outSNRmat(:,i), plotPattern{i}, 'LineWidth',2);
        hold on;
    end
    legend(noiseReductionMethods);
end
title('nips4b_birds_trainfile001')
xlabel('Input SNR (dB)');
% xticklabels({'-10', '', '-5', '', '0', '', '5'});
xlabel('Window size (sample)');
xticklabels({'256', '', '512', '', '1024'});
ylabel('Output SNR (dB)');

print(fig,'nips4b_birds_trainfile001','-dpng')

%[EOF]
