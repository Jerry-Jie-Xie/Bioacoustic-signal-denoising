% Compare multiple noise reduction methods
close all; clear; clc;
%======================%
noiseReductionMethods = {'Berouti', 'Boll',...
    'Towsey', 'MultiBandKamath', 'MMSESTSA', 'Wiener',...
    'Wavelet-Jie'};

% noiseReductionMethods = {'Wavelet-Jie'};
%=========Whistle=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Eastern Whipbird\02.wav';

%=========Click=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Eastern Yellow Robin\02.wav';

%=========Slurs=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Little Friarbird\04.wav';

%=========Warble=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Indian Peafowl\04.wav';

%=========Block=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Australian Masked Owl\03.wav';

%=========Stacked Harmonic=============%
% audioPath = 'C:\Jie\work\data\noiseReduction\Australian Owlet_Nightjar\02.wav';

%=========Oscillation=============%
audioPath = 'C:\Jie\work\data\noiseReduction\Rufous Whistler\07.wav';

%======================%
nMethods = length(noiseReductionMethods);
improvedSNR = cell(1, nMethods);
for iMethod = 1:nMethods
    
    disp(iMethod);
    %======================%
    [cleanSignal, fs] = audioread(audioPath);
    cleanSignal = cleanSignal(:,1); % mono to stero
    S1 = mean(cleanSignal.^2); %// actual signal power
    
    % normSignal = cleanSignal / sum(cleanSignal);
    % S2 = mean(normSignal.^2); %// actual signal power
    
    % ADD noise of specified SNR
    snrValue = -5;
    
    % add factory1 noise
    wavin3 = load('.\noise\factory1.mat'); % load babble noise
    wavin3 = wavin3.factory1;
    [noisedSignal,~] = addnoise(cleanSignal, wavin3(1:length(cleanSignal)), snrValue);
    
    % add white noise
    %noisedSignal_white = awgn(cleanSignal, snrValue, 'measured');
    
    awgnNoise = cleanSignal - noisedSignal;
    
    snr_value_1 = GetSNR(cleanSignal, awgnNoise);
    snr_value_2 = snr(cleanSignal, awgnNoise);
        
    %======================%
    noisedSignal = noisedSignal(:,1); % mono to stero
    
    % noiseReductionMethod = 'MultiBandKamath';
    noiseReductionMethod = noiseReductionMethods{iMethod};
    
    IS = 0.25; % seconds of Initial silence
    
    plot_type = 'no';
    %======================%
    %winLenArray = 2.^[8:12];
    winLenArray = 2.^12;
    nWinLen = length(winLenArray);
    enhancedSNR_array = zeros(nWinLen,1);
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
                
            case 'Towsey'
                W = winLen; % samples
                nfft=W;
                wnd=hamming(W);
                
                NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments
                Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
                
                y=segment(noisedSignal,W,SP,wnd);
                Y=fft(y,nfft);
                YPhase=angle(Y(1:floor(end/2)+1,:)); %Noisy Speech Phase
                Y=abs(Y(1:floor(end/2)+1,:)).^Gamma; %Specrogram
                numberOfFrames=size(Y,2);
                FreqResol=size(Y,1);
                
                T1 = length(noisedSignal)/fs;
                %[I1,F1,T1] = wav_to_spec(noisedSignal,fs,winLen,1-SP);
                rT = max(T1)/numberOfFrames;
                
                X = TowseyQUT(Y,rT,IS);
                
                % figure;
                % subplot(211); imagesc(X); axis xy;
                % subplot(212); imagesc(Y); axis xy;
                               
                noiseReducedSignal=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(SP*W));
                
            case 'MultiBandKamath'
                noiseReducedSignal = SSMultibandKamath02(noisedSignal,fs,IS,winLen,SP);
                
            case 'MMSESTSA'
                noiseReducedSignal = MMSESTSA84(noisedSignal,fs,IS,winLen,SP);
                
            case 'Wiener'
                noiseReducedSignal = WienerScalart96(noisedSignal,fs,IS,winLen,SP);
                
            case 'Wavelet-Jie'
                
                % https://ww2.mathworks.cn/help/wavelet/ref/wden.html
                lev = 5;
                wname = 'sym8';                          
                %wname = 'haar';

                %tempSignal = wden(noisedSignal,'rigrsure','h','mln',lev,wname);
                %tempSignal = wden(noisedSignal,'minimaxi','h','mln',lev,wname);
                %tempSignal = wden(noisedSignal,'sqtwolog','h','mln',lev,wname);
                tempSignal = wden(noisedSignal,'modwtsqtwolog','s','mln',lev,wname);
                
                if size(tempSignal,1) < size(tempSignal,2)
                    noiseReducedSignal = tempSignal';
                else
                    noiseReducedSignal = tempSignal;
                end
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
                
                enhancedSNR_array(iWin) = GetSNR(cleanSignal, cleanSignal - noiseReducedSignal);
        end
    end
    %======================%
    % disp(enhancedSNR_array - orgSNR_array);
    diffSNR_array = enhancedSNR_array - orgSNR_array;
    
    improvedSNR{iMethod} = diffSNR_array;
    %======================%
    % figure;
    % plot(diffSNR_array);
    % ylabel('SNR improvement (dB)');
    % title(noiseReductionMethod);
end

improvedSNRmat = cell2mat(improvedSNR);

disp(improvedSNRmat);

if nWinLen == 1
    figure;
    plot(improvedSNRmat, 'LineWidth',2);
    xticklabels(noiseReductionMethods);  
    xtickangle(45);
else    
    figure;
    [M,N] = size(improvedSNRmat);
    for i = 1:N
        plot(improvedSNRmat(:,i), 'LineWidth',2);
        hold on;
    end
    legend(noiseReductionMethods);   
end

%[EOF]
