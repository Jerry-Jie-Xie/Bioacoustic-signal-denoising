function [enhancedSignal] = SSTowsey(initialSignal, fs, IS, param)
%SSTowsey Summary of this function goes here
%   Detailed explanation goes here

W = param.winLen; % samples
nfft=W;
wnd=hamming(W);
NIS=fix((IS*fs-W)/(param.SP*W) +1);%number of initial silence segments
Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
y=segment(initialSignal,W,param.SP,wnd);
Y=fft(y,nfft);
YPhase=angle(Y(1:floor(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:floor(end/2)+1,:)).^Gamma; %Specrogram
numberOfFrames=size(Y,2);
T1 = length(initialSignal)/fs;
secPerFrame = max(T1)/numberOfFrames;

figure; imagesc(Y); axis xy;
title('Initial Spectrogram');
colorbar;

%--apply spectral subtraction
X = TowseyQUT(Y, secPerFrame);

figure; imagesc(X); axis xy;
title('Denoised Spectrogram - TowseyQUT');
colorbar;

enhancedSignal=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(param.SP*W));

end

