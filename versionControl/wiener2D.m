function [enhancedSignal] = wiener2D(initialSignal,fs,IS,param)
% wiener2D Summary of this function goes here
%   Detailed explanation goes here
%--1D signal to 2D
% W = param.winLen; % samples
W = 512; % samples

nfft=W;
wnd=hamming(W);
NIS=fix((IS*fs-W)/(param.SP*W) +1); % number of initial silence segments
Gamma=1; % Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
y=segment(initialSignal,W,param.SP,wnd);
Y=fft(y,nfft);
YPhase=angle(Y(1:floor(end/2)+1,:)); % Noisy Speech Phase
Y=abs(Y(1:floor(end/2)+1,:)).^Gamma; % Specrogram

figure; imagesc(Y); axis xy;
title('Initial Spectrogram');
colorbar;

%--apply Wiener filter
[X, noise] = wiener2(Y, [5 5]);

figure; imagesc(X); axis xy;
title('Denoised Spectrogram - Wiener');
colorbar;

enhancedSignal=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(param.SP*W));

end

