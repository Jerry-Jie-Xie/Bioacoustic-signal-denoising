function Y = calculateSpectrogram(initialSignal, param)

W = param.winLen; % samples
nfft=W;
wnd=hamming(W);
Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)
y=segment(initialSignal,W,param.SP,wnd);
Y=fft(y,nfft);
Y=abs(Y(1:floor(end/2)+1,:)).^Gamma; %Specrogram

end