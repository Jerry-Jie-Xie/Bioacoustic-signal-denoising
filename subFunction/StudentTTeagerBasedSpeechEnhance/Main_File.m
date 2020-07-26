% This is the main MATLAB code file of implementation of algorithm presented in
% "Speech Enhancement Based on Student t Modeling of Teager Energy Operated
% Perceptual Wavelet Packet Coefficients and a Custom Thresholding Function"
% by M. T. Islam, C. Shahnaz, W. P. Zhu and M. O. Ahmad submitted to
% IEEE TSALP.

% Developed by Md Tauhidul Islam, Department of Biomedical Engineering,
% Texas A&M University,College Station, Texas,USA.

% Copyright (c) 2015. Prof Celia Shahnaz.
% All rights reserved.

filename='sp05_babble_sn15.wav';
[noisy, Srate] = audioread(filename);
ms=0.064;                        % Frame Size in Second For Rayleigh
len=fix(ms*Srate);                  % Window length is 25 ms
nfft=len;
NIS  = 5;                        % Number of Silent Frames
SP=.5;                           % Shift percentage is 50% (10ms) Overlap-Add method works good with this value
win=hamming(len);                % Hanning Window Preparation
win2=kaiser(len,1);
y=segment(noisy,len,SP,win);     % Chopping the data into segements for only use in VAD (Voice Activity Detection)
Y=fft(y,nfft);                   % FFT of Chopped data for only use in VAD (Voice Activity Detection)
N=mean(Y(:,1:NIS)')';            % Initial Noise Power Spectrum mean will be used in VAD
[r,c]=size(Y);                   % Size of Y
if rem(len,2)==1                 % Making the length of signal segment even
    len=len+1;
end
len1 = floor(len*SP);            % Progress sample point from frame to frame
len2 = len-len1;
wv='db10';                      % Type of Wavelet Used for Decomposition

%%   Initialisation
Nframes=floor(length(noisy)/len1)-1; % Number of Frames
n_p=ifft(N,nfft);                % Initializing time domain noise signal taken from initial noise frames
Beta=.03;                       % The multiplier for silence regions(if noiseflag=1 in VAD)
k=1;                            % Counter of sample number in frame to frame

% For use in Overlap and Add Method
x_old      = zeros(len1,1);
xfin       = zeros(length(noisy),1); % Final array to save output

%% Main Code
% Estimating noise using IMCRA method
audiowrite('noisy.wav', noisy, 8000);
noise=omlsa('noisy','outomlsa');
xf=noisy;
xf(1:length(noise))=noise(1:length(noise));

for n = 1:Nframes
    % At frame number equal or less than 5
    if n<NIS+1
        insign=win.*noisy(k:k+len-1);  % Windowing
        xi_w=Beta*insign;    % Multiply by beta (A small value) if a silent frame comes
    else
        % Now at frame number more than 5
        insign = win.*noisy(k:k+len-1);
        n_p=win2.*xf(k:k+len-1);
        frames_re=WPT_Percept_TEO(insign,n_p,wv);
        n_p=insign-frames_re;    % Noise Update will be used in next frame
        xi_w=frames_re;               % Enhanced frame will go to overlap add method
    end
    % Overlap Add
    xfin(k:k+len1-1) = x_old(1:len1) + xi_w(1:len1);
    x_old = xi_w(len1+1:len);
    k = k + len1;
    if k>length(noisy)-(len+3*len1)        % Condition to reach at the last frame
        break;
    end
end

enhanced=xfin;  % Save the Enhanced Signal



