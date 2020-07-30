function output=WaveletSingleFrame_2(signal,fs,IS,len,SP)

pre_emph=0;
signal=filter([1 -pre_emph],1,signal);

%%
nfft=len;

% SP=.5;                           % Shift percentage is 50% (10ms) Overlap-Add method works good with this value
% NIS = 5;
NIS=fix((IS*fs-len)/(SP*len) +1); %number of initial silence segments

win=hamming(len);                % Hamming Window Preparation
win2=kaiser(len,1);
y=segment(signal,len,SP,win);     % Chopping the data into segements for only use in VAD (Voice Activity Detection)
Y=fft(y,nfft);                   % FFT of Chopped data for only use in VAD (Voice Activity Detection)
N=mean(Y(:,1:NIS)')';            % Initial Noise Power Spectrum mean will be used in VAD
[r,c]=size(Y);                   % Size of Y
if rem(len,2)==1                 % Making the length of signal segment even
    len=len+1;
end
len1 = floor(len*SP);            % Progress sample point from frame to frame
len2 = len-len1;
wv='db10';                      % Type of Wavelet Used for Decomposition

% wv = 'sym8';
%%  Initialisation
Nframes=floor(length(signal)/len1)-1; % Number of Frames
n_p=ifft(N,nfft);                % Initializing time domain noise signal taken from initial noise frames
Beta=.03;                       % The multiplier for silence regions(if noiseflag=1 in VAD)
k=1;                            % Counter of sample number in frame to frame

% For use in Overlap and Add Method
x_old      = zeros(len1,1);
xfin       = zeros(length(signal),1); % Final array to save output

%% Main Code
% Estimating noise using IMCRA method
% audiowrite('noisy.wav', noisy, 8000, 16);
audiowrite('noisy.wav', signal, 8000, 'BitsPerSample', 16);
noise=omlsa('noisy.wav','outomlsa');
xf=signal;
xf(1:length(noise))=noise(1:length(noise));

for n = 1:Nframes
    % At frame number equal or less than 5
    if n<NIS+1
        insign=win.*signal(k:k+len-1);  % Windowing
        xi_w=Beta*insign;    % Multiply by beta (A small value) if a silent frame comes
    else
        % Now at frame number more than 5
        insign = win.*signal(k:k+len-1);
        n_p=win2.*xf(k:k+len-1);
        
        %frames_re=WPT_Percept_TEO(insign,n_p,wv);
       
        lev = 13;
        wname = 'db10';
        frames_re = wden(insign,'rigrsure','s','mln',lev,wname);
        
        n_p=insign-frames_re;    % Noise Update will be used in next frame
        xi_w=frames_re;               % Enhanced frame will go to overlap add method
    end
    % Overlap Add
    xfin(k:k+len1-1) = x_old(1:len1) + xi_w(1:len1);
    x_old = xi_w(len1+1:len);
    k = k + len1;
    if k>length(signal)-(len+3*len1)        % Condition to reach at the last frame
        break;
    end
end

output=xfin;  % Save the Enhanced Signal
output=filter(1,[1 -pre_emph],output); %Undo the effect of Pre-emphasis
