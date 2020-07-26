function output=SSMultibandKamath02(signal,fs,IS,winLen,SP)

% OUTPUT=SSMULTIBANDKAMATH02(S,FS,IS)
% Multi-band Spectral subtraction [Kamath2002]
% subtraction with adjusting subtraction factor. the adjustment is
% according to local a postriori SNR and the frequency band.
% S is the noisy signal, FS is the sampling frequency and IS is the initial
% silence (noise only) length in seconds (default value is .25 sec)
%
% April-04
% Esfandiar Zavarehei

if (nargin<3 | isstruct(IS))
    IS=.025; %seconds
end
% W=fix(.025*fs); %Window length is 25 ms
W = winLen; % samples
nfft=W;
% SP=.4; %Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
wnd=hamming(W);

% IGNORE THIS SECTION FOR CAMPATIBALITY WITH ANOTHER PROGRAM FROM HERE.....
if (nargin>=3 & isstruct(IS))%This option is for compatibility with another programme
    W=IS.windowsize;
    SP=IS.shiftsize/W;
    nfft=IS.nfft;
    wnd=IS.window;
    if isfield(IS,'IS')
        IS=IS.IS;
    else
        IS=.25;
    end
end
% .......IGNORE THIS SECTION FOR CAMPATIBALITY WITH ANOTHER PROGRAM T0 HERE

NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments
Gamma=2;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)

y=segment(signal,W,SP,wnd);
Y=fft(y,nfft);
YPhase=angle(Y(1:fix(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:fix(end/2)+1,:)).^Gamma;%Specrogram
numberOfFrames=size(Y,2);
FreqResol=size(Y,1);

% N=mean(Y(:,1:NIS)')'; %initial Noise Power Spectrum mean
N=mean(Y(:,1:NIS),2); %initial Noise Power Spectrum mean

NoiseCounter=0;
NoiseLength=9;%This is a smoothing factor for the noise updating

Beta=.03;
minalpha=1;
maxalpha=5;
minSNR=-5;
maxSNR=20;
alphaSlope=(minalpha-maxalpha)/(maxSNR-minSNR);
alphaShift=maxalpha-alphaSlope*minSNR;

BN=Beta*N;

%Delta is a frequency dependent coefficient
Delta=1.5*ones(size(BN));
Delta(1:fix((-2000+fs/2)*FreqResol*2/fs))=2.5; %if the frequency is lower than FS/2 - 2KHz
Delta(1:fix(1000*FreqResol*2/fs))=1; %if the frequency is lower than 1KHz


for i=1:numberOfFrames   
    [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i).^(1/Gamma),N.^(1/Gamma),NoiseCounter); %Magnitude Spectrum Distance VAD
    if SpeechFlag==0
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1); %Update and smooth noise
        BN=Beta*N;
    end
    
    SNR=10*log(Y(:,i)./N);
    alpha=alphaSlope*SNR+alphaShift;
    alpha=max(min(alpha,maxalpha),minalpha);
    
    D=Y(:,i)-(Delta.*alpha.*N); %Nonlinear (Non-uniform) Power Specrum Subtraction
    
    X(:,i)=max(D,BN); %if BY>D X=BY else X=D which sets very small values of subtraction result to an attenuated 
                      %version of the input power spectrum.
end

output=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(SP*W));