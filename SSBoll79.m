function output=SSBoll79(signal,fs,IS,param)

% OUTPUT=SSBOLL79(S,FS,IS)
% Spectral Subtraction based on Boll 79. Amplitude spectral subtraction 
% Includes Magnitude Averaging and Residual noise Reduction
% S is the noisy signal, FS is the sampling frequency and IS is the initial
% silence (noise only) length in seconds (default value is .25 sec)
%
% April-05
% Esfandiar Zavarehei

winLen = param.winLen;
SP =  param.SP;

if (nargin<3 || isstruct(IS))
    IS=.25; %seconds
end
% W=fix(.025*fs); %Window length is 25 ms
W = winLen; % samples
nfft=W;
% SP=.4; %Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
wnd=hamming(W);

% % IGNORE THIS SECTION FOR CAMPATIBALITY WITH ANOTHER PROGRAM FROM HERE.....
% if (nargin>=3 && isstruct(IS))%This option is for compatibility with another programme
%     W=IS.windowsize;
%     SP=IS.shiftsize/W;
%     nfft=IS.nfft;
%     wnd=IS.window;
%     if isfield(IS,'IS')
%         IS=IS.IS;
%     else
%         IS=.25;
%     end
% end
% %--IGNORE THIS SECTION FOR CAMPATIBALITY WITH ANOTHER PROGRAM T0 HERE

NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments
Gamma=1;%Magnitude Power (1 for magnitude spectral subtraction 2 for power spectrum subtraction)

y=segment(signal,W,SP,wnd);
Y=fft(y,nfft);
YPhase=angle(Y(1:fix(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:fix(end/2)+1,:)).^Gamma;%Specrogram
numberOfFrames=size(Y,2);
% FreqResol=size(Y,1);

N=mean(Y(:,1:NIS),2); %initial Noise Power Spectrum mean
NRM=zeros(size(N));% Noise Residual Maximum (Initialization)
NoiseCounter=0;
NoiseLength=9; %This is a smoothing factor for the noise updating

Beta=.03;

YS=Y; %Y Magnitude Averaged
for i=2:(numberOfFrames-1)
    YS(:,i)=(Y(:,i-1)+Y(:,i)+Y(:,i+1))/3;
end

X=zeros(size(Y)); % Initialize X (memory allocation)
for i=1:numberOfFrames
    [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i).^(1/Gamma),N.^(1/Gamma),NoiseCounter); %Magnitude Spectrum Distance VAD
    if SpeechFlag==0
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1); %Update and smooth noise
        NRM=max(NRM,YS(:,i)-N);%Update Maximum Noise Residue
        X(:,i)=Beta*Y(:,i);
    else
        D=YS(:,i)-N; % Specral Subtraction
        if i>1 && i<numberOfFrames %Residual Noise Reduction            
            for j=1:length(D)
                if D(j)<NRM(j)
                    D(j)=min([D(j) YS(j,i-1)-N(j) YS(j,i+1)-N(j)]);
                end
            end
        end
        X(:,i)=max(D,0);
    end
end

output=OverlapAdd2(X.^(1/Gamma),YPhase,W,fix(SP*W));

end