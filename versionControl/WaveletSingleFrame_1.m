function output=WaveletSingleFrame(signal,fs,IS,winLen,SP)

%--wavelet based denosing for each frame

if (nargin<3 || isstruct(IS))
    IS=.25; %Initial Silence or Noise Only part in seconds
end
W = winLen; % samples
nfft=W;
% SP=.4; %Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
wnd=hamming(W);

pre_emph=0;
signal=filter([1 -pre_emph],1,signal);

NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments

y=segment(signal,W,SP,wnd); % This function chops the signal into frames
Y=fft(y,nfft); %????????
YPhase=angle(Y(1:fix(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:fix(end/2)+1,:));%Specrogram
numberOfFrames=size(Y,2);
FreqResol=size(Y,1);

% N=mean(Y(:,1:NIS)')'; %initial Noise Power Spectrum mean
N=mean(Y(:,1:NIS),2);
LambdaD=mean((Y(:,1:NIS)').^2)';%initial Noise Power Spectrum variance
alpha=.99; %used in smoothing xi (For Deciesion Directed method for estimation of A Priori SNR)

NoiseCounter=0;
NoiseLength=9;%This is a smoothing factor for the noise updating

G=ones(size(N));%Initial Gain used in calculation of the new xi
Gamma=G;

Gamma1p5=gamma(1.5); %Gamma function at 1.5
X=zeros(size(Y)); % Initialize X (memory allocation)

h=waitbar(0,'Wait...');

lev = 5; % level of wavelet transform
wname = 'sym8'; % type of wavelets

for i=1:numberOfFrames
     
    %--VAD and Noise Estimation START
    if i<=NIS % If initial silence ignore VAD
        SpeechFlag=0;
        NoiseCounter=100;
    else % Else Do VAD
        [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i),N,NoiseCounter); %Magnitude Spectrum Distance VAD
    end
    
    if SpeechFlag==0 % If not Speech Update Noise Parameters
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1); %Update and smooth noise mean
        LambdaD=(NoiseLength*LambdaD+(Y(:,i).^2))./(1+NoiseLength); %Update and smooth noise variance
    end
    %--VAD and Noise Estimation END
    
    gammaNew=(Y(:,i).^2)./LambdaD; %A postiriori SNR
    xi=alpha*(G.^2).*Gamma+(1-alpha).*max(gammaNew-1,0); %Decision Directed Method for A Priori SNR
    
    Gamma=gammaNew;
    nu=Gamma.*xi./(1+xi); % A Function used in Calculation of Gain
    G=(Gamma1p5*sqrt(nu)./Gamma).*exp(-nu/2)...
        .*((1+nu).*besseli(0,nu/2) +nu.*besseli(1,nu/2)); %MMSE Gain
    
    Indx=find(isnan(G) | isinf(G)); %Sometimes (occasionally) Gain values become 
                                    %infinite which are replaced by Wiener
                                    %Gain in that case, that is:
    G(Indx)=xi(Indx)./(1+xi(Indx)); %Wiener Fileter for out of range values  
    
    frameY = wden(Y(:,i),'modwtsqtwolog','s','mln',lev,wname);

    X(:,i)=frameY; %Obtain the new Cleaned value
    
    waitbar(i/numberOfFrames,h,num2str(fix(100*i/numberOfFrames)));
end

close(h);
output=OverlapAdd2(X,YPhase,W,fix(SP*W)); %Overlap-add Synthesis of speech
output=filter(1,[1 -pre_emph],output); %Undo the effect of Pre-emphasis

end