function output=waveletDenoisingFrame(signal, ~, ~, param)

%--wavelet based denosing for each frame
winLen = param.winLen;
SP =  param.SP;

W = winLen; % samples
wnd=hamming(W);

pre_emph=0;
signal=filter([1 -pre_emph],1,signal);

y=segment(signal,W,SP,wnd); % This function chops the signal into frames

numberOfFrames = size(y, 2);

h=waitbar(0,'Wait...');

lev = 8; % level of wavelet transform
wname = 'dmey'; % type of wavelets

X=zeros(size(y)); % Initialize X (memory allocation)
for i=1:numberOfFrames

    frameY = wden(y(:,i),'modwtsqtwolog','s','mln',lev,wname);    
    X(:,i) = frameY; %Obtain the new Cleaned value
    
    waitbar(i/numberOfFrames,h,num2str(fix(100*i/numberOfFrames)));
end

close(h);

output=OverlapAdd1(X,W,fix(SP*W)); %Overlap-add Synthesis of speech
output=filter(1,[1 -pre_emph],output); %Undo the effect of Pre-emphasis

end
