function y=semisoft_pdf(x,t,n,s)

% This function performs thresholding operation on an input signal
% x---Signal
% t---First Threshold
% n=Noise Power
% s=Signal Power
% Noise Power for Further Improvement
t2=(2)*t; % Take lambda2=lambda*2;
[pr,pn]=Probability_Calc(s,n); % Calculation of probabilty of speech and noise
for k=1:length(pr)
    alph(k)=(1+pn(k))./(2*(1+pr(k))); % Shape Parameter alpha comes from Speech and noise probablity
    beta(k)=(4*(1+pn(k))./(1+pr(k)));% Shape Parameter beta comes from Speech and noise probablity
end
for kk=1:length(alph)
    if isnan(alph(kk))
        alph(kk)=0.1;
    end
    if isnan(beta(kk))
        beta(kk)=2;
    end
end

for i=1:length(x)   
    alpha=alph(i);
    bet=beta(i);
    if abs(x(i))<t        % First Region where Y<lambda      
        y(i)=(1-alpha)*sign(x(i))*((abs(x(i)))^bet)/(t^(bet-1));        
    elseif abs(x(i))>t2    % Second Region where Y>lambda        
        y(i)=x(i);        
    else        
        y1(i)=x(i); % One threshold Modified Hard Thresholding
        y2(i)=sign(x(i)) .* t2/(t2-t) .* (abs(x(i))-t); % Semisoft threshold Function
        y(i)=alpha*y2(i)+(1-alpha)*y1(i); % Linear Combination of MH and SS        
    end    
end
