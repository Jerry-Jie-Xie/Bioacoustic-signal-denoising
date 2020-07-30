function t=thresh(n,s,m)
% This function calculates the threshold
% n=Noise power
% s= signal Power
% m= Mean of the signal
v=2; %Degree of freedom is taken as 2
n=mean(n);
s=mean(s);
if s>n
r=abs(s-n); % Signal Power=Noisy Signal Power-Noise Power
else
r=0.0001*n; % For noise power higher than signal power difference is taken 
% a small fraction of noise power
end
k=(v+1)/(2*v);
gm=r/n;
t1=sqrt(n)*(n+r);
t2=sqrt(1+gm)+2+gm;
t3=t1/(t2*k);
t=sqrt(t3)*sqrt(1/gm);
t=abs((t));