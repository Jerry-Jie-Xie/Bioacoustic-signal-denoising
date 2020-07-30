function [n1,s1]=npow_teo(sn,n)

sn_t=teager(sn); % TE Operation on Noisy Signal PWPT Coefficients 
n_t=teager(n);   % TE Operation on Noise PWPT Coefficients 

% n1=mean(abs(n_t)); % Noise Power
% s1=mean(abs(sn_t)); % Noisy Signal Power
n1=(abs(n_t)); % Noise Power
s1=(abs(sn_t)); % Noisy Signal Power
% n1=var(n);
% s1=var(sn);
% 
% if s1>n1
% r1=abs(s1-n1); % Signal Power=Noisy Signal Power-Noise Power
% else
% %     r1=n1;% For Rayleigh
% r1=0.0001*n1; % For index cal in student
end


