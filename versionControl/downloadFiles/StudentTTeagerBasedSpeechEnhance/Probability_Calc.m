function [PH1,q]=Probability_Calc(y,n)

% This function calculates the speech and noise presence probability 
% based on speech and noise power
y=[y' zeros(1,length(y))];
n=[n' zeros(1,length(n))];

% ***************************************************************@
% Inputs:
% %    y,  input signal power
% %    n,  noise power
% % Output:
% %    PH1, speech signal probability
% %    out, noise probabilty
% %
% % ***************************************************************@

% % 1) Parameters of Short Time Fourier Analysis:

M_ref=length(y);		% 1.2) Size of analysis window
Mo_ref=M_ref;	% 1.3) Number of overlapping samples in consecutive frames
Fs=8000;
% % 2) Parameters of Noise Spectrum Estimate
% % 3) Parameters of a Priori Probability for Signal-Absence Estimate
alpha_xi_ref=0.7;	% 3.1) Recursive averaging parameter

w_xi_local=1; 	% 3.2) Size of frequency local smoothing window function
w_xi_global=15; 	% 3.3) Size of frequency local smoothing window function
f_u=10e3; 		% 3.4) Upper frequency threshold for global decision
f_l=50; 		% 3.5) Lower frequency threshold for global decision
P_min=0.005; 		% 3.6) Lower bound constraint
xi_lu_dB=-5; 	% 3.7) Upper threshold for local decision
xi_ll_dB=-10; 	% 3.8) Lower threshold for local decision
xi_gu_dB=-5; 	% 3.9) Upper threshold for global decision
xi_gl_dB=-10; 	% 3.10) Lower threshold for global decision
xi_fu_dB=-5; 	% 3.11) Upper threshold for local decision
xi_fl_dB=-10; 	% 3.12) Lower threshold for local decision
xi_mu_dB=10; 	% 3.13) Upper threshold for xi_m
xi_ml_dB=0; 		% 3.14) Lower threshold for xi_m
q_max=0.998; 		% 3.15) Upper limit constraint
%
% % 4) Parameters of "Decision-Directed" a Priori SNR Estimate
alpha_eta_ref=0.95;	% 4.1) Recursive averaging parameter
eta_min_dB=-18;	% 4.2) Lower limit constraint
%
eta_min=10^(eta_min_dB/10);
M=M_ref;
Mo=Mo_ref;
alpha_eta=alpha_eta_ref;
alpha_xi=alpha_xi_ref;

% window function
win=hamming(M);
win2=win.^2;
Mno=M-Mo;
W0=win2(1:Mno);
for k=Mno:Mno:M-1
    swin2=lnshift(win2,k);
    W0=W0+swin2(1:Mno);
end
b_xi_local=hanning(2*w_xi_local+1);
b_xi_local=b_xi_local/sum(b_xi_local);  % normalize the window function
b_xi_global=hanning(2*w_xi_global+1);
b_xi_global=b_xi_global/sum(b_xi_global);   % normalize the window function
M21=M;
k_u=round(f_u/Fs*M+1);  % Upper frequency bin for global decision
k_l=round(f_l/Fs*M+1);  % Lower frequency bin for global decision
k_u=min(k_u,M21);
k2_local=round(500/Fs*M+1);
k3_local=round(3500/Fs*M+1);
eta_2term=1; xi=0; xi_frame=0;
Ya2=teager(y);
lambda_d=teager(n);
gamma=Ya2./max(lambda_d,1e-10);
eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
eta=max(eta,eta_min);
% 4. A Priori Probability for Signal-Absence Estimate
xi=alpha_xi*xi+(1-alpha_xi)*eta;
xi_local=conv(xi,b_xi_local);
xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
xi_global=conv(xi,b_xi_global);
xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
dxi_frame=xi_frame;
xi_frame=mean(xi(k_l:k_u));
dxi_frame=xi_frame-dxi_frame;
if xi_local>0, xi_local_dB=10*log10(xi_local); else xi_local_dB=-100; end
if xi_global>0, xi_global_dB=10*log10(xi_global); else xi_global_dB=-100; end
if xi_frame>0, xi_frame_dB=10*log10(xi_frame); else xi_frame_dB=-100; end

P_local=ones(M21,1);
P_local(xi_local_dB<=xi_ll_dB)=P_min;
idx=find(xi_local_dB>xi_ll_dB & xi_local_dB<xi_lu_dB);
P_local(idx)=P_min+(xi_local_dB(idx)-xi_ll_dB)/(xi_lu_dB-xi_ll_dB)*(1-P_min);

P_global=ones(M21,1);
P_global(xi_global_dB<=xi_gl_dB)=P_min;
idx=find(xi_global_dB>xi_gl_dB & xi_global_dB<xi_gu_dB);
P_global(idx)=P_min+(xi_global_dB(idx)-xi_gl_dB)/(xi_gu_dB-xi_gl_dB)*(1-P_min);

m_P_local=mean(P_local(3:(k2_local+k3_local-3)));    % average probability of speech presence
if m_P_local<0.25
    P_local(k2_local:k3_local)=P_min;    % reset P_local (frequency>500Hz) for low probability of speech presence
end

tone_flag=0;
if tone_flag               
    if (m_P_local<0.5) && (l>120)
        idx=find( lambda_dav_long(8:(M21-8)) > 2.5*(lambda_dav_long(10:(M21-6))+lambda_dav_long(6:(M21-10))) );
        P_local([idx+6;idx+7;idx+8])=P_min;   % remove interfering tonals
    end
end            

if xi_frame_dB<=xi_fl_dB
    P_frame=P_min;
elseif dxi_frame>=0
    xi_m_dB=min(max(xi_frame_dB,xi_ml_dB),xi_mu_dB);
    P_frame=1;
elseif xi_frame_dB>=xi_m_dB+xi_fu_dB
    P_frame=1;
elseif xi_frame_dB<=xi_m_dB+xi_fl_dB
    P_frame=P_min;
else
    P_frame=P_min+(xi_frame_dB-xi_m_dB-xi_fl_dB)/(xi_fu_dB-xi_fl_dB)*(1-P_min);
end
broad_flag =1;
if broad_flag
    q=1-P_global.*P_local*P_frame;
else
    q=1-P_local*P_frame;
end
q=min(q,q_max);

gamma=Ya2./max(lambda_d,1e-10);
eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
eta=max(eta,eta_min);
v=gamma.*eta./(1+eta);
PH1=zeros(M21,1);
idx=find(q<0.9);
q=q';
PH1(idx)=1./(1+q(idx)./(1-q(idx)).*(1+eta(idx)).*exp(-v(idx)));


