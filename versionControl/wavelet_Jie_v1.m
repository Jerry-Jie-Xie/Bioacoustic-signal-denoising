clc; close all; clear;
%Compile with MATLAB Wavelet Toolbox 3.0
%
%De-Noising Speech Signal Using Discrete Wavelet Transform
%
%The de-noising procedure proceeds in three steps:
%[1]
%Decomposition. Choose a wavelet, and choose a level N. Compute the wavelet
%decomposition of the signal s at level N.
%[2]
%Detail coefficients thresholding. For each level from 1 to N, select
%a threshold and apply soft thresholding to the detail coefficients.
%[3]
%Reconstruction. Compute wavelet reconstruction based on the original 
%approximation coefficients of level N and the modified detail coefficients 
%of levels from 1 to N.
%
%Programmed by Student of Lambung Mangkurat University
%Nadya Amalia

fprintf('--- De-noising Speech Signal Using Discrete Wavelet Transform ---\n\n');

%load example of speech sound
fprintf('-> Step 1/6: Load "female.wav" - ');
audioPath = 'C:\Jie\work\data\noiseReduction\Rufous Whistler\07.wav';
[truesignal,Fs] = audioread(audioPath);    
truesignal = truesignal(:,1); % mono to stero
% amp = 20;
% truesignal = amp*truesignal;               
N = length(truesignal);
fprintf('OK\n');

%----------------------------%
%  add white Gaussian noise  %
%----------------------------%
fprintf('-> Step 2/6: Add white Gaussian noise - ');
%the scalar SNR specifies the signal-to-noise ratio per sample, in dB
sn = -5;
%add white Gaussian noise to a signal
truesignalN = awgn(truesignal,sn,'measured');              
fprintf('OK\n');

%---------------------------%
%           DWT :           %
%   wavelet decomposition   %
%---------------------------%
fprintf('-> Step 3/6: Decompose speech signal - ');
level = 3;
% fprintf('\n   Input the number of specific wavelet: (1) db13, (2) db40, (3) sym13 or (4) sym21');
% wname = input('\n   wname = ');
wname = 1;
if wname == 1
    wt = 'db13';
elseif wname == 2
    wt = 'db40';
elseif wname == 3
    wt = 'sym13';
elseif wname == 4
    wt = 'sym21';
end
%computes four filters
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wt);        
[C,L] = wavedec(truesignalN,level,Lo_D,Hi_D);
cA3 = appcoef(C,L,wt,level);         
%extract the levels 3, 2, and 1 detail coefficients from C
[cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);  
%reconstruct the level 3 approximation from C
A3 = wrcoef('a',C,L,Lo_R,Hi_R,level);
%reconstruct the details at levels 1, 2, and 3, from C
D1 = wrcoef('d',C,L,Lo_R,Hi_R,1);
D2 = wrcoef('d',C,L,Lo_R,Hi_R,2);
D3 = wrcoef('d',C,L,Lo_R,Hi_R,3);
%a = approximation
%d = detail
fprintf('OK\n');

%---------------------------%
%        thresholding       %
%---------------------------%
fprintf('-> Step 4/6: Thresholding - ');
%TPTR = 'rigrsure', adaptive threshold selection using principle of Stein's
%Unbiased Risk Estimate
%TPTR = 'heursure', heuristic variant of the first option
%TPTR = 'sqtwolog', threshold is sqrt(2*log(length(X)))
%TPTR = 'minimaxi', minimax thresholding
% fprintf('\n   Input the number of threshold selection rule : (1) heursure, (2) rigrsure, (3) minimaxi or (4) sqtwolog');
% tr = input('\n   threshold selection rule = ');
tr = 4;
if tr == 1
    tptr = 'heursure';
elseif tr == 2
    tptr = 'rigrsure';
elseif tr == 3
    tptr = 'minimaxi';
elseif tr == 4
    tptr = 'sqtwolog'; 
end
thr_D1 = thselect(D1,tptr);
thr_D2 = thselect(D2,tptr);
thr_D3 = thselect(D3,tptr);

%Hard thresholding is the simplest method but soft thresholding has nice
%mathematical properties. Hard threshold signal is x if x>thr, and is 0 if
%x<=thr. And the soft threshold signal is sign(x)(x-thr) if x>thr and is 0
%if x<=thr
% fprintf('\n   Input the number of threshold type: (1) soft or (2) hard');
% sh = input('\n   threshold = ');
sh = 1;
if sh == 1
    sorh = 's';
elseif sh == 2
    sorh = 'h';
end
%threshold coefficient of details
tD1 = wthresh(D1,sorh,thr_D1);
tD2 = wthresh(D2,sorh,thr_D2);
tD3 = wthresh(D3,sorh,thr_D3);
fprintf('OK\n');

%--------------------------%
%   compute Inverse DWT    %
%--------------------------%
fprintf('-> Step 5/6: Compute Inverse DWT - ');
denoised = A3 + tD1 + tD2 + tD3;
err = max(abs(truesignalN-denoised));
fprintf('OK\n');

%---------------------------%
%       compute SNR         %
%---------------------------%
fprintf('-> Step 6/6: Compute SNR - ');
%SNR - Signal to Noise Ratio
% SNR = snr(truesignal,truesignalN);
NoisySNR = 20*log10(norm(truesignal(:)) / norm (truesignal(:)-truesignalN(:)));
% SNR = snr(truesignal,denoised);
DenoisedSNR = 20*log10(norm(truesignal(:)) / norm (truesignal(:)-denoised(:)));
fprintf('OK\n');


%-----------------    Listen Result   ------------------------------------      

% fprintf('\n   Play the Original Sound:');
% wavplay(truesignal,Fs,'sync')
% fprintf(' OK');
% fprintf('\n   Play the Noisy Sound:');
% wavplay(truesignalN,Fs,'sync')
% fprintf(' OK');
% fprintf('\n   Play the Denoised Sound:');
% wavplay(denoised,Fs,'sync')
% fprintf(' OK\n');


%-----------------    Display Figure   ------------------------------------      
figure(1)
subplot(3,1,1); plot(truesignal); title('True Speech Signal'); 
subplot(3,2,3); plot(A3); title('Approximation A3')
subplot(3,2,4); plot(D3); title('Detail D3')
subplot(3,2,5); plot(D2); title('Detail D2')
subplot(3,2,6); plot(D1); title('Detail D1')
figure(2)
subplot(3,1,1); plot(truesignal); title('True Speech Signal'); 
subplot(3,2,3); plot(A3); title('Approximation A3')
subplot(3,2,4); plot(tD3); title('Denoised Detail D3')
subplot(3,2,5); plot(tD2); title('Denoised Detail D2')
subplot(3,2,6); plot(tD1); title('Denoised Detail D1')
%display the comparison of truesignal, noisy signal, and denoised signal
figure(3)
subplot(3,1,1); plot(truesignal); title('True Speech Signal');
xlabel('Samples'); ylabel('Amplitude');
subplot(3,1,2); plot(truesignalN); title('Noisy Speech Signal');
xlabel('Samples'); ylabel('Amplitude');
subplot(3,1,3); plot(denoised); title('De-noised Speech Signal');
xlabel('Samples'); ylabel('Amplitude');
figure(4)
subplot(1,3,1); specgram(truesignal,512,Fs); title('True Speech Signal');
subplot(1,3,2); specgram(truesignalN,512,Fs); title('Noisy Speech Signal');
subplot(1,3,3); specgram(denoised,512,Fs); title('De-noised Speech Signal');


