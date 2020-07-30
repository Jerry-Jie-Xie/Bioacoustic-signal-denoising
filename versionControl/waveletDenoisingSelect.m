function enhancedSignal=waveletDenoisingSelect(signal)

% N = length(signal);

%--Decompose speech signal
level = 3;

% fprintf('\n  Input the number of specific wavelet: (1) db13, (2) db40, (3) sym13 or (4) sym21');
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

% computes four filters
[Lo_D, Hi_D, Lo_R, Hi_R] = wfilters(wt);        
[C,L] = wavedec(signal, level, Lo_D, Hi_D);
cA3 = appcoef(C,L,wt,level);         
%extract the levels 3, 2, and 1 detail coefficients from C
% [cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);  
%reconstruct the level 3 approximation from C
A3 = wrcoef('a',C,L,Lo_R,Hi_R,level);
%reconstruct the details at levels 1, 2, and 3, from C
D1 = wrcoef('d',C,L,Lo_R,Hi_R,1);
D2 = wrcoef('d',C,L,Lo_R,Hi_R,2);
D3 = wrcoef('d',C,L,Lo_R,Hi_R,3);

%--thresholding
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

sh = 1;
if sh == 1
    sorh = 's';
elseif sh == 2
    sorh = 'h';
end

% threshold coefficient of details
tD1 = wthresh(D1,sorh,thr_D1);
tD2 = wthresh(D2,sorh,thr_D2);
tD3 = wthresh(D3,sorh,thr_D3);

% compute Inverse DWT 
enhancedSignal = A3 + tD1 + tD2 + tD3;

end
