
%****************************************************
% Function to Perform Perceptual Wavelet Packet Transform on a signal
% Ref:Perceptually motivated wavelet packet transform for bioacoustic signal enhancement
% Yao Ren,Michael T. Johnson, and Jidong Tao
% J. Acoust. Soc. July 2008 Acoustical Society of America
%***************************************************
function frames_re=WPT_Percept_TEO(insign,noise,wv)


%% PWPT Decomposition of Noisy Speech Frame

[coef,len]=wavedec(insign,6,wv);
lenm=len;
frames_re=insign;
coef_re=coef;
s_no=cumsum(len);

[coef1,len1]=wavedec(coef(s_no(2):s_no(3)),1,wv);
coef1_re=coef1;
len1m=len1;
len1=cumsum(len1);

[coef2,len2]=wavedec(coef((s_no(3)+1):s_no(4)),2,wv);
len2m=len2;
coef2_re=coef2;
len2=cumsum(len2);

[coef3,len3]=wavedec(coef2((len2(2)+1):len2(3)),1,wv);
len3m=len3;
coef3_re=coef3;
len3=cumsum(len3);


% all nodes at level six are over

[coef4,len4]=wavedec(coef(s_no(4):s_no(5)),2,wv);
len4m=len4;
coef4_re=coef4;
len4=cumsum(len4);

[coef5,len5]=wavedec(coef4((len4(2)+1):len4(3)),1,wv);
len5m=len5;
coef5_re=coef5;
len5=cumsum(len5);

[coef7,len7]=wavedec(coef(s_no(5):s_no(6)),3,wv);
len7m=len7;
coef7_re=coef7;
len7=cumsum(len7);

[coef8,len8]=wavedec(coef7((len7(2)+1):len7(3)),1,wv);
len8m=len8;
coef8_re=coef8;
len8=cumsum(len8);

[coef9,len9]=wavedec(coef7((len7(3)+1):len7(4)),2,wv);
len9m=len9;
coef9_re=coef9;
len9=cumsum(len9);

[coef10,len10]=wavedec(coef((s_no(6)+1):s_no(7)),3,wv);
len10m=len10;
coef10_re=coef10;
len10=cumsum(len10);

[coef11,len11]=wavedec(coef10((len(3)+1):len10(4)),1,wv);
len11m=len11;
len11=cumsum(len11);

%% PWPT Decomposition of Noise Frame

[coefn,lenn]=wavedec(noise,6,wv);
lenmn=lenn;
frames_ren=noise;
coef_ren=coefn;
s_non=cumsum(lenn);

[coef1n,len1n]=wavedec(coefn(s_non(2):s_non(3)),1,wv);
coef1_ren=coef1n;
len1mn=len1n;
len1n=cumsum(len1n);

[coef2n,len2n]=wavedec(coefn((s_non(3)+1):s_non(4)),2,wv);
len2mn=len2n;
coef2_ren=coef2n;
len2n=cumsum(len2n);

[coef3n,len3n]=wavedec(coef2n((len2n(2)+1):len2n(3)),1,wv);
len3mn=len3n;
coef3_ren=coef3n;
len3n=cumsum(len3n);


% all nodes at level six are over
[coef4n,len4n]=wavedec(coefn(s_non(4):s_non(5)),2,wv);
len4mn=len4n;
coef4_ren=coef4n;
len4n=cumsum(len4n);

[coef5n,len5n]=wavedec(coef4n((len4n(2)+1):len4n(3)),1,wv);
len5mn=len5n;
coef5_ren=coef5n;
len5n=cumsum(len5n);

[coef7n,len7n]=wavedec(coefn(s_non(5):s_non(6)),3,wv);
len7mn=len7n;
coef7_ren=coef7n;
len7n=cumsum(len7n);

[coef8n,len8n]=wavedec(coef7n((len7n(2)+1):len7n(3)),1,wv);
len8mn=len8n;
coef8_ren=coef8n;
len8n=cumsum(len8n);

[coef9n,len9n]=wavedec(coef7n((len7n(3)+1):len7n(4)),2,wv);
len9mn=len9n;
coef9_ren=coef9n;
len9n=cumsum(len9n);

[coef10n,len10n]=wavedec(coefn((s_non(6)+1):s_non(7)),3,wv);
len10mn=len10n;
coef10_ren=coef10n;
len10n=cumsum(len10n);

[coef11n,len11n]=wavedec(coef10n((lenn(3)+1):len10n(4)),1,wv);
len11mn=len11n;
len11n=cumsum(len11n);


%% Seperating The PWPT Coefficients of Noisy Speech at 24 Subbands

energy1=(coef(1:s_no(1)));
energy2=(coef((s_no(1)+1):s_no(2)));
energy3=(coef1(1:len1(1)) );
energy4=(coef1((len1(1)+1):len1(2)));
energy5=(coef2(1:len2(1)));
energy6=(coef2( (len2(1)+1):len2(2)));
energy7=(coef3(1:len3(1)));
energy8=(coef3((len3(1)+1):len3(2)));
energy9=(coef4(1:len4(1)));
energy10=(coef4((len4(1)+1):len4(2)));
energy11=(coef5(1:len5(1)));
energy12=(coef5((len5(1)+1):len5(2)));
energy13=(coef7(1:len7(1)));
energy14=(coef7((len7(1)+1):len7(2)));
energy15=(coef8(1:len8(1)));
energy16=(coef8((len8(1)+1):len8(2)));
energy17=(coef9(1:len9(1)));
energy18=(coef9((len9(1)+1):len9(2)));
energy19=(coef9((len9(2)+1):len9(3)));
energy20=(coef10(1:len10(1)));
energy21=(coef10((len10(1)+1):len10(2)));
energy22=(coef10((len10(2)+1):len10(3)));
energy23=(coef11(1:len11(1)));
energy24=(coef11((len11(1)+1):len11(2)));

%% Seperating The PWPT Coefficients of Noise at 24 Subbands

energy1x=(coefn(1:s_no(1)));
energy2x=(coefn((s_no(1)+1):s_no(2)));
energy3x=(coef1n(1:len1(1)) );
energy4x=(coef1n((len1(1)+1):len1(2)));
energy5x=(coef2n(1:len2(1)));
energy6x=(coef2n((len2(1)+1):len2(2) ) );
energy7x=(coef3n(1:len3(1)));
energy8x=(coef3n((len3(1)+1):len3(2)));
energy9x=(coef4n(1:len4(1)));
energy10x=(coef4n((len4(1)+1):len4(2)));
energy11x=(coef5n(1:len5(1)));
energy12x=(coef5n((len5(1)+1):len5(2)));
energy13x=(coef7n(1:len7(1)));
energy14x=(coef7n((len7(1)+1):len7(2)));
energy15x=(coef8n(1:len8(1)));
energy16x=(coef8n((len8(1)+1):len8(2)));
energy17x=(coef9n(1:len9(1)));
energy18x=(coef9n((len9(1)+1):len9(2)));
energy19x=(coef9n((len9(2)+1):len9(3)));
energy20x=(coef10n(1:len10(1)));
energy21x=(coef10n((len10(1)+1):len10(2)));
energy22x=(coef10n((len10(2)+1):len10(3)));
energy23x=(coef11n(1:len11(1)));
energy24x=(coef11n((len11(1)+1):len11(2)));


%% Determining the Noise and Signal Power Using Teager Energy Operator
[n1,r1]=npow_teo(energy1,energy1x);m1=mean(energy1);
[n2,r2]=npow_teo(energy2,energy2x);m2=mean(energy2);
[n3,r3]=npow_teo(energy3,energy3x);m3=mean(energy3);
[n4,r4]=npow_teo(energy4,energy4x);m4=mean(energy4);
[n5,r5]=npow_teo(energy5,energy5x);m5=mean(energy5);
[n6,r6]=npow_teo(energy6,energy6x);m6=mean(energy6);
[n7,r7]=npow_teo(energy7,energy7x);m7=mean(energy7);
[n8,r8]=npow_teo(energy8,energy8x);m8=mean(energy8);
[n9,r9]=npow_teo(energy9,energy9x);m9=mean(energy9);
[n10,r10]=npow_teo(energy10,energy10x);m10=mean(energy10);
[n11,r11]=npow_teo(energy11,energy11x);m11=mean(energy11);
[n12,r12]=npow_teo(energy12,energy12x);m12=mean(energy12);
[n13,r13]=npow_teo(energy13,energy13x);m13=mean(energy13);
[n14,r14]=npow_teo(energy14,energy14x);m14=mean(energy14);
[n15,r15]=npow_teo(energy15,energy15x);m15=mean(energy15);
[n16,r16]=npow_teo(energy16,energy16x);m16=mean(energy16);
[n17,r17]=npow_teo(energy17,energy17x);m17=mean(energy17);
[n18,r18]=npow_teo(energy18,energy18x);m18=mean(energy18);
[n19,r19]=npow_teo(energy19,energy19x);m19=mean(energy19);
[n20,r20]=npow_teo(energy20,energy20x);m20=mean(energy20);
[n21,r21]=npow_teo(energy21,energy21x);m21=mean(energy21);
[n22,r22]=npow_teo(energy22,energy22x);m22=mean(energy22);
[n23,r23]=npow_teo(energy23,energy23x);m23=mean(energy23);
[n24,r24]=npow_teo(energy24,energy24x);m24=mean(energy24);


%% Determining Threshold Using The Formula determined analytically by pdf
%  Modelling of TE Operated PWPT coefficients
cfs1t=thresh(n1,r1,m1);
cfs2t=thresh(n2,r2,m2);
cfs3t=thresh(n3,r3,m3);
cfs4t=thresh(n4,r4,m4);
cfs5t=thresh(n5,r5,m5);
cfs6t=thresh(n6,r6,m6);
cfs7t=thresh(n7,r7,m7);
cfs8t=thresh(n8,r8,m8);
cfs9t=thresh(n9,r9,m9);
cfs10t=thresh(n10,r10,m10);
cfs11t=thresh(n11,r11,m11);
cfs12t=thresh(n12,r12,m12);
cfs13t=thresh(n13,r13,m13);
cfs14t=thresh(n14,r14,m14);
cfs15t=thresh(n15,r15,m15);
cfs16t=thresh(n16,r16,m16);
cfs17t=thresh(n17,r17,m17);
cfs18t=thresh(n18,r18,m18);
cfs19t=thresh(n19,r19,m19);
cfs20t=thresh(n20,r20,m20);
cfs21t=thresh(n21,r21,m21);
cfs22t=thresh(n22,r22,m22);
cfs23t=thresh(n23,r23,m23);
cfs24t=thresh(n24,r24,m24);


%% Thresholding Operation.....
coef1s = semisoft_pdf(energy1,cfs1t,n1,r1);
coef2s = semisoft_pdf(energy2,cfs2t,n2,r2);
coef3s = semisoft_pdf(energy3,cfs3t,n3,r3);
coef4s = semisoft_pdf(energy4,cfs4t,n4,r4);
coef5s = semisoft_pdf(energy5,cfs5t,n5,r5);
coef6s = semisoft_pdf(energy6,cfs6t,n6,r6);
coef7s = semisoft_pdf(energy7,cfs7t,n7,r7);
coef8s = semisoft_pdf(energy8,cfs8t,n8,r8);
coef9s = semisoft_pdf(energy9,cfs9t,n9,r9);
coef10s = semisoft_pdf(energy10,cfs10t,n10,r10);
coef11s = semisoft_pdf(energy11,cfs11t,n11,r11);
coef12s = semisoft_pdf(energy12,cfs12t,n12,r12);
coef13s = semisoft_pdf(energy13,cfs13t,n13,r13);
coef14s = semisoft_pdf(energy14,cfs14t,n14,r14);
coef15s = semisoft_pdf(energy15,cfs15t,n15,r15);
coef16s = semisoft_pdf(energy16,cfs16t,n16,r16);
coef17s = semisoft_pdf(energy17,cfs17t,n17,r17);
coef18s = semisoft_pdf(energy18,cfs18t,n18,r18);
coef19s = semisoft_pdf(energy19,cfs19t,n19,r19);
coef20s = semisoft_pdf(energy20,cfs20t,n20,r20);
coef21s = semisoft_pdf(energy21,cfs21t,n21,r21);
coef22s = semisoft_pdf(energy22,cfs22t,n22,r22);
coef23s = semisoft_pdf(energy23,cfs23t,n23,r23);
coef24s = semisoft_pdf(energy24,cfs24t,n24,r24);

% coef10s=energy10;

%% Reconstruction Process Starts here.......

coef11_re((len11(1)+1):len11(2))=coef24s;
coef11_re(1:len11(1))=coef23s;
coef10_re((len(3)+1):len10(4))= waverec(coef11_re,len11m,wv);
coef10_re((len10(2)+1):len10(3))=coef22s;
coef10_re((len10(1)+1):len10(2))=coef21s;
coef10_re(1:len10(1))=coef20s;
coef_re((s_no(6)+1):s_no(7))= waverec(coef10_re,len10m,wv);
coef9_re(1:len9(1))=coef17s;
coef9_re((len9(1)+1):len9(2))=coef18s;
coef9_re((len9(2)+1):len9(3))=coef19s;
coef7_re((len7(3)+1):len7(4))= waverec(coef9_re,len9m,wv);
coef8_re(1:len8(1))=coef15s;
coef8_re((len8(1)+1):len8(2))=coef16s;
coef7_re((len7(2)+1):len7(3))= waverec(coef8_re,len8m,wv);
coef7_re(1:len7(1))=coef13s;
coef7_re((len7(1)+1):len7(2))=coef14s;
coef_re(s_no(5):s_no(6))= waverec(coef7_re,len7m,wv);
coef5_re(1:len5(1))=coef11s;
coef5_re((len5(1)+1):len5(2))=coef12s;
coef4_re((len4(2)+1):len4(3))=waverec(coef5_re,len5m,wv);
coef4_re((len4(1)+1):len4(2))=coef10s;
coef4_re(1:len4(1))=coef9s;
coef_re(s_no(4):s_no(5))=waverec(coef4_re,len4m,wv);
coef3_re(1:len3(1))=coef7s;
coef3_re((len3(1)+1):len3(2))=coef8s;
coef2_re((len2(2)+1):len2(3))=waverec(coef3_re,len3m,wv);
coef2_re(1:len2(1))=coef5s;
coef2_re( (len2(1)+1):len2(2) ) =coef6s;
coef_re((s_no(3)+1):s_no(4))=waverec(coef2_re,len2m,wv);
coef1(1:len1(1))=coef3s;
coef1((len1(1)+1):len1(2))=coef4s;
coef_re(s_no(2):s_no(3))=waverec(coef1,len1m,wv);
coef_re(1:s_no(1))=coef1s;
coef_re((s_no(1)+1):s_no(2))=coef2s;

%% Reconstruction of full frame

frames_re=waverec(coef_re,lenm,wv);

end




