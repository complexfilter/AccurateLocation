[a,b]=size(I);

load('filter_pool/hh_high1.mat');
load('filter_pool/hh_high2.mat');
load('filter_pool/hh_high3.mat');
load('filter_pool/hh_low3.mat');

N1=length(hh_high1);
i=0:N1-1;

phase_A3=exp(-1j*pi/2*(i-(N1-1)/2));
phase_A2=exp(-1j*pi*(i-(N1-1)/2));
phase_A1=exp(-1j*pi*3/2*(i-(N1-1)/2));

hh_high_A1=phase_A1.*hh_high1;
hh_high_A2=phase_A2.*hh_high1;
hh_high_A3=phase_A3.*hh_high1;
hh_high_A4=hh_high1;

N2=length(hh_high2);
i=0:N2-1;

phase_B3=exp(-1j*pi/4*(i-(N2-1)/2));
phase_B2=exp(-1j*pi/2*(i-(N2-1)/2));
phase_B1=exp(-1j*3*pi/4*(i-(N2-1)/2));

hh_high_B1=phase_B1.*hh_high2;
hh_high_B2=phase_B2.*hh_high2;
hh_high_B3=phase_B3.*hh_high2;
hh_high_B4=hh_high2;

N3=length(hh_high3);
i=0:N3-1;

phase_C1=exp(-1j*3*pi/8*(i-(N3-1)/2));
phase_C2=exp(-1j*pi/4*(i-(N3-1)/2));
phase_C3=exp(-1j*pi/8*(i-(N3-1)/2));
phase_C5=exp(+1j*pi/8*(i-(N3-1)/2));


hh_high_C1=phase_C1.*hh_high3;
hh_high_C2=phase_C2.*hh_high3;
hh_high_C3=phase_C3.*hh_high3;
hh_high_C4=hh_high3;
hh_high_C5=phase_C5.*hh_high3;


N3=length(hh_high3);
i=0:N3-1;
load('filter_pool/hh_low3.mat');
phase_C00=exp(-1j*pi/23*(i-(N3-1)/2));
hh_high_C00=hh_high3.*phase_C00;
%fvtool(hh_high_C00);


A1 = cconv2(-hh_high_A2, hh_high_A4,I);
A2 = cconv2( hh_high_A4, hh_high_A1,I);
A3 = cconv2( hh_high_A4,-hh_high_A2,I);
A4 = cconv2( hh_high_A4,-hh_high_A3,I);
A5 = cconv2( hh_high_A4, hh_high_A4,I);
A6 = cconv2(-hh_high_A2, hh_high_A1,I);

B1 = cconv2(-hh_high_B2, hh_high_B4,I);
B2 = cconv2( hh_high_B4, hh_high_B1,I);
B3 = cconv2( hh_high_B4,-hh_high_B2,I);
B4 = cconv2( hh_high_B4,-hh_high_B3,I);
B5 = cconv2( hh_high_B4, hh_high_B4,I);
B6 = cconv2(-hh_high_B2, hh_high_B1,I);

C0 = cconv2(hh_high_C2,hh_high_C3,I);
C1 = cconv2(hh_high_C2,hh_high_C4,I);
C2 = cconv2(hh_high_C4,hh_high_C1,I);
C3 = cconv2(hh_high_C4,hh_high_C2,I);
C4 = cconv2(hh_high_C4,hh_high_C3,I);
C5 = cconv2(hh_high_C4,hh_high_C4,I);
C6 = cconv2(hh_high_C2,hh_high_C1,I);
C7 = cconv2(hh_high_C2,hh_high_C2,I);

C00 = cconv2(hh_high_C2,hh_high_C00,I);


C52= cconv2(hh_high_C4,hh_high_C5,I);
