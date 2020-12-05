%Equations for LCL S Compensation
%Note: Equations were test and compared with Pspice LCL_S file.
%Results perfectly agree.

clear; clc;
Vp = 200 + 1i*0;
w = 2*pi*38400;
Lp = 120e-6;
Ls = 190e-6;
Rp = 0.3;
Rs = 0.5;
Cp = 32.813e-9;
Cs = 90.412e-9;
ZL = 4;
%k = 0.4;
%M = k*sqrt(Lp*Ls);
M = 60.4e-6;
La = 200e-6;
rint = 10e-3;

Za = rint + 1i*w*La;
t = ((Za^2)/(1/(1i*w*Cp)+ Za) - Za);
c1 = Rp + 1i*w*(Lp-M)+1i*w*M - t;
m = 1 - (Za/(1/(1i*w*Cp)+ Za));
c2 = Rs + 1i*w*(Ls-M)+1/(1i*w*Cs)+ZL+1i*w*M;
Ip = (m*Vp)/(c1 + (w^2*M^2)/c2);


Is = 1i*w*M*Ip/(c2);
disp('Peak Is value (A):');
disp(abs(Is));

disp('Peak Ip value (A):');
disp(abs(Ip));

IL = Is;

SL = ZL*IL*conj(IL);
disp('Load Apparent Power (VA):');
disp(SL);



