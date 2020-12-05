%Equations for PP Compensation
%Note: Equations were test and compared with Pspice TmodelPP file.
%Results perfectly agree.

clear; clc;
Vp = 200 + 1i*0;
w = 2*pi*40000;
rint = 4.3;
Lp = 120e-6;
Ls = 80e-6;
Rp = 1.29;
Rs = 0.59;
Cp = 0.31e-6;
Cs = 0.46e-6;
ZL = 4 + 1i*0;
%k = 0.4;
%M = k*sqrt(Lp*Ls);
M = 52e-6;
t = ((rint^2)/(1/(1i*w*Cp)+ rint) - rint);
n = 1 - (rint/(1/(1i*w*Cp)+ rint));
c = Rp + 1i*w*(Lp-M)+1i*w*M - t;
d = Rs + 1i*w*(Ls-M)+ZL/(1+1i*w*ZL*Cs)+1i*w*M;
Ip = (n*Vp)/(c + (w^2*M^2)/d);
disp('Peak Ip value (A):');
disp(abs(Ip));

Is = 1i*w*M*Ip/(d);
disp('Peak Is value (A):');
disp(abs(Is));

Vbc = (ZL*Is)/(1 + 1i*w*ZL*Cs);
IL = Is - 1i*w*Cs*Vbc;

SL = ZL*IL*conj(IL);
disp('Load Apparent Power (VA):');
disp(SL);