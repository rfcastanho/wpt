%Equations for PP Compensation with Current AC Source
%Note: Equations were test and compared with Pspice TmodelPP_CurrentFed file.
%Results perfectly agree.

clear; clc;
Isource = 10 +1i*0;
w = 2*pi*40e3;
Lp = 120e-6;
Ls = 80e-6;
Rp = 1.29;
Rs = 0.57;
Cp = 0.31e-6;
Cs = 0.46e-6;
ZL = 4;
%k = 0.4;
%M = k*sqrt(Lp*Ls);
M = 52e-6;


Zp = Rp + 1i*w*(Lp-M);
Zs = Rs + 1i*w*(Ls-M);
Zm = 1i*w*M;
r = (1/ZL + 1i*w*Cs);
m = (1 + 1/(1i*w*Cp*Zm) + Zp/Zm);
n = (r/(1i*w*Cp) + Zp*r)/(1 + Zs*r);

Ip = Isource*(((r/(1i*w*Cp))/(1 + Zs*r))+ 1/(1i*w*Cp* Zm))/(m+n);

disp('Peak Ip value (A):');
disp(abs(Ip));

Is = (Isource*r/(1i*w*Cp) - (r/(1i*w*Cp) + Zp*r)*Ip)/(1 + Zs*r);
disp('Peak Is value (A):');
disp(abs(Is));

Im = Ip - Is;
Vcs = Zm*Im - Zs*Is;
IL = Vcs/ZL;

SL = ZL*IL*conj(IL);
disp('Load Apparent Power (VA):');
disp(SL);




