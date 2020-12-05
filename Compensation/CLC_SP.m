%Equations for CLC SP Compensation
%Note: Equations were test and compared with Pspice CLC_SP file.
%Results perfectly agree.

clear; clc;
Vp = 200 + 1i*0;
w = 2*pi*38400;
Lp = 120e-6;
Ls = 190e-6;
Rp = 0.3;
Rs = 0.5;
Cp = 0.170e-6;
Cs = 90.412e-9;
ZL = 200;
%k = 0.4;
%M = k*sqrt(Lp*Ls);
M = 60.4e-6;
La = 200e-6;
Ca = 100e-9;
Cb = 100e-9;
Zp = Rp + 1i*w*(Lp-M) + 1/(1i*w*Cp);
Zs = Rs + 1i*w*(Ls-M) + 1/(1i*w*Cs);
Zm = 1i*w*M;

t = 1 + (1i*w*La)/(Zp*(1-w^2*Cb*La));
b = w^2*Cb*La*Vp/(Zp*(1-w^2*Cb*La));
a = (Zs - 1/(1i*w*Cs - w^2*Cs^2*ZL))/Zm;

Is = (Vp/(Zp*t) + b/t)/(1 + a + Zm*a/(Zp*t));

disp('Peak Is value (A):');
disp(abs(Is));

Ip = (1 + a)*Is;

disp('Peak Ip value (A):');
disp(abs(Ip));

IL = Is/(1 + 1i*w*ZL*Cs);

SL = ZL*IL*conj(IL);
disp('Load Apparent Power (VA):');
disp(SL);



