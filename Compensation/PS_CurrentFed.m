%Equations for PS Compensation with Current AC Source
%Note: Equations were test and compared with Pspice TmodelPS_CurrentFed file.
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
Zm = 1i*w*M;
beta =  Rs + 1i*w*(Ls-M) + 1/(1i*w*Cs)+ZL;
a = ((1 + beta/Zm)*(Isource/(1i*w*Cp)))/beta;

Ip = a/(1 + ((1 + beta/Zm)*(1/(1i*w*Cp)+Zp))/beta);

disp('Peak Ip value (A):');
disp(abs(Ip));

Is = Ip/(1 + beta/Zm);
disp('Peak Is value (A):');
disp(abs(Is));

IL = Is;

SL = ZL*IL*conj(IL);
disp('Load Apparent Power (VA):');
disp(SL);




