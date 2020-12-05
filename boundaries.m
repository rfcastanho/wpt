
clear; clc;

Voc = 60; Sutg = 400; Lp = 50e-6;
Vin = 250; w = 2*pi*38400;
Ip = 21.2; Isctg = 6;
Rs = 1.0;
n = 0;


kMin = sqrt(Sutg/(w*Lp*(Ip/sqrt(2))^2));  %This is the minimum k that satisfies Sutg.
cV = Lp*Voc^2/(Vin^2*kMin^2);             %For Voc to be achieved, Ls must be greater or equal to cV;

for Ls = 1e-6:1e-6:150e-6
   
    n = n + 1;
    ind(n) = Ls;
    M(n) = sqrt(Sutg*Ls/(w*(Ip/sqrt(2))^2));
    
    k(n) = M(n)/(sqrt(Lp*Ls));
    
    cS(n) = sqrt(Sutg*Ls/(w*(Ip/sqrt(2))^2));  %For Su_target to be achieved, M must be greater or equal to cS;
    cI(n) = abs(Isctg*(Rs + 1i*w*Ls)/(1i*w*Ip));    %For Isc_target to be achieved, M must be greater or equal to cI;
    
end

plot(ind,cS,'r--')
hold on
plot(ind,cI,'k--')

cVvecY = [0 1.5*max(M)]; cVvecX = [cV cV];
plot(cVvecX,cVvecY,'b--')





