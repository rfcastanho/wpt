
clear; clc;

Voc = 300; Sutg = 400; Lp = 50e-6;
Vin = 200; w = 2*pi*38400;
Ip = 21.21; Isctg = 6;
Rs = 2.0;
n = 0;
coluna = 0; linha = 0;

% for coupling = 0.0025:2.5e-3:0.5
%     
%     coluna = 0;
%     linha = linha + 1;
%     k(linha) = coupling;
%     
%     for emiInd = 1e-6:1e-6:200e-6
%         
%         coluna = coluna + 1;
%         Lp(coluna) = emiInd;
%         Ls(linha,coluna) = Lp(coluna)*Voc^2/(k(linha)^2*Vin^2);
%     end
%     
%     
%     
% end
% 
% Lpalt = Lp*Voc^2/Vin^2;
% %k = k.^2;
% [K LP] = meshgrid(k,Lpalt);
% LS = LP*pinv(K.^2);
% colormap(jet);
% surf(K,LP,LS);

kMin = sqrt(Sutg/(w*Lp*(Ip/sqrt(2))^2));  %This is the minimum k that satisfies Sutg.
cV = Lp*Voc^2/(Vin^2*kMin^2);             %For Voc to be achieved, Ls must be greater or equal to cV;

for Ls = 1e-6:1e-6:1000e-6
   
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





