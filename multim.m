%##########################################################################
%File: multim.m                                                           #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 13, 2013                                           #
%Description: Calculates mutual inductance for multiple files previously  #
%generated and plots a 3D surface of (emitter radius, receiver radius and #
%mutual inductance.                                                       #
%##########################################################################

clearvars -except Mprt recEmiRatio
clc
global Mprt

m = 0;                  %Contadores Auxiliares.
n = 0;
k = 0;
p = 0;

[nRow,nCol] = size(Mprt);

for m = 1:nCol
    for n = 1:nCol
   
        openfemm;
        str1 = sprintf('geometry_%d%d.fem',m,n);
        str2 = sprintf('/geometry_%d%d.fem',m,n);
        try
            vv=ver;
            opendocument([cd,str2]);  %Abre o arquivo modelo, desenvolvdido pelo autor.
        catch
            opendocument(str1);
        end
        
%         mi_selectsegment(0,height/2);
%         mi_selectsegment(length,height/2);
%         mi_selectsegment(length/2,height);
%         mi_selectsegment(length/2,-height);
%         mi_setsegmentprop('<none>',0,1,0,0);
%         mi_clearselected;
        
        mi_analyze;
        mi_loadsolution;
        r=mo_getcircuitproperties('receiver');
        emiOpenCalc = abs(r(2));
        
        v = mo_getprobleminfo;
        frequency = v(2);
        
        u=mo_getcircuitproperties('emitter');
        mutualCalc(m,n) = abs(emiOpenCalc/(2*pi*frequency*u(1)));
        emiCalc = abs(u(3)/u(1));
              
        closefemm;
        str3 = sprintf('geometry_%d%d.ans',m,n);
        delete(str3);
    end
end

for k = 1:nCol    
    emiRadius(k) = Mprt(8,k);   
    recRadius(k) = Mprt(11,k);
      
end

maxMutual = max(max(mutualCalc));
mutualCalc = mutualCalc/maxMutual;      %Normaliza indutâncoa mútua.

surf(2*emiRadius,2*recRadius,mutualCalc) %Plota considerando os diâmetros.
colormap jet
colorbar

