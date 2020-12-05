%##########################################################################
%File: multidwr.m                                                         #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 13, 2013                                           #
%Description: Generates multiple .fem files containing all combinations of#
%emitter radius and receiver radius (saved in Mprt matrix)for a given     #
%separation distance.                                                     #
%##########################################################################


clearvars -except Mprt
clc
global Mprt

k = 0;          %Contadores auxiliares.
p = 0;

openfemm;                                         %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/blank.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('blank.fem');
end
%%
disp('**** Definição do Domínio de Simulação ****');

height = input('Qual a altura da janela de simulação (mm)?:');
length = input('Qual o comprimento da janela de simulação (mm)?:');

mesh = 1;                      %Não está usando.

mi_addnode(0,-height);
mi_addnode(0,height);
mi_addnode(length,height);
mi_addnode(length,-height);
mi_addblocklabel(0.05*length,0.95*height);
mi_selectlabel(0.05*length,0.95*height);
mi_setblockprop('Air',1,mesh,0,0,0,0)
mi_clearselected;
mi_addsegment(0,-height,0,height);
mi_addsegment(0,-height,length,-height);
mi_addsegment(length,-height,length,height);
mi_addsegment(length,height,0,height);
%mi_selectsegment(0,height/2);                      %Não coloca condição de contorno para a linha de simetria.
mi_selectsegment(length,height/2);
mi_selectsegment(length/2,height);
mi_selectsegment(length/2,-height);
mi_setsegmentprop('<None>',0,1,0,0);
mi_clearselected;

%%
disp('**** Separação Emissor-Receptor ****');
dist = input('Qual a distância entre emissor e receptor (mm)?:');

%% Definição de Geometrias

[nRow,nCol] = size(Mprt);

for k = 1:nCol
    
    frequency =     Mprt(1,k);
    frequency =     frequency*1e3;
    wirep =         Mprt(2,k);
    nlitzp =        Mprt(3,k);
    emiTurn =       Mprt(6,k);
    emiCoilHeight = Mprt(7,k);
    emiRadius =     Mprt(8,k);
    
    mi_addnode(emiRadius,0);
    mi_addnode(emiRadius,emiCoilHeight);
    mi_addarc(emiRadius,0,emiRadius,emiCoilHeight,180,1);
    mi_addarc(emiRadius,emiCoilHeight,emiRadius,0,180,1);
    mi_addblocklabel(emiRadius,emiCoilHeight/2);

    rwirep = (8.2514694*exp(-0.115943*wirep))/2;  
    mi_addmaterial('Litz Emissor',1,1,0,0,58,0,0,1,5,0,0,nlitzp,2*rwirep); %Cria condutor Litz com as especificações da iteração atual.

    mi_selectlabel(emiRadius,emiCoilHeight/2);
    mi_setblockprop('Litz Emissor',0,0.1,0,0,0,0)
    mi_clearselected;

    mi_selectlabel(emiRadius,emiCoilHeight/2);
    mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,emiTurn);  
    mi_clearselected;

        for p = 1:nCol

            wires =         Mprt(4,p);
            nlitzs =        Mprt(5,p);
            recTurn =       Mprt(9,p);
            recCoilHeight = Mprt(10,p);
            recRadius =     Mprt(11,p);
            
            mi_addnode(recRadius,emiCoilHeight+dist);
            mi_addnode(recRadius,emiCoilHeight+dist+recCoilHeight);
            mi_addarc(recRadius,emiCoilHeight+dist,recRadius,emiCoilHeight+dist+recCoilHeight,180,1);
            mi_addarc(recRadius,emiCoilHeight+dist+recCoilHeight,recRadius,emiCoilHeight+dist,180,1);
            mi_addblocklabel(recRadius,emiCoilHeight+dist+recCoilHeight/2);

            rwires = (8.2514694*exp(-0.115943*wires))/2;  
            mi_addmaterial('Litz Receptor',1,1,0,0,58,0,0,1,5,0,0,nlitzs,2*rwires); %Cria condutor Litz com as especificações desejadas.

            mi_selectlabel(recRadius,emiCoilHeight+dist+recCoilHeight/2);
            mi_setblockprop('Litz Receptor',0,0.1,0,0,0,0)
            mi_clearselected;
            mi_selectlabel(recRadius,emiCoilHeight+dist+recCoilHeight/2);
            mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,1,recTurn);  
            mi_clearselected;

            mi_probdef(frequency,'millimeters','axi',1e-8,0,20,(0));
            fileName = sprintf('geometry_%d%d.fem',k,p);
            mi_saveas(fileName);
            recEmiRatio(k,p) = recRadius/emiRadius; 
            
            mi_selectnode(recRadius,emiCoilHeight+dist);
            mi_selectnode(recRadius,emiCoilHeight+dist+recCoilHeight);
            mi_deleteselectednodes;
            mi_selectlabel(recRadius,emiCoilHeight+dist+recCoilHeight/2);
            mi_deleteselectedlabels;
        end
    
    mi_selectnode(emiRadius,0);
    mi_selectnode(emiRadius,emiCoilHeight);
    mi_deleteselectednodes;
    mi_selectlabel(emiRadius,emiCoilHeight/2);
    mi_deleteselectedlabels;        
        
end

closefemm;

