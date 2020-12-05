%##########################################################################
%File: drawer.m                                                           #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 13, 2013                                           #
%Description: Planar and Axissymmetric fast drawer.                       #
%##########################################################################

clear
clc

lateralDistance = 0;
angles = 0;
windType = 0;

openfemm;                                         %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/blank.fem']);  %Abre arquivo em branco.
catch
    opendocument('blank.fem');
end

disp('**** Definição do Tipo de Simulação ****');

planarSim = input ('Digite 0 para simulação axisimétrica ou 1 para planar:');

disp('**** Definição do Domínio de Simulação ****');

height = input('Qual a altura da janela de simulação (mm)?:');
length = input('Qual o comprimento da janela de simulação (mm)?:');

mesh = 1;                      %Não está usando.

if planarSim == 0
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
    mi_selectsegment(0,height/2);
    mi_selectsegment(length,height/2);
    mi_selectsegment(length/2,height);
    mi_selectsegment(length/2,-height);
    mi_setsegmentprop('<None>',0,1,0,0);
    mi_clearselected;
    
    disp('**** Definição do Tipo de Enrolamento ****');

    windType = input('Enrolamento (1) Concentrado ou (2) "Panqueca":');

    disp('**** Definição do Emissor ****');

    if windType == 1
        rp = input('Raio do emissor (mm):');
    elseif windType == 2
        rp = input('Raio interno do emissor (mm):');  
    end

elseif planarSim ==1
    mi_addnode(-length/2,-height/2);
    mi_addnode(-length/2,height/2);
    mi_addnode(length/2,height/2);
    mi_addnode(length/2,-height/2);
    mi_addblocklabel(0.95*(-length/2),0.95*height/2);
    mi_selectlabel(0.05*(-length/2),0.95*height/2);
    mi_setblockprop('Air',1,mesh,0,0,0,0)
    mi_clearselected;
    mi_addsegment(-length/2,-height/2,-length/2,height/2);
    mi_addsegment(-length/2,-height/2,length/2,-height/2);
    mi_addsegment(length/2,-height/2,length/2,height/2);
    mi_addsegment(length/2,height/2,-length/2,height/2);
    mi_selectsegment(-length/2,0);
    mi_selectsegment(length/2,0);
    mi_selectsegment(0,height/2);
    mi_selectsegment(0,-height/2);
    mi_setsegmentprop('<None>',0,1,0,0);
    mi_clearselected;     
    
    disp('**** Definição do Emissor ****');
    rp = input('Raio do emissor (mm):');
    
end

Ep = input('Espessura do emissor (mm):');
wirep = input('Condutor do emissor (AWG):');
nlitzp = input('Número de condutores em paralelo:');
emiTurn = input('Número de espiras:');

mi_addnode(rp,0);
mi_addnode(rp,Ep);
mi_addarc(rp,0,rp,Ep,180,1);
mi_addarc(rp,Ep,rp,0,180,1);
mi_addblocklabel(rp,Ep/2);

rwirep = (8.2514694*exp(-0.115943*wirep))/2;  
mi_addmaterial('Litz Emissor',1,1,0,0,58,0,0,1,5,0,0,nlitzp,2*rwirep); %Cria condutor Litz com as especificações da iteração atual.

mi_selectlabel(rp,Ep/2);
mi_setblockprop('Litz Emissor',0,0.1,0,0,0,0)
mi_clearselected;

if windType == 1
    mi_selectlabel(rp,Ep/2);
    mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,emiTurn);  
    mi_clearselected;

elseif windType == 2
    mi_selectlabel(rp,Ep/2);
    mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,1);  
    mi_clearselected;
    
    mi_selectnode(rp,Ep);
    mi_selectnode(rp,0);
    mi_setgroup(1);
    mi_clearselected;
    
    mi_selectarcsegment(rp-1,Ep/2);
    mi_selectarcsegment(rp+1,Ep/2);
    mi_setgroup(1);
    mi_clearselected;
    
    mi_selectlabel(rp,Ep/2);
    mi_setgroup(1);
    mi_clearselected;
    
    mi_selectgroup(1);
    mi_copytranslate(Ep+0.1,0,emiTurn-1);
    mi_clearselected;
    
end

if (planarSim == 1)

    mi_addnode(-rp,0);
    mi_addnode(-rp,Ep);
    mi_addarc(-rp,0,-rp,Ep,180,1);
    mi_addarc(-rp,Ep,-rp,0,180,1);
    mi_addblocklabel(-rp,Ep/2);

    mi_selectlabel(-rp,Ep/2);
    mi_setblockprop('Litz Emissor',0,0.1,0,0,0,0)
    mi_clearselected;

    mi_selectlabel(-rp,Ep/2);
    mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,-emiTurn);  
    mi_clearselected;

end

windType = 0;
clc;
%%
disp('**** Separação Emissor-Receptor ****');
dist = input('Qual a distância entre emissor e receptor (mm)?:');

%%
disp('**** Definição do Receptor ****');

if planarSim == 0
    
     disp('**** Definição do Tipo de Enrolamento ****');

    windType = input('Enrolamento (1) Concentrado ou (2) "Panqueca":');

    if windType == 1
        rs = input('Raio do receptor (mm):');
    elseif windType == 2
        rs = input('Raio interno do receptor (mm):');  
    end
    
elseif planarSim == 1
    
    rs = input('Raio do receptor (mm):');
    
end

Es = input('Espessura do receptor (mm):');
wires = input('Condutor do receptor (AWG):');
nlitzs = input('Número de condutores em paralelo:');
recTurn = input('Número de espiras:');
mi_addnode(rs,Ep+dist);
mi_addnode(rs,Ep+dist+Es);
mi_addarc(rs,Ep+dist,rs,Ep+dist+Es,180,1);
mi_addarc(rs,Ep+dist+Es,rs,Ep+dist,180,1);
mi_addblocklabel(rs,Ep+dist+Es/2);

rwires = (8.2514694*exp(-0.115943*wires))/2;  
mi_addmaterial('Litz Receptor',1,1,0,0,58,0,0,1,5,0,0,nlitzs,2*rwires); %Cria condutor Litz com as especificações desejadas.

mi_selectlabel(rs,Ep+dist+Es/2);
mi_setblockprop('Litz Receptor',0,0.1,0,0,0,0)
mi_clearselected;

if windType == 1
    mi_selectlabel(rs,Ep+dist+Es/2);
    mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recTurn);  
    mi_clearselected;
    mi_selectnode(rs,Ep+dist);
    mi_selectnode(rs,Ep+dist+Es);
    mi_setgroup(2);
    mi_clearselected;
    
    mi_selectarcsegment(rs-1,(Ep+dist+Es)/2);
    mi_selectarcsegment(rs+1,(Ep+dist+Es)/2);
    mi_setgroup(2);
    mi_clearselected;

elseif windType == 2

    mi_selectlabel(rs,Ep+dist+Es/2);
    mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recTurn);  
    mi_clearselected;
    mi_selectnode(rs,Ep+dist);
    mi_selectnode(rs,Ep+dist+Es);
    mi_setgroup(2);
    mi_clearselected;
    
    mi_selectarcsegment(rs-1,(Ep+dist+Es)/2);
    mi_selectarcsegment(rs+1,(Ep+dist+Es)/2);
    mi_setgroup(2);
    mi_clearselected;
    
%     mi_selectlabel(rs,Ep+dist+Es/2);
%     mi_setgroup(2);
%     mi_clearselected;
    
    mi_selectgroup(2);
    mi_copytranslate(Es+0.1,0,recTurn-1);
    mi_clearselected;
    
end

if (planarSim == 1)
    mi_addnode(-rs,Ep+dist);
    mi_addnode(-rs,Ep+dist+Es);
    mi_addarc(-rs,Ep+dist,-rs,Ep+dist+Es,180,1);
    mi_addarc(-rs,Ep+dist+Es,-rs,Ep+dist,180,1);
    mi_addblocklabel(-rs,Ep+dist+Es/2);

    mi_selectlabel(-rs,Ep+dist+Es/2);
    mi_setblockprop('Litz Receptor',0,0.1,0,0,0,0)
    mi_clearselected;

    mi_selectlabel(-rs,Ep+dist+Es/2);
    mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,1,-recTurn);  
    mi_clearselected;
    
    angles = input('Qual a inclinação entre emissor e receptor (em graus):?');
    lateralDistance = input('Qual a distância entre centros do emissor e receptor (mm)?:');    
    
    mi_selectnode(-rs,Ep+dist);
    mi_selectnode(-rs,Ep+dist+Es);
    mi_selectnode(rs,Ep+dist);
    mi_selectnode(rs,Ep+dist+Es);
    mi_setgroup(2);
    mi_clearselected;
    
    mi_selectarcsegment(-rs-Es/2,Ep+dist);
    mi_selectarcsegment(-rs+Es/2,Ep+dist);
    mi_selectarcsegment(rs-Es/2,Ep+dist);
    mi_selectarcsegment(rs+Ep/2,Ep+dist);
    mi_setgroup(2);
    mi_clearselected;
    
    mi_selectlabel(-rs,Ep+dist+Es/2);
    mi_selectlabel(rs,Ep+dist+Es/2);
    mi_setgroup(2);
    mi_clearselected;
    
    mi_selectgroup(2);
    mi_movetranslate(lateralDistance,0);
    mi_clearselected;
    
    mi_selectgroup(2);
    mi_moverotate(-rs+lateralDistance,Ep+dist+Es/2,angles);
    mi_clearselected;
 end

frequency = input ('Freqüência de operação (kHz):');
frequency = frequency*1e3;

if planarSim == 0
    mi_probdef(frequency,'millimeters','axi',1e-8,0,20,(0));

elseif planarSim == 1
    mi_probdef(frequency,'millimeters','planar',1e-8,2*rp,20,(0));
    
end





projectName = input('Digite o nome do projeto: ', 's');
fileName = sprintf('geometry_%s.fem',projectName);
mi_saveas(fileName);
closefemm;

