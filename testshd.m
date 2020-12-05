%teste shield:

clear all
clc

x2 = 3.3;                                    %Posições Iniciais dos dos enrolamentos.
y1 = -8.05;                                  %Estas coordenadas devem-se ao arquivo inicial desenvolvido
x2 = 3.3;                                    %no FEMM.
y2 = 41.95;
yfill2 = -8.2375;
y2_limite_inf = 41.7625;
xfill2 = 3.3;
yfill2 = 45.8625;
xfill1 = 3.3;
yfill1 = -4.1375;
shieldDist = 0.5;


openfemm;                                    %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/ax_magnetico.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('ax_magnetico.fem');
end
mi_saveas('temp_ind.fem');                   %Aqui salva como um arquivo temporário, para não
             
recShdDesired = input('Digite (1) para usar blindagem no receptor ou (0) para não blindar:');

if recShdDesired == 1

    recMatName1 = input('Qual o material da blindagem interna?: ', 's');
    recThicknessShield1 = input('Qual a espessura da blindagem (mm)?:');
    
mi_addnode(0,yfill2+shieldDist);
mi_addnode(0,yfill2+recThicknessShield1+shieldDist);
mi_addnode(x2,yfill2+shieldDist);
mi_addnode(x2,yfill2+recThicknessShield1+shieldDist);
mi_addsegment(0,yfill2+shieldDist,x2,yfill2+shieldDist);
mi_addsegment(0,yfill2+recThicknessShield1+shieldDist,x2,yfill2+recThicknessShield1+shieldDist);
mi_addsegment(x2,yfill2+recThicknessShield1+shieldDist,x2,yfill2+shieldDist);
mi_addblocklabel(x2/2,yfill2+shieldDist+(recThicknessShield1)/2);
mi_selectlabel(x2/2,yfill2+shieldDist+(recThicknessShield1)/2);
mi_setblockprop(recMatName1,0,0.5,'<None>',0,4,0);
mi_clearselected;
mi_selectnode(0,yfill2+shieldDist);
mi_selectnode(0,yfill2+recThicknessShield1+shieldDist);
mi_selectnode(x2,yfill2+shieldDist);
mi_selectnode(x2,yfill2+recThicknessShield1+shieldDist);
mi_selectsegment(x2/2,yfill2+shieldDist);
mi_selectsegment(x2/2,yfill2+recThicknessShield1+shieldDist);
mi_selectsegment(x2,(yfill2+recThicknessShield1+shieldDist+(yfill2-shieldDist))/2);
mi_setgroup(4);
mi_clearselected;

    recShdDesired2 = input('Digite (1) para adicionar segunda blindagem ao receptor ou (0) para não adicionar:');
    if recShdDesired2 == 1
        
        emiMatName2 = input('Qual o material da blindagem externa?: ', 's');
        recThicknessShield2 = input('Qual a espessura da segunda blindagem (mm)?:');
        
        mi_selectnode(0,yfill2+recThicknessShield1+shieldDist);
        mi_selectnode(x2,yfill2+recThicknessShield1+shieldDist);
        mi_selectlabel(x2/2,((yfill2+recThicknessShield1+shieldDist)+(yfill2+shieldDist))/2);
        mi_copytranslate(0,recThicknessShield2,1);
        mi_clearselected;
        mi_addblocklabel(x2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_selectlabel(x2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_setblockprop(emiMatName2,0,0.5,'<None>',0,4,0);
        mi_addsegment(0,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,x2,yfill2+recThicknessShield1+shieldDist+recThicknessShield2);
        mi_addsegment(x2,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,x2,yfill2+recThicknessShield1+shieldDist);
        mi_selectsegment(x2/2,yfill2+shieldDist+recThicknessShield1+recThicknessShield2);
        mi_selectsegment(x2,(yfill2+recThicknessShield1+recThicknessShield2+shieldDist+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_setgroup(4);
        mi_clearselected;
    end
    
end
