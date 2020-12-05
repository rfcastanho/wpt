%##########################################################################
%File: padfix.m                                                           #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Axissymmetric representation of distributed coils (pancake) #
%with maximum diameter restriction.                                       #
%##########################################################################


function panc_dfix = padfix(fileName,logDesired,frequency,emiMaxDiameter,recMaxDiameter,emiInduct,recInduct,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent)


if nargin <19 
    error('Faltam argumentos de entrada')
end


% clear
% clc
% 
% emiInduct = 166.5*1e-6;
% recInduct = 46.5*1e-6;
% emiMaxDiameter = 200;
% recMaxDiameter = 140;
% x1 = 3.3;
% x2 = 3.3;
% y1 = -8.05;
% y2 = 41.95;
% y1_limite_inf = -8.2375;
% y2_limite_inf = 41.7625;
% xfill2 = 3.3;
% yfill2 = 43.8455;
% xfill1 = 3.3;
% yfill1 = -5.7392;


iterCounter1 = 0;                                  %Contador de Iterações 1, inicia com valor nulo.
testCurrent  = 2;                                  %Corrente de teste dos enrolamentos (A);
emiTurn = 1;                                       %Inicia com número de espiras igual a 1.
recTurn = emiTurn;

femTime1 = 0;
mesh = 1;

emiMaxRadius = emiMaxDiameter/2;
recMaxRadius = recMaxDiameter/2;
emiQFactor = 0;
recQFactor = 0;

emiErrorPrev = 1e6;     %Erro anterior para o calculo do emissor. Inicia com valor muito alto.
recErrorPrev = 1e6;

shieldDist = 0.5;
emiShdDesired1 = 0;
emiShdDesired2 = 0;
recShdDesired1 = 0;
recShdDesired2 = 0;
emiThicknessShield1 = 0;
emiThicknessShield2 = 0;
recThicknessShield1 = 0;
recThicknessShield2 = 0;
emiLayerCounter = 1;
recLayerCounter = 1;
emiMaxLayer = 1000;
recMaxLayer = 1000;

x1prev = x1;
x2prev = x2;
usingShield = 0;
maxError = 0;

openfemm;                                    %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/result_wind.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('result_wind.fem');
end
mi_saveas('temp_dfix.fem');                    %Aqui salva como um arquivo temporário, para não
                                               %destruir o orignal.

%% Análise Inicial dos Enrolamentos                                % Esta rotina ajusta o fator de utilização dos enrolamentos
                                                          
mi_selectlabel(x1,y1);                                            
mi_setblockprop('Litz Emissor',0,mesh,'emitter',0,1,emiTurn);
mi_clearselected;

mi_selectlabel(x2,y2);
mi_setblockprop('Litz Receptor',0,mesh,'receiver',0,2,recTurn);
mi_clearselected;
mi_setcurrent('emitter',testCurrent);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',testCurrent); 
mi_analyze;
mi_loadsolution;
emiWindFactor = mo_getfill(x1,y1);
recWindFactor = mo_getfill(x2,y2);

recCoilHeight = (yfill2 - y2_limite_inf);
emiCoilHeight = (yfill1 - y1_limite_inf);

    while ((recWindFactor < 0.7) || (recWindFactor > 0.8))                   %Ajusta espessura do enrolamento sec. até que
                                                         %o fator de utilização seja maior que 70%.
                                                   
        if recWindFactor > 0.8
            windStep = 0.2;
        elseif recWindFactor < 0.7
            windStep = -0.15;
        end
    
        if ((xfill2 - recCoilHeight/2) <= xfill2/10)
            mi_selectgroup(2);
            mi_movetranslate(recCoilHeight,0);
            xfill2 = xfill2 + recCoilHeight;
            x2 = x2 + recCoilHeight;
        end  
    
        mi_selectnode(xfill2,yfill2);
        mi_movetranslate(0,windStep)
        yfill2 = yfill2 + windStep;
        recCoilHeight = (yfill2 - y2_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        recWindFactor = mo_getfill(x2,y2);
    end

    recCoilHeight = (yfill2 - y2_limite_inf);

    while ((emiWindFactor < 0.7) || (emiWindFactor > 0.8))                 %Ajusta espessura do enrolamento prim. até que
                                                       %o fator de utilização seja maior que 70%.  
                                  
        if emiWindFactor > 0.8
            windStep = 0.2;
        elseif emiWindFactor < 0.7
            windStep = -0.15;
        end                          
                                  
        if ((xfill1 - emiCoilHeight/2) <= xfill1/10)
            mi_selectgroup(1);
            mi_movetranslate(emiCoilHeight,0);
            xfill1 = xfill1 + emiCoilHeight;
            x1 = x1 + emiCoilHeight;
        end  
    
        mi_selectnode(xfill1,yfill1);
        mi_movetranslate(0,windStep)
        yfill1 = yfill1 + windStep;
        emiCoilHeight = (yfill1 - y1_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        emiWindFactor = mo_getfill(x1,y1); 
    end

    emiCoilHeight = (yfill1 - y1_limite_inf);

    y1_limite_sup = yfill1;  %Enrolamento sec. não pode chegar nesta altura ou tocará o enrolamento primário.


%%
disp('Iniciando Fase 1 - Geometria do Enrolamento Primário...');

% mi_selectgroup(1);                              %Coloca primário e secundário nos respectivos diâmetros máximos.
% mi_movetranslate(emiRadius-x1,0);
% x1 = emiRadius;
% xfill1 = x1;
% mi_clearselected;

emiRadius = x1;
initialX = x1;

%###################################

%Add Emitter Shield
emiShdDesired1 = input('Digite (1) para usar blindagem no emissor ou (0) para não blindar:');

if emiShdDesired1 == 1

    emiMatName1 = input('Qual o material da blindagem interna?: ', 's');
    emiThicknessShield1 = input('Qual a espessura da blindagem (mm)?:');
    usingShield = 1; %Indicates that ths loosely coupled system has at least one shielding layer.
    
    emiXShd = x1 + emiCoilHeight/2 +emiCoilHeight*0.05;
    
    mi_addnode(0,y1_limite_inf-shieldDist);
    mi_addnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addnode(emiXShd,y1_limite_inf-shieldDist);
    mi_addnode(emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addsegment(0,y1_limite_inf-shieldDist,emiXShd,y1_limite_inf-shieldDist);
    mi_addsegment(0,y1_limite_inf-emiThicknessShield1-shieldDist,emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addsegment(emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist,emiXShd,y1_limite_inf-shieldDist);
    mi_addblocklabel(emiXShd/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_selectlabel(emiXShd/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_setblockprop(emiMatName1,0,0.5,'<None>',0,3,0);
    mi_clearselected;
    mi_selectnode(0,y1_limite_inf-shieldDist);
    mi_selectnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectnode(emiXShd,y1_limite_inf-shieldDist);
    mi_selectnode(emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectsegment(emiXShd/2,y1_limite_inf-shieldDist);
    mi_selectsegment(emiXShd/2,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectsegment(emiXShd,(y1_limite_inf-emiThicknessShield1-shieldDist+(y1_limite_inf-shieldDist))/2);
    mi_setgroup(3);
    mi_clearselected;

        emiShdDesired2 = input('Digite (1) para adicionar segunda blindagem ao emissor ou (0) para não adicionar:');
        if emiShdDesired2 == 1
            emiMatName2 = input('Qual o material da blindagem externa?: ', 's');
            emiThicknessShield2 = input('Qual a espessura da segunda blindagem (mm)?:');
        
            mi_selectnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_selectnode(emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_copytranslate(0,-emiThicknessShield2,1);
            mi_clearselected;
            mi_addblocklabel(emiXShd/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_selectlabel(emiXShd/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_setblockprop(emiMatName2,0,0.5,'<None>',0,3,0);
            mi_addsegment(0,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2,emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2);
            mi_addsegment(emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2,emiXShd,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_selectsegment(emiXShd/2,y1_limite_inf-shieldDist-emiThicknessShield1-emiThicknessShield2);
            mi_selectsegment(emiXShd,(y1_limite_inf-emiThicknessShield1-emiThicknessShield2-shieldDist+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_setgroup(3);
            mi_clearselected;
        end
end
%###################################

%Add Emitter "Winglets"

emiWing = input('Type "1" to add shielding "winglets" or "0" if it is not necessary:');

if emiWing == 1
   %emiWingHeight = input('Enter "winglet" height, in mm:');
   
   if emiShdDesired1 == 1
        mi_selectnode(emiXShd,y1_limite_inf-shieldDist);
        mi_copytranslate(0,shieldDist+emiCoilHeight,1);
        mi_selectnode(emiXShd,shieldDist+emiCoilHeight);
        mi_copytranslate(emiThicknessShield1,0,1);
        mi_selectnode(emiXShd+emiThicknessShield1,shieldDist+emiCoilHeight);
        mi_copytranslate(0,-shieldDist-emiCoilHeight-emiThicknessShield1,1);
        mi_clearselected;
        mi_addsegment(emiXShd,y1_limite_inf-shieldDist,emiXShd,y1_limite_inf+emiCoilHeight);
        mi_addsegment(emiXShd,y1_limite_inf+emiCoilHeight,emiXShd+emiThicknessShield1,y1_limite_inf+emiCoilHeight);
        mi_addsegment(emiXShd+emiThicknessShield1,y1_limite_inf+emiCoilHeight,emiXShd+emiThicknessShield1,y1_limite_inf-emiCoilHeight-shieldDist-emiThicknessShield1);
        mi_addsegment(emiXShd+emiThicknessShield1,y1_limite_inf-emiCoilHeight-shieldDist-emiThicknessShield1,emiXShd,y1_limite_inf-shieldDist-emiThicknessShield1);
        
        mi_selectnode(emiXShd,y1_limite_inf+shieldDist+emiCoilHeight);
        mi_selectnode(emiXShd+emiThicknessShield1,y1_limite_inf+shieldDist+emiCoilHeight);
        mi_selectnode(emiXShd+emiThicknessShield1,y1_limite_inf-shieldDist-emiCoilHeight-emiThicknessShield1);
        mi_selectsegment(emiXShd,y1_limite_inf+emiCoilHeight/2);
        mi_selectsegment(emiXShd+emiThicknessShield1/2,y1_limite_inf+emiCoilHeight);
        mi_selectsegment(emiXShd+emiThicknessShield1,y1_limite_inf+emiCoilHeight/2);
        mi_selectsegment(emiXShd+emiThicknessShield1/2,y1_limite_inf-emiThicknessShield1-shieldDist);
        mi_setgroup(33);
        mi_clearselected;
        
        mi_addblocklabel(emiXShd+emiThicknessShield1/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
        mi_selectlabel(emiXShd+emiThicknessShield1/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
        mi_setblockprop(emiMatName1,0,0.5,'<None>',0,33,0);
        mi_clearselected;
        
   end
   
   if emiShdDesired2 == 1
        mi_selectnode(emiXShd,y1_limite_inf-shieldDist-emiThicknessShield1-emiThicknessShield2);
        mi_copytranslate(emiThicknessShield1+emiThicknessShield2,0,1);
        mi_selectnode(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf-shieldDist-emiThicknessShield1-emiThicknessShield2);
        mi_copytranslate(0,emiThicknessShield1+emiThicknessShield2+shieldDist+emiCoilHeight,1);    
        mi_addsegment(emiXShd,y1_limite_inf-shieldDist-emiXShd-emiThicknessShield1-emiThicknessShield2,emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf-shieldDist-emiXShd-emiThicknessShield1-emiThicknessShield2);
        mi_addsegment(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf-shieldDist-emiXShd-emiThicknessShield1-emiThicknessShield2,emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf+emiCoilHeight);
        mi_addsegment(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf+emiCoilHeight,emiXShd+emiThicknessShield1,y1_limite_inf+emiCoilHeight);
        
        mi_selectnode(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf-shieldDist-emiXShd-emiThicknessShield1-emiThicknessShield2);
        mi_selectnode(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf+emiCoilHeight);     
        mi_selectsegment(emiXShd+emiThicknessShield1+emiThicknessShield2/2,y1_limite_inf-shieldDist-emiXShd-emiThicknessShield1-emiThicknessShield2);
        mi_selectsegment(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf-emiThicknessShield1);
        mi_selectsegment(emiXShd+emiThicknessShield1+emiThicknessShield2/2,y1_limite_inf+emiCoilHeight);
        mi_setgroup(33);
        mi_clearselected;
        
        mi_addblocklabel(emiXShd+emiThicknessShield1+emiThicknessShield2/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
        mi_selectlabel(emiXShd+emiThicknessShield1+emiThicknessShield2/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
        mi_setblockprop(emiMatName2,0,0.5,'<None>',0,33,0);
        mi_clearselected;
       
   end
    
end



disp('It is possible to design windings with minimum diameter but higher coil height. ');
disp('It is also possible to design windings with minimum height but higher coil diameter.');
emiOptim = input('Type "1" to optimize diameter or "0" to optimize height:');

if emiOptim == 1
    emiMaxHeight = input('Enter maximum emitter height (including shields, if any), in mm:');
    emiMaxLayer = floor((emiMaxHeight-emiThicknessShield1-emiThicknessShield2)/emiCoilHeight);

    if emiWing ==1
        mi_selectnode(emiXShd,y1_limite_inf+shieldDist+emiCoilHeight);
        mi_selectnode(emiXShd+emiThicknessShield1,y1_limite_inf+shieldDist+emiCoilHeight);
        mi_selectnode(emiXShd+emiThicknessShield1+emiThicknessShield2,y1_limite_inf+emiCoilHeight); 
        mi_movetranslate(0,emiMaxHeight-(emiThicknessShield1+emiThicknessShield2+emiCoilHeight+shieldDist));
        mi_clearselected;
    end

end


a = 1;
mi_setcurrent('emitter',emiOpenCurrent);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',0); 
mi_analyze;
mi_loadsolution;
r=mo_getcircuitproperties('emitter');          %Primeira estimativa do erro de emiInduct.
emiCalc = abs(r(3)/r(1));   
emiCalcError = abs((emiInduct - emiCalc)/emiInduct)*100;
emiNewTurn = 0;
emiShdTop = y1_limite_inf - shieldDist;

emiClock = tic;
    while (emiCalcError < emiErrorPrev)                %Enquanto erro atual for menor que erro da iteração anterior...
    
        emiErrorPrev = emiCalcError;                   % Salva resultado para erro anterior.
        iterCounter1 = iterCounter1 + 1; 
             
        if emiOptim == 0 || emiLayerCounter == emiMaxLayer    

                if emiShdDesired1 == 1                                            
                    mi_selectnode(emiXShd,emiShdTop);
                    mi_selectnode(emiXShd,emiShdTop-emiThicknessShield1);
                    mi_selectgroup(33);     
                    if emiShdDesired2 == 1
                        mi_selectnode(emiXShd,emiShdTop-emiThicknessShield1-emiThicknessShield2);
                    end  
                    mi_movetranslate((emiCoilHeight+emiCoilHeight*0.01),0);
                    emiXShd = emiXShd + emiCoilHeight + emiCoilHeight*0.01;
                end
                       
            mi_selectnode(x1,yfill1);
            mi_selectnode(x1,y1_limite_inf);
            mi_copytranslate((emiCoilHeight+emiCoilHeight*0.01),0,1);
            mi_selectarcsegment((x1-emiCoilHeight/2),yfill1);
            mi_selectarcsegment((x1+emiCoilHeight/2),yfill1);
            mi_copytranslate2((emiCoilHeight+emiCoilHeight*0.01),0,1,3);
            mi_selectlabel(x1,y1_limite_inf);
            mi_copytranslate2((emiCoilHeight+emiCoilHeight*0.01),0,1,2);
            x1 = x1+(emiCoilHeight+emiCoilHeight*0.01);
            emiRadius = emiRadius + (emiCoilHeight+emiCoilHeight*0.01);
            emiNewTurn = emiNewTurn + 1;
                if emiLayerCounter == emiMaxLayer
                    emiLayerCounter = 1;
                    a = -a;
                else
                    emiLayerCounter = 0;
                end
        elseif emiOptim == 1
            mi_selectnode(x1,yfill1);
            mi_selectnode(x1,y1_limite_inf);
            mi_copytranslate(0,a*(emiCoilHeight+emiCoilHeight*0.01),1);
            mi_selectarcsegment((x1-emiCoilHeight/2),yfill1);
            mi_selectarcsegment((x1+emiCoilHeight/2),yfill1);
            mi_copytranslate2(0,a*(emiCoilHeight+emiCoilHeight*0.01),1,3);
            mi_selectlabel(x1,y1_limite_inf);
            mi_copytranslate2(0,a*(emiCoilHeight+emiCoilHeight*0.01),1,2);
            yfill1 = yfill1+a*(emiCoilHeight+emiCoilHeight*0.01);
            y1_limite_inf = y1_limite_inf+a*(emiCoilHeight+emiCoilHeight*0.01);
            emiNewTurn = emiNewTurn + 1;
            emiLayerCounter = emiLayerCounter + 1;
        end
        
         
        
        if emiRadius > emiMaxRadius
            clc;
            maxError = 1;
            disp('Maximum Diameter reached. Emitter inductance will be used "as is".');
            disp('Press any key to continue');
            pause;
            break
        end

        femClock = tic;
        mi_analyze;                                     %Analisa nova geometria e recomeça.
        mi_loadsolution;
        femTime1(iterCounter1) = toc(femClock);

        r=mo_getcircuitproperties('emitter');
        emiCalc = abs(r(3)/r(1));                                 %Recalcula a indutância emiInduct.
        emiCalcError = abs((emiInduct - emiCalc)/emiInduct)*100                                  
        
        emiVecError(iterCounter1) = ((emiInduct - emiCalc)/emiInduct)*100;  %Expressa erro percentual para gráfico.
        emiVecIter(iterCounter1) = iterCounter1;
        
    end
emiTime = toc (emiClock);
 

%Undo Last Iteration

if maxError == 0
    mi_selectnode(x1,yfill1);
    mi_selectnode(x1,y1_limite_inf);
    mi_selectarcsegment((x1-emiCoilHeight/2),yfill1);
    mi_selectarcsegment((x1+emiCoilHeight/2),yfill1);
    mi_selectlabel(x1,y1_limite_inf);
    mi_deleteselected;

    x1 = x1 - (emiCoilHeight+emiCoilHeight*0.01);
    emiRadius = emiRadius - (emiCoilHeight+emiCoilHeight*0.01);
    emiNewTurn = emiNewTurn - 1;
    iterCounter1 = iterCounter1 - 1;
    emiVecIter(:,length(emiVecIter)) = [];   %Remove last element of this vector.
    emiVecError(:,length(emiVecError)) = [];
    femTime1(:,length(femTime1)) = [];

    mi_analyze(1);
    mi_loadsolution;
    r = mo_getcircuitproperties('emitter');
    emiResistance = real(r(2)/r(1));  
    emiCalc = abs(r(3)/r(1));                                 %Recalcula a indutância emiInduct.
    emiQFactor = (2*pi*frequency*emiCalc)/emiResistance;    
    emiCalcError = abs((emiInduct - emiCalc)/emiInduct)*100; 
end

%Emitter Results
figure;
plot(emiVecIter,emiVecError);
xlabel('Iteration'); ylabel('Error (%)');
title('Emitter Error vs. Iteration');

disp('*************************')
disp('Indutância Primária desejada (H):');
disp(emiInduct);
disp('Indutância obtida (H):');
disp(emiCalc);
disp('Número de espiras enrolamento primário:');
disp(emiNewTurn);
disp('Fator de Qualidade do enrolamento primário:');
disp(emiQFactor);
disp('Número de iterações');
disp(iterCounter1);
disp('*****************************************************')

if emiShdDesired1 == 1 
    mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_setblockprop('Air',0,0,'<None>',0,3,0);
    mi_clearselected;

    if emiShdDesired2 == 1 
        mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
        mi_setblockprop('Air',0,0,'<None>',0,3,0);
        mi_clearselected;
    end
    
end
iterCounter2 = 0;
femTime2 = 0;
maxError = 0;

%%
disp('Iniciando Fase 2 - Geometria do Enrolamento Secundário...');

recRadius = x2;

%##########################################################
%Add Receiver shield:

recShdDesired1 = input('Digite (1) para usar blindagem no receptor ou (0) para não blindar:');

    if recShdDesired1 == 1

        recMatName1 = input('Qual o material da blindagem interna?: ', 's');
        recThicknessShield1 = input('Qual a espessura da blindagem (mm)?:');
        usingShield = 1; %Indicates that ths loosely coupled system has at least one shielding layer.
        
        recXShd = x2 + recCoilHeight/2 +recCoilHeight*0.05;
        
        mi_addnode(0,yfill2+shieldDist);
        mi_addnode(0,yfill2+recThicknessShield1+shieldDist);
        mi_addnode(recXShd,yfill2+shieldDist);
        mi_addnode(recXShd,yfill2+recThicknessShield1+shieldDist);
        mi_addsegment(0,yfill2+shieldDist,recXShd,yfill2+shieldDist);
        mi_addsegment(0,yfill2+recThicknessShield1+shieldDist,recXShd,yfill2+recThicknessShield1+shieldDist);
        mi_addsegment(recXShd,yfill2+recThicknessShield1+shieldDist,recXShd,yfill2+shieldDist);
        mi_addblocklabel(recXShd/2,yfill2+shieldDist+(recThicknessShield1)/2);
        mi_selectlabel(recXShd/2,yfill2+shieldDist+(recThicknessShield1)/2);
        mi_setblockprop(recMatName1,0,0.5,'<None>',0,4,0);
        mi_clearselected;
        mi_selectnode(0,yfill2+shieldDist);
        mi_selectnode(0,yfill2+recThicknessShield1+shieldDist);
        mi_selectnode(recXShd,yfill2+shieldDist);
        mi_selectnode(recXShd,yfill2+recThicknessShield1+shieldDist);
        mi_selectsegment(recXShd/2,yfill2+shieldDist);
        mi_selectsegment(recXShd/2,yfill2+recThicknessShield1+shieldDist);
        mi_selectsegment(recXShd,(yfill2+recThicknessShield1+shieldDist+(yfill2-shieldDist))/2);
        mi_setgroup(4);
        mi_clearselected;

    recShdDesired2 = input('Digite (1) para adicionar segunda blindagem ao receptor ou (0) para não adicionar:');
        if recShdDesired2 == 1
        
            recMatName2 = input('Qual o material da blindagem externa?: ', 's');
            recThicknessShield2 = input('Qual a espessura da segunda blindagem (mm)?:');
        
            mi_selectnode(0,yfill2+recThicknessShield1+shieldDist);
            mi_selectnode(recXShd,yfill2+recThicknessShield1+shieldDist);
            mi_copytranslate(0,recThicknessShield2,1);
            mi_clearselected;
            mi_addblocklabel(recXShd/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
            mi_selectlabel(recXShd/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
            mi_setblockprop(recMatName2,0,0.5,'<None>',0,4,0);
            mi_addsegment(0,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,recXShd,yfill2+recThicknessShield1+shieldDist+recThicknessShield2);
            mi_addsegment(recXShd,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,recXShd,yfill2+recThicknessShield1+shieldDist);
            mi_selectsegment(recXShd/2,yfill2+shieldDist+recThicknessShield1+recThicknessShield2);
            mi_selectsegment(recXShd,(yfill2+recThicknessShield1+recThicknessShield2+shieldDist+(yfill2+shieldDist+recThicknessShield1))/2);
            mi_setgroup(4);
            mi_clearselected;
        end
    
    end

%##########################################################

%Add Receiver "Winglets"

recWing = input('Type "1" to add shielding "winglets" or "0" if it is not necessary:');

if recWing == 1
   %emiWingHeight = input('Enter "winglet" height, in mm:');
   
   if recShdDesired1 == 1
        mi_selectnode(recXShd,yfill2+shieldDist);
        mi_copytranslate(0,-shieldDist-recCoilHeight,1);
        mi_selectnode(recXShd,yfill2-shieldDist-recCoilHeight);
        mi_copytranslate(recThicknessShield1,0,1);
        mi_selectnode(recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight);
        mi_copytranslate(0,shieldDist+recCoilHeight+recThicknessShield1,1);
        mi_clearselected;
        mi_addsegment(recXShd,yfill2+shieldDist,recXShd,yfill2-shieldDist-recCoilHeight);
        mi_addsegment(recXShd,yfill2-shieldDist-recCoilHeight,recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight);
        mi_addsegment(recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight,recXShd+recThicknessShield1,yfill2+shieldDist+recThicknessShield1);
        mi_addsegment(recXShd+recThicknessShield1,yfill2+shieldDist+recThicknessShield1,recXShd,yfill2+shieldDist+recThicknessShield1);
        
        mi_selectnode(recXShd,yfill2-shieldDist-recCoilHeight);
        mi_selectnode(recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight);
        mi_selectnode(recXShd+recThicknessShield1,yfill2+shieldDist+recThicknessShield1);
        mi_selectsegment(recXShd,y2_limite_inf+recCoilHeight/2);
        mi_selectsegment(recXShd+recThicknessShield1/2,y2_limite_inf);
        mi_selectsegment(recXShd+recThicknessShield1,y2_limite_inf+recCoilHeight);
        mi_selectsegment(recXShd+recThicknessShield1/2,yfill2+recThicknessShield1+shieldDist);
        mi_setgroup(44);
        mi_clearselected;
        
        mi_addblocklabel(recXShd+recThicknessShield1/2,yfill2+shieldDist+(recThicknessShield1)/2);
        mi_selectlabel(recXShd+recThicknessShield1/2,yfill2+shieldDist+(recThicknessShield1)/2);
        mi_setblockprop(recMatName1,0,0.5,'<None>',0,44,0);
        mi_clearselected;
        
   end
   
   if recShdDesired2 == 1
        mi_selectnode(recXShd+recThicknessShield1,y2_limite_inf);
        mi_copytranslate(recThicknessShield2,0,1);
        mi_selectnode(recXShd+recThicknessShield1+recThicknessShield2,y2_limite_inf);
        mi_copytranslate(0,recThicknessShield1+recThicknessShield2+shieldDist+recCoilHeight,1);    
        mi_addsegment(recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight,recXShd+recThicknessShield1+recThicknessShield2,y2_limite_inf);
        mi_addsegment(recXShd+recThicknessShield1+recThicknessShield2,yfill2-shieldDist-recCoilHeight,recXShd+recThicknessShield1+recThicknessShield2,yfill2+recThicknessShield1+recThicknessShield2+shieldDist+recCoilHeight);
        mi_addsegment(recXShd+recThicknessShield1+recThicknessShield2,yfill2+recThicknessShield1+recThicknessShield2+shieldDist+recCoilHeight,recXShd,yfill2+recThicknessShield1+recThicknessShield2+shieldDist+recCoilHeight);
        
        mi_selectnode(recXShd+recThicknessShield1+recThicknessShield2,yfill2-shieldDist-recCoilHeight);
        mi_selectnode(recXShd+recThicknessShield1+recThicknessShield2,yfill2+recThicknessShield1+recThicknessShield2+shieldDist);
        
        mi_selectsegment(recXShd+recThicknessShield1+recThicknessShield2/2,y2_limite_inf);
        mi_selectsegment(recXShd+recThicknessShield1+recThicknessShield2,yfill2+recThicknessShield1+recThicknessShield2/2);
        mi_selectsegment(recXShd+recThicknessShield1+recThicknessShield2/2,yfill2+recThicknessShield1+recThicknessShield2);
        mi_setgroup(44);
        mi_clearselected;
        
        mi_addblocklabel(recXShd+recThicknessShield1+recThicknessShield2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_selectlabel(recXShd+recThicknessShield1+recThicknessShield2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_setblockprop(recMatName2,0,0.5,'<None>',0,44,0);
        mi_clearselected;
       
   end
    
end

mi_setcurrent('emitter',0);             %Força uma corrente de teste no circuito receptor do modelo do FEMM.
mi_setcurrent('receiver',recShortCurrent); 
mi_analyze;
mi_loadsolution;
r=mo_getcircuitproperties('receiver');          %Primeira estimativa do erro de emiInduct.
recCalc = abs(r(3)/r(1));   
recCalcError = abs((recInduct - recCalc)/recInduct)*100;
recNewTurn = 0;
recOptim = input('Type "1" to optimize diameter or "0" to optimize height:');

if recOptim == 1
    recMaxHeight = input('Enter maximum emitter height (including shields, if any), in mm:');
    recMaxLayer = floor((recMaxHeight-recThicknessShield1-recThicknessShield2)/recCoilHeight);
    
    if recWing ==1
        mi_selectnode(recXShd,yfill2-shieldDist-recCoilHeight);
        mi_selectnode(recXShd+recThicknessShield1,yfill2-shieldDist-recCoilHeight);
        mi_selectnode(recXShd+recThicknessShield1+recThicknessShield2,yfill2-shieldDist-recCoilHeight);
        mi_movetranslate(0,-recMaxHeight+(recThicknessShield1+recThicknessShield2+recCoilHeight+shieldDist));
        mi_clearselected;
    end
    
end

a = 1;
recShdTop = yfill2+shieldDist;

recClock = tic;
    while (recCalcError < recErrorPrev)                %Enquanto erro atual for menor que erro da iteração anterior...
    
        recErrorPrev = recCalcError;                   % Salva resultado para erro anterior.
        iterCounter2 = iterCounter2 + 1;     
        
        if recOptim == 0 || recLayerCounter == recMaxLayer
            
            if recShdDesired1 == 1                                            
               mi_selectnode(recXShd,recShdTop);
               mi_selectnode(recXShd,recShdTop+recThicknessShield1);
               mi_selectgroup(44);     
               if recShdDesired2 == 1
                  mi_selectnode(recXShd,recShdTop+recThicknessShield1+recThicknessShield2);
               end  
               mi_movetranslate((recCoilHeight+recCoilHeight*0.01),0);
               recXShd = recXShd + recCoilHeight + recCoilHeight*0.01;
            end
            
            mi_selectnode(x2,yfill2);
            mi_selectnode(x2,y2_limite_inf);
            mi_copytranslate((recCoilHeight+recCoilHeight*0.01),0,1);
            mi_selectarcsegment((x2-recCoilHeight/2),yfill2);
            mi_selectarcsegment((x2+recCoilHeight/2),yfill2);
            mi_copytranslate2((recCoilHeight+recCoilHeight*0.01),0,1,3);
            mi_selectlabel(x2,y2_limite_inf);
            mi_copytranslate2((recCoilHeight+recCoilHeight*0.01),0,1,2);
            x2 = x2 + (recCoilHeight+recCoilHeight*0.01);
            recRadius = recRadius + (recCoilHeight+recCoilHeight*0.01);
            recNewTurn = recNewTurn + 1;
                if recLayerCounter == recMaxLayer
                    recLayerCounter = 1;
                    a = -a;
                else
                    recLayerCounter = 0;
                end
        elseif recOptim == 1
            mi_selectnode(x2,yfill2);
            mi_selectnode(x2,y2_limite_inf);
            mi_copytranslate(0,-a*(recCoilHeight+recCoilHeight*0.01),1);
            mi_selectarcsegment((x2-recCoilHeight/2),yfill2);
            mi_selectarcsegment((x2+recCoilHeight/2),yfill2);
            mi_copytranslate2(0,-a*(recCoilHeight+recCoilHeight*0.01),1,3);
            mi_selectlabel(x2,y2_limite_inf);
            mi_copytranslate2(0,-a*(recCoilHeight+recCoilHeight*0.01),1,2);
            yfill2 = yfill2-a*(recCoilHeight+recCoilHeight*0.01);
            y2_limite_inf = y2_limite_inf-a*(recCoilHeight+recCoilHeight*0.01);
            recNewTurn = recNewTurn + 1;
            recLayerCounter = recLayerCounter + 1;
        end
        
        if recRadius > recMaxRadius
            clc;
            maxError = 1;
            disp('Maximum Diameter reached. Receiver inductance will be used "as is".');
            disp('Press any key to continue');
            pause;
            break
        end

        femClock = tic;
        mi_analyze;                                     %Analisa nova geometria e recomeça.
        mi_loadsolution;
        femTime2(iterCounter2) = toc(femClock);
        
        r=mo_getcircuitproperties('receiver');
        recCalc = abs(r(3)/r(1));                                 %Recalcula a indutância emiInduct.
        recCalcError = abs((recInduct - recCalc)/recInduct)*100                                  
        
        recVecError(iterCounter2) = ((recInduct - recCalc)/recInduct)*100;  %Expressa erro percentual para gráfico.
        recVecIter(iterCounter2) = iterCounter2;
        
    end
recTime = toc (recClock);


%Undo Last Iteration
if maxError == 0
    mi_selectnode(x2,yfill2);
    mi_selectnode(x2,y2_limite_inf);
    mi_selectarcsegment((x2-recCoilHeight/2),yfill2);
    mi_selectarcsegment((x2+recCoilHeight/2),yfill2);
    mi_selectlabel(x2,y2_limite_inf);
    mi_deleteselected;

    x2 = x2 - (recCoilHeight+recCoilHeight*0.01);
    recRadius = recRadius - (recCoilHeight+recCoilHeight*0.01);
    recNewTurn = recNewTurn - 1;       
    iterCounter2 = iterCounter2 - 1;
    recVecIter(:,length(recVecIter)) = [];   %Remove last element of this vector.
    recVecError(:,length(recVecError)) = [];
    femTime2(:,length(femTime2)) = [];

    mi_analyze(1);
    mi_loadsolution;
    r = mo_getcircuitproperties('receiver');
    recResistance = real(r(2)/r(1));  
    recCalc = abs(r(3)/r(1));                                 %Recalcula a indutância emiInduct.
    recQFactor = (2*pi*frequency*recCalc)/recResistance;    
    recCalcError = abs((recInduct - recCalc)/recInduct)*100;
end
    
figure;
plot(recVecIter,recVecError);
xlabel('Iteration'); ylabel('Error (%)');
title('Receiver Error vs. Iteration');

figure;
plot(emiVecIter,femTime1,'r');
emiLeg = ('Emitter');
hold on
plot(recVecIter,femTime2,'k');
recLeg = ('Receiver');
xlabel('Iteration'); ylabel('Simulation Time, s');
title('Simulation Time Analysis');
legend(emiLeg,recLeg);
hold off

disp('*************************')
disp('Indutância Secundária desejada (H):');
disp(recInduct);
disp('Indutância obtida (H):');
disp(recCalc);
disp('Número de espiras enrolamento secundário:');
disp(recNewTurn);
disp('Fator de Qualidade do enrolamento secundário:');
disp(recQFactor);
disp('Número de iterações');
disp(iterCounter2);
disp('*************************')

pause;

if emiShdDesired1 == 1 
    mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_setblockprop(emiMatName1,0,0.5,'<None>',0,3,0);
    mi_clearselected;

    if emiShdDesired2 == 1 
        mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
        mi_setblockprop(emiMatName2,0,0.5,'<None>',0,3,0);
        mi_clearselected;
    end
end

emiTurn = emiNewTurn;                                   %Atualiza recTurn e emiTurn com os novos valores determinados.
recTurn = recNewTurn;
    
panc_dfix(1) = y1_limite_sup;
panc_dfix(2) = emiCalc;
panc_dfix(3) = recCalc;
panc_dfix(4) = emiRadius;
panc_dfix(5) = recRadius;
panc_dfix(6) = usingShield;
panc_dfix(7) = emiTurn;
panc_dfix(8) = recTurn;   
    
mi_saveas('result_dimensions.fem'); %Aqui salva arquivo .fem com representação dos indutores.
closefemm

if logDesired == 1
      
    datafile = fopen(fileName,'at');
    fprintf(datafile,'\n \n ****** Resultados para Sistema Espiral com diâmetro fixo ******');
    fprintf(datafile,'\n Resultados para o Emissor: ');
    fprintf(datafile,'\n Número de espiras: %d espiras', emiTurn);
    fprintf(datafile,'\n Indutância obtida: %d H', emiCalc);
    fprintf(datafile,'\n Espessura: %d mm', emiCoilHeight);
    fprintf(datafile,'\n Raio: %d mm', emiRadius);    
    fprintf(datafile,'\n Resistência Série: %d Ohms',emiResistance);
    fprintf(datafile,'\n Fator de Qualidade: %d',emiQFactor);
    fprintf(datafile,'\n Tempo para determinação da geometria: %d segundos',emiTime);
    
    fprintf(datafile,'\n \n Resultados para o Receptor: ');
    fprintf(datafile,'\n Número de espiras: %d espiras', recTurn);
    fprintf(datafile,'\n Indutância obtida: %d H', recCalc);
    fprintf(datafile,'\n Espessura: %d mm', recCoilHeight);
    fprintf(datafile,'\n Raio: %d mm', recRadius);    
    fprintf(datafile,'\n Resistência Série: %d Ohms',recResistance);
    fprintf(datafile,'\n Fator de Qualidade: %d',recQFactor);
    fprintf(datafile,'\n Tempo para determinação da geometria: %d segundos',recTime);
    fclose(datafile);
    
end

delete('temp_dfix.fem');
delete('temp_dfix.ans');
delete('result_wind.fem');
end
