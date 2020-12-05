%##########################################################################
%File: codfix.m                                                           #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Axissymmetric representation of concentrated coils          #
%with maximum diameter restriction.                                       #
%##########################################################################


function conc_dfix = codfix(fileName,logDesired,frequency,emiMaxDiameter,recMaxDiameter,emiInduct,recInduct,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent)

if nargin <19 
    error('Faltam argumentos de entrada')
end

% clear
% clc

% emiInduct = 166.5*1e-6;
% recInduct = 46.7*1e-6;
% emiMaxDiameter = 120;
% recMaxDiameter = 90;
% x1 = 3.3;
% x2 = 3.3;
% y1 = -8.05;
% y2 = 41.95;
% y1_limite_inf = -8.2375;
% y2_limite_inf = 41.7625;
% xfill2 = 3.3;
% yfill2 = 45.7125;
% xfill1 = 3.3;
% yfill1 = -4.5875;

iterCounter1 = 0;                                  %Contador de Iterações 1, inicia com valor nulo.
emiTurn = 1;                                      %Inicia com número de espiras igual a 1.
recTurn = emiTurn;

femTime1 = 0;

emiRadius = emiMaxDiameter/2;
recRadius = recMaxDiameter/2;
emiNewTurn = emiTurn;
recNewTurn = recTurn;

emiErrorPrev = 1e6;     %Erro anterior para o calculo do emissor. Inicia com valor muito alto.
recErrorPrev = 1e6;

emiShdDesired1 = 0;
emiShdDesired2 = 0;
recShdDesired1 = 0;
recShdDesired2 = 0;
shieldDist = 0.5;
usingShield = 0;

openfemm;                                    %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/result_wind.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('result_wind.fem');
end
mi_saveas('temp_dfix.fem');                    %Aqui salva como um arquivo temporário, para não
                                               %destruir o orignal.

%% Fase 1 - Determinação da Geometria do Enrolamento Primário.
disp('Iniciando Fase 1 - Geometria do Enrolamento Primário...');

mi_selectgroup(1);                              %Coloca primário e secundário nos respectivos diâmetros máximos.
mi_movetranslate(emiRadius-x1,0);
x1 = emiRadius;
xfill1 = x1;
mi_clearselected;
initialX = x1;

%###################################

%Add Emitter Shield
emiShdDesired1 = input('Digite (1) para usar blindagem no emissor ou (0) para não blindar:');

if emiShdDesired1 == 1

    emiMatName1 = input('Qual o material da blindagem interna?: ', 's');
    emiThicknessShield1 = input('Qual a espessura da blindagem (mm)?:');
    usingShield = 1; %Indicates that ths loosely coupled system has at least one shielding layer.
    
    mi_addnode(0,y1_limite_inf-shieldDist);
    mi_addnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addnode(x1,y1_limite_inf-shieldDist);
    mi_addnode(x1,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addsegment(0,y1_limite_inf-shieldDist,x1,y1_limite_inf-shieldDist);
    mi_addsegment(0,y1_limite_inf-emiThicknessShield1-shieldDist,x1,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_addsegment(x1,y1_limite_inf-emiThicknessShield1-shieldDist,x1,y1_limite_inf-shieldDist);
    mi_addblocklabel(x1/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_selectlabel(x1/2,((y1_limite_inf-emiThicknessShield1-shieldDist)+(y1_limite_inf-shieldDist))/2);
    mi_setblockprop(emiMatName1,0,0.5,'<None>',0,3,0);
    mi_clearselected;
    mi_selectnode(0,y1_limite_inf-shieldDist);
    mi_selectnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectnode(x1,y1_limite_inf-shieldDist);
    mi_selectnode(x1,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectsegment(x1/2,y1_limite_inf-shieldDist);
    mi_selectsegment(x1/2,y1_limite_inf-emiThicknessShield1-shieldDist);
    mi_selectsegment(x1,(y1_limite_inf-emiThicknessShield1-shieldDist+(y1_limite_inf-shieldDist))/2);
    mi_setgroup(3);
    mi_clearselected;

        emiShdDesired2 = input('Digite (1) para adicionar segunda blindagem ao emissor ou (0) para não adicionar:');
        if emiShdDesired2 == 1
            emiMatName2 = input('Qual o material da blindagem externa?: ', 's');
            emiThicknessShield2 = input('Qual a espessura da segunda blindagem (mm)?:');
        
            mi_selectnode(0,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_selectnode(x1,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_copytranslate(0,-emiThicknessShield2,1);
            mi_clearselected;
            mi_addblocklabel(x1/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_selectlabel(x1/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_setblockprop(emiMatName2,0,0.5,'<None>',0,3,0);
            mi_addsegment(0,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2,x1,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2);
            mi_addsegment(x1,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2,x1,y1_limite_inf-emiThicknessShield1-shieldDist);
            mi_selectsegment(x1/2,y1_limite_inf-shieldDist-emiThicknessShield1-emiThicknessShield2);
            mi_selectsegment(x1,(y1_limite_inf-emiThicknessShield1-emiThicknessShield2-shieldDist+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
            mi_setgroup(3);
            mi_clearselected;
        end
end


%###################################

mi_selectlabel(x1,y1);                                       %Os comandos iniciados por "mi_" referem-se à discretização de elementos magnéticos (m)
mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,emiTurn);      %na etapa inicial (i), ou seja, no pré-processamento.
mi_clearselected;

mi_selectlabel(x2,y2);
mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recTurn);
mi_clearselected;
mi_setcurrent('emitter',emiOpenCurrent);             %Força uma corrente no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',0); 
mi_analyze;
mi_loadsolution;

r=mo_getcircuitproperties('emitter');          %Primeira estimativa do erro de emiInduct.
emiCalc = abs(r(3)/r(1));   
emiCalcError = abs((emiInduct - emiCalc)/emiInduct)*100;

emiClock = tic;
while (emiCalcError < emiErrorPrev)                %Enquanto erro atual for menor que erro da iteração anterior...
    
        emiErrorPrev = emiCalcError;                   % Salva resultado para erro anterior.
        iterCounter1 = iterCounter1 + 1;
        
        
        if (emiCalcError > 0)                      
            emiNewTurn = emiNewTurn + 1;          % Se o erro é positivo, emiCalc < Lp_desejado. Então aumenta uma espira.     
        elseif (emiCalcError < 0)
            emiNewTurn = emiNewTurn - 1;          % Se o erro é negativo, emiCalc > Lp_desejado. Então diminui uma espira. 
        end

        mi_setgroup(1);
        mi_selectlabel(emiRadius,y1);
        mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,emiNewTurn);
        mi_clearselected;
    
        femClock = tic;
        mi_analyze;                                     %Analisa nova geometria e recomeça.
        mi_loadsolution;
        femTime1(iterCounter1) = toc(femClock);
        
        
        emiWindFactor = mo_getfill(x1,y1);

            while ((emiWindFactor < 0.7) || (emiWindFactor > 0.8))                 %Ajusta espessura do enrolamento prim. até que
                                                                                   %o fator de utilização seja maior que 70%.                                   
                if emiWindFactor > 0.8
                    windStep = 0.2;
                elseif emiWindFactor < 0.7
                    windStep = -0.15;
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
        
    r=mo_getcircuitproperties('emitter');
    emiCalc = abs(r(3)/r(1));   %Recalcula a indutância emiInduct.
    emiCalcError = abs((emiInduct - emiCalc)/emiInduct)*100;                                  
    
    emiVecError(iterCounter1) =((emiInduct - emiCalc)/emiInduct)*100;  %Expressa erro percentual para gráfico.
    emiVecIter(iterCounter1) = iterCounter1;
    
end

emiTime = toc (emiClock); %Total time for emitter calculation.

%Repete a iteração de menor erro:
emiNewTurn = emiNewTurn - 1; %Reduz uma espira, pois o menor erro aconteceu para a iteração anterior ao fim do while.
iterCounter1 = iterCounter1 - 1;
emiVecIter(:,length(emiVecIter)) = [];   %Remove last element of this vector.
emiVecError(:,length(emiVecError)) = [];
femTime1(:,length(femTime1)) = [];

mi_setgroup(1);
mi_selectlabel(emiRadius,y1);
mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,1,emiNewTurn);   
mi_clearselected;
mi_analyze;                                     %Analisa nova geometria e recomeça.
mi_loadsolution;
    
    emiWindFactor = mo_getfill(x1,y1);

            while ((emiWindFactor < 0.7) || (emiWindFactor > 0.8))                 %Ajusta espessura do enrolamento prim. até que
                                                                                   %o fator de utilização seja maior que 70%.                                   
                if emiWindFactor > 0.8
                    windStep = 0.2;
                elseif emiWindFactor < 0.7
                    windStep = -0.15;
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
    y1_limite_sup = yfill1;
        
    r=mo_getcircuitproperties('emitter');
    emiCalc = abs(r(3)/r(1));   %Recalcula a indutância emiInduct.
    
r = mo_getcircuitproperties('emitter');
emiResistance = real(r(2)/r(1));  
emiQFactor = (2*pi*frequency*emiCalc)/emiResistance;

figure;
plot(emiVecIter,emiVecError);
xlabel('Iteration'); ylabel('Error (%)');
title('Emitter Error vs. Iteration');

disp('*************************')
disp('Indutância Primária desejada (H):');
disp(emiInduct);
disp('Indutância obtida (H):');
disp(emiCalc);
emiRadius = x1;
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
end

if emiShdDesired2 == 1 
    mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
    mi_setblockprop('Air',0,0,'<None>',0,3,0);
    mi_clearselected;
end

iterCounter2 = 0; 
femTime2 = 0;

%% Fase 2 - Determinação da Geometria do Enrolamento Secundário.
disp('Iniciando Fase 2 - Geometria do Enrolamento Secundário...');

mi_selectgroup(2);
mi_movetranslate(recRadius-x2,0);
x2 = recRadius;
xfill2 = x2;
mi_clearselected;

%##########################################################
%Add Receiver shield:

recShdDesired1 = input('Digite (1) para usar blindagem no receptor ou (0) para não blindar:');

    if recShdDesired1 == 1

        recMatName1 = input('Qual o material da blindagem interna?: ', 's');
        recThicknessShield1 = input('Qual a espessura da blindagem (mm)?:');
        usingShield = 1; %Indicates that ths loosely coupled system has at least one shielding layer.
        
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
        
        recMatName2 = input('Qual o material da blindagem externa?: ', 's');
        recThicknessShield2 = input('Qual a espessura da segunda blindagem (mm)?:');
        
        mi_selectnode(0,yfill2+recThicknessShield1+shieldDist);
        mi_selectnode(x2,yfill2+recThicknessShield1+shieldDist);
        mi_copytranslate(0,recThicknessShield2,1);
        mi_clearselected;
        mi_addblocklabel(x2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_selectlabel(x2/2,((yfill2+recThicknessShield1+shieldDist+recThicknessShield2)+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_setblockprop(recMatName2,0,0.5,'<None>',0,4,0);
        mi_addsegment(0,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,x2,yfill2+recThicknessShield1+shieldDist+recThicknessShield2);
        mi_addsegment(x2,yfill2+recThicknessShield1+shieldDist+recThicknessShield2,x2,yfill2+recThicknessShield1+shieldDist);
        mi_selectsegment(x2/2,yfill2+shieldDist+recThicknessShield1+recThicknessShield2);
        mi_selectsegment(x2,(yfill2+recThicknessShield1+recThicknessShield2+shieldDist+(yfill2+shieldDist+recThicknessShield1))/2);
        mi_setgroup(4);
        mi_clearselected;
    end
    
    end

%##########################################################

mi_selectlabel(x2,y2);
mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recTurn);
mi_clearselected;
mi_setcurrent('emitter',0);             %Força uma corrente no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',recShortCurrent); 
mi_analyze;
mi_loadsolution;


r=mo_getcircuitproperties('receiver');          %Primeira estimativa do erro de emiInduct.
recCalc = abs(r(3)/r(1));   
recCalcError = abs((recInduct - recCalc)/recInduct)*100;

recClock = tic;
    while (recCalcError < recErrorPrev)
    
        recErrorPrev = recCalcError;
        iterCounter2 = iterCounter2 + 1;
   
        if (recCalcError > 0)                      
            recNewTurn = recNewTurn + 1;          % Se o erro é positivo, emiCalc < Lp_desejado. Então aumenta uma espira.     
        elseif (recCalcError < 0)
            recNewTurn = recNewTurn - 1;          % Se o erro é negativo, emiCalc > Lp_desejado. Então diminui uma espira. 
        end

    mi_setgroup(2);
    mi_selectlabel(recRadius,y2);
    mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recNewTurn);
    mi_clearselected;
    
    femClock = tic;
    mi_analyze;                                     %Analisa nova geometria e recomeça.
    mi_loadsolution;
    femTime2(iterCounter2) = toc(femClock);
    
    recWindFactor = mo_getfill(x2,y2);

            while ((recWindFactor < 0.7) || (recWindFactor > 0.8))                 %Ajusta espessura do enrolamento prim. até que
                                                       %o fator de utilização seja maior que 70%.                                   
                if recWindFactor > 0.8
                    windStep = 0.2;
                elseif recWindFactor < 0.7
                    windStep = -0.15;
                end                          
        
                %###################
                mi_selectgroup(4);
                mi_movetranslate(0,windStep);
                %##################
                
                mi_selectnode(xfill2,yfill2);
                mi_movetranslate(0,windStep)
                yfill2 = yfill2 + windStep;
                recCoilHeight = (yfill2 - y2_limite_inf);
         
                mi_analyze;
                mi_loadsolution;
                recWindFactor = mo_getfill(x2,y2); 
            end

    recCoilHeight = (yfill2 - y2_limite_inf);
    y2_limite_sup = yfill2;
        
    r=mo_getcircuitproperties('receiver');
    recCalc = abs(r(3)/r(1));   
    recCalcError = abs((recInduct - recCalc)/recInduct)*100;                                   

    recVecError(iterCounter2) =((recInduct - recCalc)/recInduct)*100;  %Expressa erro percentual para gráfico.
    recVecIter(iterCounter2) = iterCounter2;

    end
    
recTime = toc (recClock);

%Undo Last Iteration (Minimum Error)
recNewTurn = recNewTurn - 1;
iterCounter2 = iterCounter2 - 1;
recVecIter(:,length(recVecIter)) = [];   %Remove last element of this vector.
recVecError(:,length(recVecError)) = [];
femTime2(:,length(femTime2)) = [];

    mi_setgroup(2);
    mi_selectlabel(recRadius,y2);
    mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,2,recNewTurn);
    mi_clearselected;
    mi_analyze;                                    
    mi_loadsolution;
    
    recWindFactor = mo_getfill(x2,y2);

            while ((recWindFactor < 0.7) || (recWindFactor > 0.8))                 %Ajusta espessura do enrolamento prim. até que
                                                       %o fator de utilização seja maior que 70%.                                   
                if recWindFactor > 0.8
                    windStep = 0.2;
                elseif recWindFactor < 0.7
                    windStep = -0.15;
                end                          
        
                %###################
                mi_selectgroup(4);
                mi_movetranslate(0,windStep);
                %##################
                
                mi_selectnode(xfill2,yfill2);
                mi_movetranslate(0,windStep)
                yfill2 = yfill2 + windStep;
                recCoilHeight = (yfill2 - y2_limite_inf);
         
                mi_analyze;
                mi_loadsolution;
                recWindFactor = mo_getfill(x2,y2); 
            end

    recCoilHeight = (yfill2 - y2_limite_inf);
    y2_limite_sup = yfill2;
        
    r=mo_getcircuitproperties('receiver');
    recCalc = abs(r(3)/r(1)); 
    

recResistance = real(r(2)/r(1));
recQFactor = (2*pi*frequency*emiCalc)/recResistance;

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
recRadius = x2;
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
end

if emiShdDesired2 == 1 
    mi_selectlabel(initialX/2,((y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2)+(y1_limite_inf-shieldDist-emiThicknessShield1))/2);
    mi_setblockprop(emiMatName2,0,0.5,'<None>',0,3,0);
    mi_clearselected;
end

emiTurn = emiNewTurn;                                   %Atualiza recTurn e emiTurn com os novos valores determinados.
recTurn = recNewTurn;

conc_dfix(1) = y1_limite_sup;                   %ToDo : verificar quais as saídas desta função.
conc_dfix(2) = emiCalc;
conc_dfix(3) = recCalc;
conc_dfix(4) = emiRadius;
conc_dfix(5) = recRadius;
conc_dfix(6) = emiTurn;
conc_dfix(7) = recTurn;
conc_dfix(8) = usingShield;

mi_saveas('result_dimensions.fem'); %Aqui salva arquivo .fem com representação dos indutores.
closefemm

if logDesired == 1
      
    datafile = fopen(fileName,'at');
    fprintf(datafile,'\n \n ****** Resultados para Sistema Concentrado com diâmetro fixo ******');
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

