%##########################################################################
%File: panfix.m                                                           #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Axissymmetric representation of distributed coils (pancake) #
%with fixed number of turns.                                              #
%##########################################################################


function panc_nfix = panfix(fileName,logDesired,frequency,emiInduct,recInduct,emiTurn,recTurn,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent)

if nargin <19 
    error('Faltam argumentos de entrada')
end

% clear
% clc
% 
% emiInduct = 120*1e-6;
% recInduct = 46.7*1e-6;
% emiTurn = 15;
% recTurn = 8;
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

testCurrent = 2;
iterCounter1 = 0;
iterCounter2 = 0;
adjFactor = 130;

femTime1 = 0;
nAdj = 0;

targetWind = 0.75;                         %Fator de utilização dos enrolamentos desejado, exemplo 70%.
recWindError = (10/100)*targetWind;
emiWindError = recWindError;
recWindStep = 0;
emiWindStep = 0;

emiError = (0.1/100);                         %Erro aceitável entre emiInduct desejado e emiInduct calculado pelo FEMM, em % de emiInduct.
recError = emiError;                           %Erro aceitável entre recInduct desejado e recInduct calculado pelo FEMM, em % de recInduct.
emiCoilStep = 0;
recCoilStep = 0;

shieldDist = 0.5;                                 %Distância entre a blindagem e o enrolamento
emiShdDesired1 = 0;
emiShdDesired2 = 0;
recShdDesired1 = 0;
recShdDesired2 = 0;
usingShield = 0;

openfemm;                                    %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/result_wind.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('result_wind.fem');
end
mi_saveas('temp_nfix.fem');                    %Aqui salva como um arquivo temporário, para não
                                               %destruir o orignal.

%% Análise Inicial dos Enrolamentos

mi_setcurrent('emitter',testCurrent);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',testCurrent); 
mi_analyze;
mi_loadsolution;
emiWindCalc = mo_getfill(x1,y1);
recWindCalc = mo_getfill(x2,y2);
emiCoilHeight = (yfill1 - y1_limite_inf);
recCoilHeight = (yfill2 - y2_limite_inf);

   recWindCalcError = (targetWind - recWindCalc)/targetWind;
    
    while (abs(recWindCalcError) > recWindError)                %Ajusta espessura do enrolamento sec. até que
                                                         %o fator de utilização seja maior que 70%.
                                                   
        if ((xfill2 - recCoilHeight/2) <= xfill2/10)
            mi_selectgroup(2);
            mi_movetranslate(recCoilHeight,0);
            xfill2 = xfill2 + recCoilHeight;
            x2 = x2 + recCoilHeight;
        end  
    
        mi_selectnode(xfill2,yfill2);
        mi_movetranslate(0,recWindStep)
        yfill2 = yfill2 + recWindStep;
        recCoilHeight = (yfill2 - y2_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        recWindCalc = mo_getfill(x2,y2);
        
        recWindCalcError = (targetWind - recWindCalc)/targetWind;
        recWindStep = -recWindCalcError*0.5;
        
    end

    recCoilHeight = (yfill2 - y2_limite_inf);

    emiWindCalcError = (targetWind - emiWindCalc)/targetWind;
    
    while (abs(emiWindCalcError) > emiWindError)         %Ajusta espessura do enrolamento sec. até que
                                                         %o fator de utilização seja maior que 70%.
                                                   
        if ((xfill1 - emiCoilHeight/2) <= xfill1/10)
            mi_selectgroup(1);
            mi_movetranslate(emiCoilHeight,0);
            xfill1 = xfill1 + emiCoilHeight;
            x1 = x1 + emiCoilHeight;
        end  
    
        mi_selectnode(xfill1,yfill1);
        mi_movetranslate(0,emiWindStep)
        yfill1 = yfill1 + emiWindStep;
        emiCoilHeight = (yfill1 - y1_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        emiWindCalc = mo_getfill(x1,y1);
        
        emiWindCalcError = (targetWind - emiWindCalc)/targetWind;
        emiWindStep = -emiWindCalcError*0.5;                   %O erro é usado como passo de ajuste
        
    end

    emiCoilHeight = (yfill1 - y1_limite_inf);
    
%% Ajuste do Enrolamento Primário

emiIniRadius = input('Enter emitter winding initial radius (mm):');

mi_selectgroup(1);                              %Coloca primário e secundário nos respectivos diâmetros iniciais desejados.
mi_movetranslate(emiIniRadius-x1,0);
x1 = emiIniRadius;
xfill1 = x1;
mi_clearselected;

recIniRadius = input('Enter receiver winding initial radius (mm):');

mi_selectgroup(2);                              %Coloca primário e secundário nos respectivos diâmetros máximos.
mi_movetranslate(recIniRadius-x2,0);
x2 = emiIniRadius;
xfill2 = x2;
mi_clearselected;


mi_selectgroup(1); 
mi_copytranslate(emiCoilHeight+emiCoilHeight*0.01,0,emiTurn-1);          %Aumenta uma espira copiando a espira anterior e deixando entre elas distancia de 0.01Ep.
x1 = x1 + (emiTurn-1)*emiCoilHeight+(emiTurn-1)*emiCoilHeight*0.01;
initialX = x1;

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

mi_setcurrent('emitter',emiOpenCurrent);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',0); 
mi_analyze;
mi_loadsolution;

r=mo_getcircuitproperties('emitter');
emiCalc = abs(r(3)/r(1));                     %Aqui obtém a indutância do primário, calculada, fazendo fluxo/corrente.
emiCalcError = (emiInduct - emiCalc)/emiInduct;

emiClock = tic;
while (abs(emiCalcError) > emiError)     % Procura o raio do enrolamento de emiInduct até que
                                         % a diferença entre o valor desejado e o
                                         % calculado seja inferior ao erro
                                         % permitido.
                                         
     if emiShdDesired1 == 1                                            
        mi_selectnode(x1,y1_limite_inf-shieldDist);
        mi_selectnode(x1,y1_limite_inf-emiThicknessShield1-shieldDist);
    
        if emiShdDesired2 == 1
            mi_selectnode(x1,y1_limite_inf-emiThicknessShield1-shieldDist-emiThicknessShield2);
        end
        
        mi_movetranslate(emiCoilStep,0);
    
    end                                    
                                         
    mi_selectgroup(1);
    mi_movetranslate(emiCoilStep,0)
    x1 = x1 + emiCoilStep;                 %Ajusta o raio do enrolamento primário, de acordo com o passo.
   
    iterCounter1 = iterCounter1 + 1;                %Atualiza o contador de iterações.
    iterCounter2 = iterCounter2 + 1;
    
    femClock = tic;
    mi_analyze;                                     %Analisa nova geometria e recomeça.
    mi_loadsolution;
    femTime1(iterCounter1) = toc(femClock);
    
    r=mo_getcircuitproperties('emitter');
    emiCalc = abs(r(3)/r(1));           %Recalcula a indutância emiInduct.
    emiCalcError = (emiInduct - emiCalc)/emiInduct
    
    emiVecError(iterCounter1) = emiCalcError*100;  %Expressa erro percentual para gráfico.
    emiVecIter(iterCounter1) = iterCounter1;
    
    emiCoilStep = emiCalcError*adjFactor;
    
    if iterCounter2 > 15
        adjFactor = adjFactor/2;
        iterCounter2 = 0;
        
        nAdj = nAdj + 1;
        
        if nAdj > 2
            emiError = 2*emiError; 
            nAdj = 0;
            adjFactor = 100;
        end
        
    end

end

emiTime = toc (emiClock);


mi_loadsolution;
v = mo_getcircuitproperties('emitter');
emiResistance = real(v(2)/v(1));  
emiQFactor = (2*pi*frequency*emiCalc)/emiResistance;

figure;
plot(emiVecIter,emiVecError);
xlabel('Iteration'); ylabel('Error (%)');
title('Emitter Error vs. Iteration');

disp('*************************')
disp('Indutância desejada (H):');
disp(emiInduct);
disp('Indutância obtida (H):');
disp(emiCalc);
emiRadius = x1;
disp('Raio do enrolamento primário(mm):');
disp(emiRadius);
disp('Fator de Qualidade do enrolamento primário:');
disp(emiQFactor);
disp('Número de iterações');
disp(iterCounter1);
disp('*************************')

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

iterCounter3 = 0;
iterCounter4 = 0;
femTime2 = 0;
adjFactor = 100;
nAdj = 0;

%% Ajuste do Enrolamento Secundário
    
mi_selectgroup(2); 
mi_copytranslate(recCoilHeight+recCoilHeight*0.01,0,recTurn-1);              %Aumenta uma espira copiando a espira anterior e deixando entre elas distancia de 0.01Ep.
x2 = x2 + (recTurn-1)*recCoilHeight+(recTurn-1)*recCoilHeight*0.01;

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

mi_setcurrent('emitter',0);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
mi_setcurrent('receiver',recShortCurrent); 
mi_analyze;
mi_loadsolution;

r=mo_getcircuitproperties('receiver');
recCalc = abs(r(3)/r(1));                     %Aqui obtém a indutância do primário, calculada, fazendo fluxo/corrente.

recCalcError = (recInduct - recCalc)/recInduct;

recClock = tic;

while (abs(recCalcError) > recError)     % Procura o raio do enrolamento de recInduct até que
                                   % a diferença entre o valor desejado e o
                                   % calculado seja inferior ao erro
                                   % permitido.
        
    if  recShdDesired1 == 1                                         
        mi_selectnode(x2,yfill2+shieldDist);
        mi_selectnode(x2,yfill2+recThicknessShield1+shieldDist);
        
        if recShdDesired2 == 1
           mi_selectnode(x2,yfill2+recThicknessShield1+shieldDist+recThicknessShield2);
        end
        
        mi_movetranslate(recCoilStep,0);
     end    
    
    mi_selectgroup(2);
    mi_movetranslate(recCoilStep,0)
    x2 = x2 + recCoilStep;
    iterCounter3 = iterCounter3 + 1;
    iterCounter4 = iterCounter4 + 1;
    
    femClock = tic;
    mi_analyze;                                     %Analisa nova geometria e recomeça.
    mi_loadsolution;
    femTime2(iterCounter3) = toc(femClock);
    
    r=mo_getcircuitproperties('receiver');
    recCalc = abs(r(3)/r(1));
    recCalcError = (recInduct - recCalc)/recInduct;

    recVecError(iterCounter3) = recCalcError;
    recVecIter(iterCounter3) = iterCounter3;
    
    recCoilStep = recCalcError*adjFactor;
    
    if iterCounter2 > 15
        adjFactor = adjFactor/2;
        iterCounter2 = 0;
        nAdj = nAdj + 1;
        
        if nAdj > 2
            recError = 2*recError; 
            nAdj = 0;
            adjFactor = 100;
        end
        
    end
    
    disp(recCalc)
      
end
recTime = toc (recClock);

mi_loadsolution;
v = mo_getcircuitproperties('receiver');
recResistance = real(v(2)/v(1));  
recQFactor = (2*pi*frequency*recCalc)/recResistance;

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
disp('Indutância desejada (H):');
disp(recInduct);
disp('Indutância obtida (H):');
disp(recCalc);
recRadius = x2;
disp('Raio do enrolamento secundário(mm):');
disp(recRadius);
disp('Fator de Qualidade do enrolamento secundário:');
disp(recQFactor);
disp('Número de iterações');
disp(iterCounter3);
disp('*************************')

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

%%
y1_limite_sup = yfill1;
panc_nfix(1) = y1_limite_sup;                   %ToDo : verificar quais as saídas desta função.
panc_nfix(2) = emiCalc;
panc_nfix(3) = recCalc;
panc_nfix(4) = emiRadius;
panc_nfix(5) = recRadius;
panc_nfix(6) = usingShield;

mi_saveas('result_dimensions.fem'); %Aqui salva arquivo .fem com representação dos indutores.
closefemm


if logDesired == 1
    
    datafile = fopen(fileName,'at');
    fprintf(datafile,'\n \n ****** Resultados para Sistema Espiral com N fixo ******');
    fprintf(datafile,'\n Resultados para o Emissor: ');
    fprintf(datafile,'\n Número de espiras: %d espiras', emiTurn);
    fprintf(datafile,'\n Indutância obtida: %d H', emiCalc);
    fprintf(datafile,'\n Erro em relação à indutância desejada: %d %',emiCalcError*100);
    fprintf(datafile,'\n Espessura: %d mm', emiCoilHeight);
    fprintf(datafile,'\n Raio: %d mm', emiRadius);    
    fprintf(datafile,'\n Resistência Série: %d Ohms',emiResistance);
    fprintf(datafile,'\n Fator de Qualidade: %d',emiQFactor);
    fprintf(datafile,'\n Tempo para determinação da geometria: %d segundos',emiTime);
    
    fprintf(datafile,'\n \n Resultados para o Receptor: ');
    fprintf(datafile,'\n Número de espiras: %d espiras', recTurn);
    fprintf(datafile,'\n Indutância obtida: %d H', recCalc);
    fprintf(datafile,'\n Erro em relação à indutância desejada: %d %',recCalcError*100);
    fprintf(datafile,'\n Espessura: %d mm', recCoilHeight);
    fprintf(datafile,'\n Raio: %d mm', recRadius);    
    fprintf(datafile,'\n Resistência Série: %d Ohms',recResistance);
    fprintf(datafile,'\n Fator de Qualidade: %d',recQFactor);
    fprintf(datafile,'\n Tempo para determinação da geometria: %d segundos',recTime);
    fclose(datafile);
    
end

delete('temp_nfix.fem');
delete('temp_nfix.ans');
delete('result_wind.fem');
end

