%##########################################################################
%File: coupling.m                                                         #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Adjust distance between emitter and receiver until correct  #
%mutual coupling is achieved (based on user inputs in mag_eqv.m).         #
%##########################################################################

function mutual = coupling(fileName,logDesired,recOpenVoltage,emiOpenCurrent,frequency,y2_limite_inf,y1_limite_sup,yfill2,emiCalc,recCalc,usingShield)

if nargin <11
    error('Faltam argumentos de entrada')
end

% clear
% clc
% recOpenVoltage = 65;
% emiOpenCurrent = 10;
% frequency = 60e3;
% y2_limite_inf = 41.7625;
% y1_limite_sup = -4.1392;
% yfill2 = 43.8455;
% emiCalc = 1.6644e-4;
% recCalc = 1.6644e-4;

adjFactor =     100;
nAdj =            0;
errorFlag =       0;
iterCounter1 =    0;
mutualStep =      0;
emiOpenError = 1/100;                              %Erro aceitável entre recOpenVoltage desejado e recOpenVoltage calculado pelo FEMM, em V.
testCurrent = 2;

openfemm;                                         %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/result_dimensions.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('result_dimensions.fem');
end
mi_saveas('temp_coupling.fem');               %Aqui salva como um arquivo temporário, para não
                                              %destruir o modelo.                                            
                                              
%% Posicionamento do Secundário para Obter o Fator de Acoplamento Correto
disp('Iniciando Fase 3 - Obtenção da Indutância Mútua...');

mi_setcurrent('receiver',0);               %Força corrente nula no circuito secundário do modelo do FEMM.
mi_setcurrent('emitter',emiOpenCurrent);    %Força corrente em circuito aberto no circuito primário do modelo do FEMM (depende do circuito simulado).
mi_analyze;
mi_loadsolution;

r=mo_getcircuitproperties('receiver'); %Obtém os resultados do circuito secundário, para avaliar a tensão induzida em circuito aberto.
emiOpenCalc = abs(r(2));
emiOpenCalcError = ((recOpenVoltage - emiOpenCalc)/recOpenVoltage);

while (abs(emiOpenCalcError) > emiOpenError)
 
    mi_selectgroup(2);
    mi_selectgroup(4);
    mi_movetranslate(0,mutualStep)
    y2_limite_inf  = y2_limite_inf  + mutualStep;    %Atualiza coordenadas do enrolamento secundário
    yfill2 = yfill2 + mutualStep;
    
    if y2_limite_inf < y1_limite_sup   %Se o secundário é aproximado demais do primário (distância e = 0),
                                       %indica erro, porque esta configuração de Ip, recOpenVoltage, Dp e Ds não existe. 
        disp('Erro: Esta geometria não é possível.')
        errorFlag = 1;
        break
    end  
      
    mi_analyze(1);
    mi_loadsolution;
    r=mo_getcircuitproperties('receiver');
    emiOpenCalc = abs(r(2));
    
    iterCounter1 = iterCounter1 + 1;
    if iterCounter1 > 15
        iterCounter1 = 0;
        adjFactor = adjFactor/2;
        nAdj = nAdj + 1;
        
        if nAdj > 1
           nAdj = 0;
           emiOpenError = emiOpenError*2;
           adjFactor = 100;
        end
    end
    
    emiOpenCalcError = ((recOpenVoltage - emiOpenCalc)/recOpenVoltage);
    mutualStep = -emiOpenCalcError*adjFactor;
    disp(emiOpenCalcError)
    
end

if errorFlag == 0      %Se não tem erro de geometria, mostra resultados.
    
dist = y2_limite_inf - y1_limite_sup;

u=mo_getcircuitproperties('emitter');
mutualCalc = abs(emiOpenCalc/(2*pi*frequency*u(1)));

    if usingShield == 1
        mi_setcurrent('receiver',testCurrent);  %If any kind of magnetic shield is used, self inductances must be calculated before
        mi_setcurrent('emitter',0);             % solving k = M/sqrt(emiCalc*recCalc).
        mi_analyze;
        mi_loadsolution;
        rr=mo_getcircuitproperties('receiver');
        recCalc = abs(rr(3)/rr(1));
    
        mi_setcurrent('emitter',emiOpenCurrent); 
        mi_setcurrent('receiver',0);
        mi_analyze;
        mi_loadsolution;
        r=mo_getcircuitproperties('emitter');
        emiCalc = abs(r(3)/r(1));
    end


couplingCalc = mutualCalc/(sqrt(emiCalc*recCalc));

proceedMutual = 1; %Aqui finalizou a representação por eixo de simetria. Pode prosseguir para a rep. planar, se for o caso.
mutual(1) = proceedMutual;
mutual(2) = yfill2;
mutual(3) = dist;
mutual(4) = y2_limite_inf;

mi_saveas('result_coupled.fem');        %Aqui salva arquivo .fem com geometria final.
end

if logDesired == 1
    
    datafile = fopen(fileName,'at');
    fprintf(datafile,'\n \n ****** Resultados do Acoplamento Magnético ******');
    fprintf(datafile,'\n Distância entre enrolamentos: %d mm', dist);
    fprintf(datafile,'\n Indutância Mútua: %d H', mutualCalc);
    fprintf(datafile,'\n Coeficiente de Acoplamento: %d', couplingCalc);
    fclose(datafile);
    
end

closefemm                          %Fecha FEMM.

delete('temp_coupling.fem');
delete('temp_coupling.ans');

end