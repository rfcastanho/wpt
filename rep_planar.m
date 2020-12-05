%##########################################################################
%File: rep_planar.m                                                       #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Iterative coil representation (radius, number of turns,     # 
%inductance) in planar approach. Provides lateral and angular misalignment#
%analysis.                                                                #
%##########################################################################


function planar = rep_planar(fileName,logDesired,frequency,y1,y2_limite_inf,emiRadius,emiTurn,emiInduct,dist)

if nargin <7
    error('Faltam argumentos de entrada')
end

% clear
% clc
% frequency = 60e3;
% y1 =-8.05;
% y2_limite_inf =114.7625;
% emiRadius = 134.7809;
% emiturn = 15;
% emiInduct = 166.5e-6;
% emiOpenCurrent = 10;
% dist = 118.9017;
% yfill2 =116.8455;

recTurn = emiTurn;
recRadius = emiRadius;   
testCurrent = 2;            %Corrente de teste dos enrolamentos (A);

emiNewTurn = emiTurn;

openfemm;
try
    vv=ver;
    opendocument([cd,'/result_coupled.fem']); %Abre geometria representada por simetria.
catch
    opendocument('result_coupled.fem');       %Abre geometria representada por simetria.
end
mi_saveas('temp_planar.fem');                   %Salva como temporário para manipular.

mi_probdef(frequency,'millimeters','planar',1e-8,2*emiRadius,20,(0)); %Configurações iniciais do problema. 
                                                                      %Profundidade igual a 2*emiRadius = Dp
%% Ajuste da Janela de Simulação                                                       
mi_selectnode(0,300);
mi_selectnode(0,-300);
mi_movetranslate(-4*emiRadius-650,0);
mi_clearselected;
mi_selectnode(650,300);
mi_selectnode(650,-300);
mi_movetranslate(4*emiRadius,0);
mi_clearselected;
mi_selectnode(2*emiRadius+650,300);
mi_selectnode(-2*emiRadius-650,300);
mi_movetranslate(0,3*dist);
mi_clearselected;

mi_selectgroup(1);
mi_copytranslate(-2*emiRadius,0,1);
mi_clearselected;

                                           %Aqui Configura o enrolamento primário, separando-o 
mi_selectlabel(-emiRadius,y1);                    %em dois grupos (grupo 3 com as espiras que "entram na página"
mi_selectnode(-emiRadius,y1);                     %e grupo 1 com as espiras que saem da página)
mi_selectnode(-emiRadius,y1+10);
mi_selectarcsegment(-emiRadius-1,y1);
mi_selectarcsegment(-emiRadius+1,y1);
mi_setgroup(3);
mi_selectlabel(-emiRadius,y1);
mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,3,-emiTurn);
mi_clearselected;

mi_selectgroup(2);
mi_copytranslate(-2*recRadius,0,1);
mi_clearselected;

                                         %Aqui faz o mesmo para os grupos 2 e 4, que representam o enrolamento
mi_selectlabel(-recRadius,y2_limite_inf);       %secundário.
mi_selectlabel(-recRadius,y2_limite_inf);
mi_selectnode(-recRadius,y2_limite_inf);
mi_selectnode(-recRadius,y2_limite_inf+10);
mi_selectarcsegment(-recRadius-1,y2_limite_inf);
mi_selectarcsegment(-recRadius+1,y2_limite_inf);
mi_setgroup(4);
mi_selectlabel(-recRadius,y2_limite_inf);
mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,4,-recTurn);
mi_clearselected;

%% Cálculo das Indutâncias próprias (emiInduct = recInduct) pela representação planar

mi_setcurrent('receiver',0);      %Força corrente nula no circuito secundário do modelo do FEMM.
mi_setcurrent('emitter',testCurrent);   %Força corrente de teste no circuito primário do modelo do FEMM.
mi_analyze;
mi_loadsolution;

h=mo_getcircuitproperties('emitter');
emiCalc = abs(h(3))/h(1);           %Aqui obtém emiInduct calculado pela representação planar.

emiError = ((emiInduct - emiCalc)/emiInduct)*100;  %Calcula erro entre a indutância desejada e a obtida.

while (abs(emiError) > 10)           %Enquanto erro for maior que 10%, ajusta número de espiras.
                                    %Obs1: O número de espiras é sempre um número inteiro.
                                    %Obs2: A tolerância de 10% é suficiente para análise do desempenho dos enrolamentos.
                                    
if (emiError > 0)                    
    
    emiNewTurn = emiNewTurn + 1;          % Se o erro é positivo, emiCalc < Lp_desejado. Então aumenta uma espira. 
    
elseif (emiError < 0)
    
    emiNewTurn = emiNewTurn - 1;          % Se o erro é negativo, emiCalc > Lp_desejado. Então diminui uma espira. 
end

mi_setgroup(3);
mi_selectlabel(-emiRadius,y1);
mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,3,-emiNewTurn);  %Atualiza o enrolamento primário com o novo número
mi_clearselected;                                     % de espiras

mi_setgroup(1);
mi_selectlabel(emiRadius,y1);
mi_setblockprop('Litz Emissor',0,0.5,'emitter',0,3,emiNewTurn);
mi_clearselected;

mi_analyze;
mi_loadsolution;

h=mo_getcircuitproperties('emitter');
emiCalc = abs(h(3))/h(1);

emiError = ((emiInduct - emiCalc)/emiInduct)*100;                   %Recalcula erro para reiniciar ajuste de espiras ou sair do loop.

end

recNewTurn = emiNewTurn;

%% Atualização do Enrolamento Secundário

mi_setgroup(4);                               %Aqui atualiza o secundário com os mesmos resultados obtidos para o primário,
mi_selectlabel(-recRadius,y2_limite_inf);            %já que são iguais.
mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,4,-recNewTurn);
mi_clearselected;

mi_setgroup(2);
mi_selectlabel(recRadius,y2_limite_inf);
mi_setblockprop('Litz Receptor',0,0.5,'receiver',0,4,recNewTurn);
mi_clearselected;

mi_saveas('result_planar.fem');     %Aqui salva arquivo .fem com geometria final.
closefemm;


if logDesired == 1
    
    datafile = fopen(fileName,'at');
    fprintf(datafile,'\n \n ****** Resultados para Representação Planar ******');
    fprintf(datafile,'\n Novo Número de espiras: %d espiras', emiNewTurn);
    fprintf(datafile,'\n Número de espiras iniciais: %d espiras', emiTurn);
    fprintf(datafile,'\n Indutância obtida na representação planar: %d H', emiCalc);
    fprintf(datafile,'\n Erro em relação à indutância obtida por eixo de simetria: %d %%', abs(emiError));
    fclose(datafile);
    
end

delete('temp_planar.fem');
delete('temp_planar.ans');

proceedMisalign = 1;                 %If everything is ok, set flag to indicate Misalignment study is possible. 

planar(1) = proceedMisalign;
planar(2) = emiCalc;

end
