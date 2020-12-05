%##########################################################################
%File: misalign.m                                                         #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: December 7, 2013                                            #
%Description: Misalignment analysis routine                               #
%##########################################################################

function misAl = misalign (frequency,emiCalc,emiOpenCurrent,yfill2,y2_limite_inf,dist,emiRadius)

flagFirst = 0;
recRadius = emiRadius;

openfemm;
try
    vv=ver;
    opendocument([cd,'/result_planar.fem']); %Recupera geometria representada na forma planar.
catch
    opendocument('result_planar.fem');     
end

mi_saveas('temp_align.fem');

disp('Iniciando Estudo de Enrolamentos Não-Coaxiais...');
mi_setcurrent('receiver',0);                  %Força corrente nula no circuito secundário do modelo do FEMM.
mi_setcurrent('emitter',emiOpenCurrent);      %Força corrente com secundário em aberto no modelo do FEMM. 

deltaDc = 10;                 %Passo de dc para realizar as simulações. Quanto maior, mais detalhado o estudo.
plotCounter = 0;               %Contador de quantas curvas (k x dc) foram traçadas
% plotNumber = floor(dist/10)-1; %O número de curvas a serem plotas é obtido de acordo com a distância entre
%                                %primário e secundário.
%                                
% if plotNumber <= 0
%    plotNumber = 1;
% end
% 
% if plotNumber > 5;
%     plotNumber = 5;
% end

plotNumber = 2;   %LIXO -- RETIRAR.
  
while plotCounter < plotNumber
  counter = 0;                                    %Reseta os contadores para reiniciar.
  counter2 = 0;  
  mi_analyze(1);
  mi_loadsolution;
    
for counter2 = 0:emiRadius/deltaDc:4*emiRadius                    %Move o secundário da esquerda para a direita.
  counter = counter + 1;
 
  h=mo_getcircuitproperties('receiver'); %Aqui calcula M e k.
  emiOpenCalc2 = abs(h(2));
  u=mo_getcircuitproperties('emitter');
  mutualCalc = abs(emiOpenCalc2/(2*pi*frequency*u(1)));
  couplingCalc(counter) = mutualCalc/emiCalc;      %sqrt(emiInduct*recInduct) = emiInduct, pois emiInduct = recInduct.
  dc(counter) = counter2; 
    
  mi_selectgroup(2);                      %Depois de calcular k para a posição espacial atual,
  mi_selectgroup(4);                      % move o secundário lateralmente e reinicia.
  mi_movetranslate(emiRadius/deltaDc,0);
  mi_clearselected;
  mi_analyze(1);
  mi_loadsolution;
    
end

if flagFirst == 0
     figure(3);
     n = length(dc);
     [rho,theta] = meshgrid(dc,linspace(0,2*pi,n));
     Z = repmat(couplingCalc(:).',n,1);
     X = rho.*cos(theta);
     Y = rho.*sin(theta);
     surf(X,Y,Z)
     colormap jet;
     colorbar;
     
     flagFirst = 1;
end
figure(4);
plot(dc,couplingCalc);                     %Plota a curva (k x dc) para uma distância e.

hold on

  mi_selectgroup(2);                                %Depois que plotou a curva para uma distância e, coloca o enrolamento
  mi_selectgroup(4);                                %secundário na posição inicial (bem à esquerda do primário), para ser usado na
  mi_movetranslate(-counter*emiRadius/deltaDc,0);   %próxima iteração.
  mi_clearselected;
   
  mi_selectgroup(2);                        %Aqui determina qual a nova distância e. Apenas se divide a distância primário-secundário
  mi_selectgroup(4);                        %(que foi obtida na representação por eixo de simetria) em partes iguais e se realiza a varredura
  mi_movetranslate(0,-(dist)/(plotNumber)); %de dc para cada uma delas.
  mi_clearselected;

plotCounter = plotCounter + 1;            %Atualiza o contador de quantas curvas foram traçadas.

end
disp('Fim Estudo de Enrolamentos Não-Coaxiais. Pressione qualquer tecla para continuar...');

delete('temp_align.fem');
delete('temp_align.ans');

pause

closefemm;


%% Estudo da Variação de alpha (enrolamentos não-paralelos)

disp('Iniciando Estudo de Enrolamentos Não-Paralelos...');

openfemm;
try
    vv=ver;
    opendocument([cd,'/result_planar.fem']); %Recupera geometria representada na forma planar.
catch
    opendocument('result_planar.fem');     
end

mi_saveas('temp_align2.fem');

mi_setcurrent('receiver',0);      %Força corrente nula no circuito secundário do modelo do FEMM.
mi_setcurrent('emitter',emiOpenCurrent);      %Força corrente com secundário em aberto no modelo do FEMM. 

mi_selectgroup(4);
mi_setgroup(2);
mi_clearselected;

rotationCenter = y2_limite_inf + (yfill2 - y2_limite_inf)/2;  %Adiciona um ponto no centro do enrolamento secundário esquerdo para ser o
                                                         %centro de giro.
mi_addnode(-recRadius,rotationCenter);
mi_selectnode(-recRadius,rotationCenter);
mi_setgroup(2);
mi_clearselected;

%Girando o secundário no sentido anti-horário
alpha = 0;       %Angulo inicial nulo, enrolamentos paralelos.
deltaAlpha = 5; %Incremento de alpha.
couplingCalc = 0;
dc = 0;

while dc <= 2*recRadius
  counter = 0;
  counter2 = 0;
  mi_analyze(1);
  mi_loadsolution;

  for counter2 = 0:deltaAlpha:(180-deltaAlpha)

  counter = counter+1;
  b=mo_getcircuitproperties('receiver'); %Aqui calcula M e k.
  emiOpenCalc = abs(b(2));
  w=mo_getcircuitproperties('emitter');
  mutualCalc = abs(emiOpenCalc/(2*pi*frequency*w(1)));
  couplingCalc(counter) = mutualCalc/emiCalc;      %sqrt(emiInduct*recInduct) = emiInduct, pois emiInduct = recInduct.
  
  alpha(counter) = counter2;
  
  mi_selectgroup(2);
  mi_moverotate(-recRadius+dc,rotationCenter,deltaAlpha)  %Gira o secundário em passos de 5 graus.
  mi_clearselected;                            %O centro de giro é o ponto central do lado esquerdo do enrolamento secundário.

  mi_analyze(1);
  mi_loadsolution;

  end

plot(alpha,couplingCalc);
hold on

mi_selectgroup(2);                                %Retorna o secundário para a posição inicial.
mi_moverotate(-recRadius+dc,rotationCenter,-counter*deltaAlpha);
mi_clearselected; 

dc = dc + recRadius/5;                          %Move o secundário para a direita e reinicia.
mi_selectgroup(2);
mi_movetranslate(recRadius/5,0);
mi_clearselected;

end

misAl(1) = 1;

delete('temp_align2.fem');
delete('temp_align2.ans');

end