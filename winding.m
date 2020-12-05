%##########################################################################
%File: winding.m                                                          #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Iterative winding representation (fill factor, height and   #
%conductors).                                                             #
%##########################################################################

function wind = winding(emiOpenCurrent,recShortCurrent,frequency)
global prt

if nargin <3
    error('Faltam argumentos de entrada')
end

% x1 = 3.3;                                    %Posições Iniciais dos dos enrolamentos.
% y1 = -8.05;                                  %Estas coordenadas devem-se ao arquivo inicial desenvolvido
% x2 = 3.3;                                    %no FEMM.
% y2 = 41.95;
% y1_limite_inf = -8.2375;
% y2_limite_inf = 41.7625;
% xfill2 = 3.3;
% yfill2 = 45.8625;
% xfill1 = 3.3;
% yfill1 = -4.1375;

x1 = 5;                                    %Posições Iniciais dos dos enrolamentos.
y1 = -19.95;                                  %Estas coordenadas devem-se ao arquivo inicial desenvolvido
x2 = 5;                                    %no FEMM.
y2 = 40.05;
y1_limite_inf = -20;
y2_limite_inf = 40;
xfill2 = x2;
yfill2 = 45;
xfill1 = x1;
yfill1 = -15;


targetWind = 0.75;                         %Fator de utilização dos enrolamentos desejado, exemplo 75%.
recWindError = (10/100)*targetWind;
emiWindError = recWindError;
recWindStep = 0;
emiWindStep = 0;
mesh = 1;
testCurrent = 2;
naux = 2;                                   %Número de espiras para primeira avaliação do enrolamento.
Jp = 3.5;                                                     %Densidade de Corrente, em A/mm2.

openfemm;                                    %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/ax_magnetico.fem']);  %Abre o arquivo modelo, desenvolvdido pelo autor.
catch
    opendocument('ax_magnetico.fem');
end
mi_saveas('temp_ind.fem');                   %Aqui salva como um arquivo temporário, para não
                                             %destruir o modelo.

%% Seção Transversal dos Condutores
disp('Análise da Seção Transversal dos Condutores (Litz)...');

mi_probdef(frequency,'millimeters','axi',1e-8,0,20,(0));

u0 = 4*pi*1e-7;                                             %Permeabilidade magnética.
ur = 1;                                                     %Permeabilidade magnética relativa.
ro = 1.72e-8;                                               %Resistividade do cobre, em Ohm.m.

delta = 1e3*sqrt(2*ro/((2*pi*frequency)*ur*u0));            %Profundidade de penetração para a freq. de operação, em mm.
awgEquivalent = round((log(2*delta/8.2514694))/(-0.115943));          %Profundidade de penetração, em AWG, arredondado para o inteiro mais próximo.
awgEquivalent = awgEquivalent + 1;                          %O condutor escolhido deve ser superior ao mínimo calculado na linha acima.

if awgEquivalent < 0
   awgEquivalent = 0;
end

awg = awgEquivalent;

counter = 0;

% Enrolamento Emissor
secaop = (emiOpenCurrent/sqrt(2))/Jp;                                  %Seção transversão do condutor, considerando valor eficaz da corrente no emissor.
rcondp = sqrt(secaop/pi);                                   %Raio do condutor circular necessário para o emissor.
% Enrolamento Receptor
Js = Jp;                                                    %Densidade de Corrente, em A/mm2.
secaos = (recShortCurrent/sqrt(2))/Js;                                  %Seção transversão do condutor, considerando valor eficaz da corrente no receptor.
rconds = sqrt(secaos/pi);                                   %Raio do condutor circular necessário para o receptor.

awgMax = awg + 6;                                           %Número máximo de sugestões de bitolas de condutor.

autoWire = input('Digite 0 para determinar a bitola do condutor ou 1 para automático:');

if autoWire == 0
    
    awgDesired = input('Bitola do condutor desejado (AWG):');
    nStrands = input('Número de condutores em paralelo:');

    awg = awgDesired;       %This garantees that the "for" below is performed only for the desired wire gauge.
    awgMax = awg;
    
end

    
for wire = awg:2:awgMax
    iterCounter1 = 0;
    adjFactor = 1;
    
    counter = counter+1;
    rskin = (8.2514694*exp(-0.115943*wire))/2;                  %Converte de awg para mm depois do arredondamento. Este é o raio do condutor que atende
                                                                %a profundidade de penetração calculada.
    nlitzp_min = ceil((rcondp/rskin)^2);                            %Número mínimo de filamentos do fio Litz, arredondado para o maior inteiro.
    nlitzs_min = ceil((rconds/rskin)^2);                            %Número mínimo de filamentos do fio Litz, arredondado para o maior inteiro.

    
    if autoWire == 1
        nlitzp = nlitzp_min;                                        %ToDo: não precisa usar o mínimo numero de filamentos. Quanto mais melhor.
        nlitzs = nlitzs_min;                                        %Aqui poderia ajustar numero de filamentos até atingir alguma restrição mecânica.
    elseif autoWire == 0
        nlitzp = nStrands;
        nlitzs = nStrands;
    end
    
    mi_deletematerial('Litz Emissor');                          %Deleta condutor Litz das iteações anteriores.
    mi_deletematerial('Litz Receptor');
    mi_addmaterial('Litz Emissor',1,1,0,0,58,0,0,1,5,0,0,nlitzp,2*rskin); %Cria condutor Litz com as especificações da iteração atual.
    mi_addmaterial('Litz Receptor',1,1,0,0,58,0,0,1,5,0,0,nlitzs,2*rskin);

%% Análise Inicial dos Enrolamentos                          % Esta rotina ajusta o fator de utilização dos enrolamentos
                                                          
    mi_selectlabel(x1,y1);                                       %Os comandos iniciados por "mi_" referem-se à discretização de elementos magnéticos (m)
    mi_setblockprop('Litz Emissor',0,mesh,'emitter',0,1,naux);      %na etapa inicial (i), ou seja, no pré-processamento.
    mi_clearselected;

    mi_selectlabel(x2,y2);
    mi_setblockprop('Litz Receptor',0,mesh,'receiver',0,2,naux);
    mi_clearselected;
    mi_setcurrent('emitter',testCurrent);             %Força uma corrente de teste no circuito primário do modelo do FEMM.
    mi_setcurrent('receiver',testCurrent); 
    mi_analyze;
    mi_loadsolution;
    emiWindCalc = mo_getfill(x1,y1);
    recWindCalc = mo_getfill(x2,y2);

    recCoilHeight = (yfill2 - y2_limite_inf);
%     expectedCoilHeight = recCoilHeight;
    emiCoilHeight = (yfill1 - y1_limite_inf);

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
        mi_movetranslate(0,recWindStep);
        yfill2 = yfill2 + recWindStep;
        recCoilHeight = (yfill2 - y2_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        recWindCalc = mo_getfill(x2,y2);
        iterCounter1 = iterCounter1 + 1;
        
        recWindCalcError = (targetWind - recWindCalc)/targetWind;
        recWindStep = -recWindCalcError*adjFactor;
        
        if iterCounter1 > 15
            adjFactor = adjFactor/2;
            iterCounter1 = 0;
        end
        
    end

    recCoilHeight = (yfill2 - y2_limite_inf);
    emiWindCalcError= (targetWind - emiWindCalc)/targetWind;
    adjFactor = 1;
    iterCounter1 = 0;
    
    while (abs(emiWindCalcError) > emiWindError)                %Ajusta espessura do enrolamento sec. até que
                                                         %o fator de utilização seja maior que 70%.
                                                   
        if ((xfill1 - emiCoilHeight/2) <= xfill1/10)
            mi_selectgroup(1);
            mi_movetranslate(emiCoilHeight,0);
            xfill1 = xfill1 + emiCoilHeight;
            x1 = x1 + emiCoilHeight;
        end  
    
        mi_selectnode(xfill1,yfill1);
        mi_movetranslate(0,emiWindStep);
        yfill1 = yfill1 + emiWindStep;
        emiCoilHeight = (yfill1 - y1_limite_inf);
         
        mi_analyze;
        mi_loadsolution;
        emiWindCalc = mo_getfill(x1,y1);
        iterCounter1 = iterCounter1 + 1;
        
        emiWindCalcError= (targetWind - emiWindCalc)/targetWind;
        emiWindStep = -emiWindCalcError*adjFactor;                   %O erro é usado como passo de ajuste.
        
         if iterCounter1 > 15
            adjFactor = adjFactor/2;
            iterCounter1 = 0;
        end
        
    end

    emiCoilHeight = (yfill1 - y1_limite_inf);

    y1_limite_sup = yfill1;  %Enrolamento sec. não pode chegar nesta altura ou tocará o enrolamento primário.

    a = 0;
    
     mi_probdef(0,'millimeters','axi',1e-8,0,20,(0));       %Realiza primeira simulação para determinar resistência DC, para o cálculo teórico.
     mi_analyze;
     mi_loadsolution;
     h = mo_getcircuitproperties('emitter');
     emidcResistance = (real(h(2)/h(1)))/naux;                 %Divide a resistência série por naux porque a simulação foi feita para naux espiras.
                                                            %Assim se obtém a resistência por espira.    
    frequencyInc = frequency/5;
    
    for testFrequency = 0.5*frequency:frequencyInc:1.5*frequency
               
        a = a+1;
        
%         delta = 1e3*sqrt(2*ro/((2*pi*testFrequency)*ur*u0));                            %Calcula prof. de penetração teórica para cada frequencia, em mm.
%         rp_ac(a) = ((pi*rskin^2)/((pi*rskin^2)- pi*(rskin-delta)^2))*emidcResistance/nlitzp;  %Cálculo teórico de Rp AC.
                
        frequencyVector(a) = testFrequency; 
        mi_probdef(testFrequency,'millimeters','axi',1e-8,0,20,(0));
        mi_analyze;
        mi_loadsolution;
        h = mo_getcircuitproperties('emitter');
        emiResistance(a) = (real(h(2)/h(1)))/(naux*2*pi*x1*1e-3);  %Divide a resistência série por naux porque a simulação foi feita para naux espiras.
                                                                 %Assim se obtém a resistência por metro.
        v = mo_getcircuitproperties('receiver');
        recResistance(a) = (real(v(2)/v(1)))/(naux*2*pi*x2*1e-3);
        
        if counter ==1
            vec_rp_awg1 = emiResistance;
            vec_rs_awg1 = recResistance;
%             vec_rpac_awg1 = rp_ac;
            option1p = sprintf('1) %d condutores %d AWG',nlitzp,wire);
            legn1p = sprintf('%d condutores %d AWG',nlitzp,wire);
            option1s = sprintf('1) %d condutores %d AWG',nlitzs,wire);
            legn1s = sprintf('%d condutores %d AWG',nlitzs,wire);
        elseif counter == 2
            vec_rp_awg2 = emiResistance;
            vec_rs_awg2 = recResistance;
%             vec_rpac_awg2 = rp_ac;
            option2p = sprintf('2) %d condutores %d AWG',nlitzp,wire);
            legn2p = sprintf('%d condutores %d AWG',nlitzp,wire);
            option2s = sprintf('2) %d condutores %d AWG',nlitzs,wire);
            legn2s = sprintf('%d condutores %d AWG',nlitzs,wire);
        elseif counter == 3
            vec_rp_awg3 = emiResistance;
            vec_rs_awg3 = recResistance;
%             vec_rpac_awg3 = rp_ac;
            option3p = sprintf('3) %d condutores %d AWG',nlitzp,wire);
            legn3p = sprintf('%d condutores %d AWG',nlitzp,wire);
            option3s = sprintf('3) %d condutores %d AWG',nlitzs,wire);
            legn3s = sprintf('%d condutores %d AWG',nlitzs,wire);
        elseif counter == 4
            vec_rp_awg4 = emiResistance;
            vec_rs_awg4 = recResistance;
%             vec_rpac_awg4 = rp_ac;
            option4p = sprintf('4) %d condutores %d AWG',nlitzp,wire);
            legn4p = sprintf('%d condutores %d AWG',nlitzp,wire);
            option4s = sprintf('4) %d condutores %d AWG',nlitzs,wire);
            legn4s = sprintf('%d condutores %d AWG',nlitzs,wire);
        end
        
    end
    
end
   hold on

   if autoWire ==1
   figure(1);
   plot(frequencyVector,vec_rp_awg1*1e3,'k',frequencyVector,vec_rp_awg2*1e3,'b',frequencyVector,vec_rp_awg3*1e3,'r',frequencyVector,vec_rp_awg4*1e3,'g');
%    plot(frequencyVector,vec_rpac_awg1*1e3,'kx',frequencyVector,vec_rpac_awg2*1e3,'bx',frequencyVector,vec_rpac_awg3*1e3,'rx',frequencyVector,vec_rpac_awg4*1e3,'gx');
   xlabel('Freqüência (Hz)')
   ylabel('Resistência Série,(mOhm/m)')
   legend(legn1p,legn2p,legn3p,legn4p)   
     
   disp('Com base no gráfico de resist. série vs. freq., qual condutor será utilizado no emissor?:')
   disp(option1p)
   disp(option2p)
   disp(option3p)
   disp(option4p)
   choicep = input('Digite o número da opção desejada:');
   
    if choicep == 1
            A      = sscanf(legn1p,['%d condutores %d AWG'])';
            nlitzp = A(1);
            wirep  = A(2);
        elseif choicep == 2
            A      = sscanf(legn2p,['%d condutores %d AWG'])';
            nlitzp = A(1);
            wirep  = A(2);
        elseif choicep == 3
            A      = sscanf(legn3p,['%d condutores %d AWG'])';
            nlitzp = A(1);
            wirep  = A(2);
        elseif choicep == 4
            A      = sscanf(legn4p,['%d condutores %d AWG'])';
            nlitzp = A(1);
            wirep  = A(2);
    end
    
   figure(2);
   plot(frequencyVector,vec_rs_awg1*1e3,'k',frequencyVector,vec_rs_awg2*1e3,'b',frequencyVector,vec_rs_awg3*1e3,'r',frequencyVector,vec_rs_awg4*1e3,'g');
%    plot(frequencyVector,vec_rpac_awg1*1e3,'kx',frequencyVector,vec_rpac_awg2*1e3,'bx',frequencyVector,vec_rpac_awg3*1e3,'rx',frequencyVector,vec_rpac_awg4*1e3,'gx');
   xlabel('Freqüência (Hz)')
   ylabel('Resistência Série,(mOhm/m))')
   legend(legn1s,legn2s,legn3s,legn4s)  
    
    
   disp('Com base no gráfico de resist. série vs. freq., qual condutor será utilizado no receptor?:')
   disp(option1s)
   disp(option2s)
   disp(option3s)
   disp(option4s)
   choices = input('Digite o número da opção desejada:');
   
    if choices == 1
            B      = sscanf(legn1s,['%d condutores %d AWG'])';
            nlitzs = B(1);
            wires  = B(2);
        elseif choices == 2
            B      = sscanf(legn2s,['%d condutores %d AWG'])';
            nlitzs = B(1);
            wires  = B(2);
        elseif choices == 3
            B      = sscanf(legn3s,['%d condutores %d AWG'])';
            nlitzs = B(1);
            wires  = B(2);
        elseif choices == 4
            B      = sscanf(legn4s,['%d condutores %d AWG'])';
            nlitzs = B(1);
            wires  = B(2);
    end
   
   else
       
       figure(1);
       plot(frequencyVector,vec_rp_awg1*1e3,'k')
       xlabel('Freqüência (Hz)')
       ylabel('Resistência Série,(mOhm/m)')
       legend(legn1p);
       
%        figure(2);
%        plot(frequencyVector,vec_rs_awg1*1e3,'k')
%        xlabel('Freqüência (Hz)')
%        ylabel('Resistência Série norm. por espira, Rs/esp (mOhms)')
%        legend(legn1s);
       
       wirep = awgDesired;
       wires = wirep;
       
   end
   
   rwirep = (8.2514694*exp(-0.115943*wirep))/2;                %AWG para mm.
   rwires = (8.2514694*exp(-0.115943*wires))/2;
   mi_deletematerial('Litz Emissor');                          %Deleta condutor Litz das iteações anteriores.
   mi_deletematerial('Litz Receptor');
   mi_addmaterial('Litz Emissor',1,1,0,0,58,0,0,1,5,0,0,nlitzp,2*rwirep); %Cria condutor Litz com as especificações da iteração atual.
   mi_addmaterial('Litz Receptor',1,1,0,0,58,0,0,1,5,0,0,nlitzs,2*rwires);

%Atualiza FEMM com opções escolhidas pelo usuário.   
   mi_probdef(frequency,'millimeters','axi',1e-8,0,20,(0));
   mi_selectlabel(x1,y1);                                       
   mi_setblockprop('Litz Emissor',0,mesh,'emitter',0,1,1);     
   mi_clearselected;

   mi_selectlabel(x2,y2);
   mi_setblockprop('Litz Receptor',0,mesh,'receiver',0,2,1);
   mi_clearselected;
   
   wind(1) = x1;
   wind(2) = x2;
   wind(3) = y1;                                                                   
   wind(4) = y2;
   wind(5) = y1_limite_inf;
   wind(6) = y2_limite_inf;
   wind(7) = xfill2;
   wind(8) = yfill2;
   wind(9) = xfill1;
   wind(10) = yfill1;
   
mi_saveas('result_wind.fem'); %Aqui salva arquivo .fem com representação dos indutores.
closefemm

prt = [frequency/1e3; wirep; nlitzp; wires; nlitzs];        %Matriz de parâmetros básicos da geometria.

delete('temp_ind.fem');
delete('temp_ind.ans');

end