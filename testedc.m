clear;
clc;


f = 2000;
rp = 50;
dist_e = 1;
Lp_calc = 1.03e-3;

openfemm;
try
    vv=ver;
    opendocument([cd,'/exp_planar.fem']); %Abre geometria representada por simetria.
catch
    opendocument('exp_planar.fem');       %Abre geometria representada por simetria.
end
mi_saveas('temp_exp.fem');               %Salva como temporário para manipular.


%% Estudo Variando dc (enrolamentos não-coaxiais)
%Este algoritmo plota curvas de fator de acoplamento em função de dc
%(distância entre os centros do primário e secundário, para diferentes
%distâncias e.

disp('Iniciando Estudo de Enrolamentos Não-Coaxiais...');
mi_setcurrent('secundario',0);      %Força corrente nula no circuito secundário do modelo do FEMM.
mi_setcurrent('primario',2);      %Força corrente com secundário em aberto no modelo do FEMM. 

delta_dc = 20;                  %Passo de dc para realizar as simulações. Quanto maior, mais detalhado o estudo.
cont_curvas = 0;               %Contador de quantas curvas (k x dc) foram traçadas
  
n_curvas = 1;     %RETIRAR ESTA LINHA. É SÓ TESTE.
  
  

  n = 0;                                    %Reseta os contadores para reiniciar.
  cont = 0;  
  mi_analyze(1);
  mi_loadsolution;
    
for cont = 0:rp/delta_dc:4*rp                    %Move o secundário da esquerda para a direita.
  n = n + 1;
 
  h=mo_getcircuitproperties('secundario'); %Aqui calcula M e k.
  v2abplan_calc = abs(h(2));
  u=mo_getcircuitproperties('primario');
  Mplan_calc = abs(v2abplan_calc/(2*pi*f*u(1)));
  kplan_calc(n) = Mplan_calc/Lp_calc;      %sqrt(Lp*Ls) = Lp, pois Lp = Ls.
  dc(n) = cont; 
    
  mi_selectgroup(2);                      %Depois de calcular k para a posição espacial atual,
  mi_selectgroup(4);                      % move o secundário lateralmente e reinicia.
  mi_movetranslate(rp/delta_dc,0);
  mi_clearselected;
  mi_analyze(1);
  mi_loadsolution;
    
end
  
plot(dc,kplan_calc);                     %Plota a curva (k x dc) para uma distância e.

hold on


disp('Fim Estudo de Enrolamentos Não-Coaxiais...');

