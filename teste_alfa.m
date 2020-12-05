clear;
clc;

Lp_calc = 18.2e-6;
rp = 118.79;
rs = rp;
np = 5;
ns = np;
diamp = 25;
diams = diamp;
y1 = -8.05;
y2_limite_inf = 32.7625;
iteste = 2;
np_novo = np;
f = 15000;
Lp = 20e-6;
iab = 27;
n =0;
dist_e = 47.9;


openfemm;
try
    vv=ver;
    opendocument([cd,'/resultado_planar_esp.fem']); %Recupera geometria representada na forma planar.
catch
    opendocument('resultado_planar_esp.fem');     
end

mi_saveas('temp_planar_esp.fem');

%%Estudo da Varia��o de alfa (enrolamentos n�o paralelos)

%Parte 1- Aqui dc = 0, portanto prim�rio e secund�rio s�o coaxiais. Alpha �
%variado de 0 a 180 graus.

disp('Iniciando Estudo de Enrolamentos N�o-Paralelos...');
mi_setcurrent('secundario',0);      %For�a corrente nula no circuito secund�rio do modelo do FEMM.
mi_setcurrent('primario',iab);      %For�a corrente com secund�rio em aberto no modelo do FEMM. 

mi_selectgroup(4);
mi_setgroup(2);
mi_clearselected;

%Girando o secund�rio no sentido anti-hor�rio
alfa = 0;   %Angulo inicial nulo, enrolamentos paralelos.
delta_alfa = 5; %Incremento de alfa.
primeiro_plot2 =1;
dc = 0;

while dc <= 2*rs
  n = 0;
  cont = 0;
  mi_analyze(1);
  mi_loadsolution;

  for cont = 0:delta_alfa:180

  n = n+1;
  b=mo_getcircuitproperties('secundario'); %Aqui calcula M e k.
  vsaplan_calc = abs(b(2));
  w=mo_getcircuitproperties('primario');
  Mplan_calc = abs(vsaplan_calc/(2*pi*f*w(1)));
  kplan_calc(n) = Mplan_calc/Lp_calc;      %sqrt(Lp*Ls) = Lp, pois Lp = Ls.
  
  alfa(n) = cont;
  
  mi_selectgroup(2);
  mi_moverotate(-rs+dc,y2_limite_inf,delta_alfa)  %Gira o secund�rio em passos de 5 graus.
  mi_clearselected;                            %O centro de giro � o ponto inferior do lado direito do enrolamento secund�rio.

  mi_analyze(1);
  mi_loadsolution;

  end

plot(alfa,kplan_calc);

if primeiro_plot2 == 1                    %Flag de primeira passada usado para fazer o "hold"
hold                                     %apenas uma vez.
primeiro_plot2 = 0;
end

mi_selectgroup(2);                                %Retorna o secund�rio para a posi��o inicial.
mi_moverotate(-rs+dc,y2_limite_inf,-n*delta_alfa);
mi_clearselected; 

dc = dc + rs/5;                          %Move o secund�rio para a direita e reinicia.
mi_selectgroup(2);
mi_movetranslate(rs/5,0);
mi_clearselected;

end









