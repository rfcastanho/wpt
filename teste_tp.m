clc;
clear;
f = 15e3;
iab = 15;
Ns = 30;

%% Análise Temporal do Sistema Fracamente Acoplado
disp('Iniciando Fase 4 - Análise Temporal do Sistema ...');

openfemm;      %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/resultado_ax_esp.fem']);  %Abre o arquivo salvo anteriormente.
catch
    opendocument('resultado_ax_esp');
end

freq = 0;
mi_probdef(freq,'millimeters','axi',1e-8,0,20,(0)); %Seta simulação com frequencia = 0, para tratar o problema
                                                    %como magnetostático.
w = 2*pi*f;                                                    
mi_setcurrent('primario',0);                        %Força correntes nulas
mi_setcurrent('secundario',0); 
n = 0;
nsim = 50;  %Número de simulações por período da forma de onda
nciclos = 1; %Número de periodos a serem simulados
fluxo_anterior_sec = 0;

%Simulação com conteúdo harmônico em ip(t).

mkdir ('Figuras');

%Componente fundamental:
for tp = 0:1/(nsim*f):nciclos/f    %Passo de simulação = período de acordo com a freq. informada no início da representação axisimétrica.
    
    n = n+1;
    ip = iab*sin(w*tp);%+(iab/3)*sin(3*w*tp)+(iab/5)*sin(5*w*tp);
    iprim(n) = ip;
    tempo(n) = tp;
    mi_setcurrent('primario',ip); 
    mi_analyze(1);
    mi_loadsolution;
    mo_maximize;
    mo_zoom(0,-60,180,90);
    mo_hidecontourplot;
    mo_showdensityplot(0,0,0.01,0,'mag');
    mo_savebitmap(sprintf('result_time%d.bmp',n));
    
    movefile(sprintf('result_time%d.bmp',n),'Figuras');
    
    r=mo_getcircuitproperties('primario');
    fluxo_real_prim(n) = real(r(3)); 
    
    s=mo_getcircuitproperties('secundario');
    fluxo_real_sec(n) = real(s(3)); 
    
    delta_fluxo_sec = fluxo_real_sec(n) - fluxo_anterior_sec;  %Cálculo da tensão induzida no secundário em aberto.
    fluxo_anterior_sec = fluxo_real_sec(n);
    vsab(n) = -delta_fluxo_sec*nsim*f;      %Lei de Faraday-Neumann-Lenz, Vsa = -dfluxo/dt;
    
end

plot(tempo,iprim);
hold;
plot(tempo,10e3*fluxo_real_prim,'r');% Fluxo real no primário em mWb.
plot(tempo,vsab,'k');

