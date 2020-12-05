%Time analysis


function [ output_args ] = time_ans( input_args )

%% An�lise Temporal do Sistema Fracamente Acoplado
disp('Iniciando Fase 4 - An�lise Temporal do Sistema ...');

openfemm;      %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,'/resultado_ax_esp.fem']);  %Abre o arquivo salvo anteriormente.
catch
    opendocument('resultado_ax_esp');
end

freq = 0;
mi_probdef(freq,'millimeters','axi',1e-8,0,20,(0)); %Seta simula��o com frequencia = 0, para tratar o problema
                                                    %como magnetost�tico.
n = 0;

for tp = 0:1/(100*f):1/f    %Passo de simula��o = per�odo de acordo com a freq. informada no in�cio da representa��o axisim�trica.
    
    n = n+1;
    ip(n) = iab*cos(2*pi*f*tempo);
    tempo(n) = tp;
end

plot (tempo,ip);

pause

closefemm                   %Fecha FEMM.


end

