%##########################################################################
%File: retriever.m                                                        #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: November 16, 2013                                           #
%Description: Data retriever, a post-processor.                           #
%##########################################################################

clear all;
clc;
cloudH = [];
cloudB = [];
n = 1;  %Figure index.
retrieveTime = 0;
Bm = [];
Hm = [];


openProject = input('Digite o nome do projeto a ser processado: ', 's');
projName = sprintf('%s.fem',openProject);
projName2 = sprintf('/%s.fem',openProject);

openfemm;                                         %Matlab abre o aplicativo FEMM
try
    vv=ver;
    opendocument([cd,projName2]); 
catch
    opendocument(projName);
end
mi_saveas('temp_ret.fem');

mi_analyze(1);
mi_loadsolution;
v = mo_getprobleminfo;
simType = v(1);
depth = v(3);

if simType == 0
    simType = 'planar';
    displacement = 0;
elseif simType == 1
    simType = 'axi';
    displacement = 1;
end

clc;
disp('***** Ponto de Operação *****')
disp('A análise de Campo Magnético pode ser feita para qualquer ponto de operação')
disp('do conversor. Entretanto, é interessante obter resultados para as condições')
disp('mais críticas.')
emiCurrent = input('Corrente no emissor (A):');
recCurrent = input('Corrente no receptor (A):');
frequency = input('Freqüência de operação (kHz):');
frequency = frequency*1e3;

mi_setcurrent('emitter',emiCurrent);
mi_setcurrent('receiver',-recCurrent); 
mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));
mi_analyze(1);
mi_loadsolution;

% Tabela de Limites de Campo Elétrico (Eref), Intensidade de campo
% magnético (Href) e densidade de fluxo magnético (Bref), de acordo com
% Tabela 7 de ICNIRP GUIDELINES FOR LIMITING EXPOSURE TO TIME VARYING 
% ELECTRIC, MAGNETIC AND ELECTROMAGNETIC FIELDS (UP TO 300 GHZ), de 1998.

if frequency < 1
        Eref = 1e10; Href = 3.2e4; Bref = 4e4;
elseif frequency >= 1 && frequency < 8
        Eref = 10e3; Href = 3.2e4/(frequency^2); Bref = (4e4/(frequency^2));
elseif frequency >= 8 && frequency < 25
        Eref = 10e3; Href = 4000/frequency; Bref = 5000/frequency;
elseif frequency >= 25 && frequency < 800
        Eref = 250/(frequency/1e3); Href = 4/(frequency/1e3); Bref = 5/(frequency/1e3); %Em todos os casos, Bref = u0*Href, com u0 = 4pi*10^-7 H/m.
elseif frequency >= 800 && frequency < 3000
         Eref = 250/(frequency/1e3); Href = 5; Bref = 6.25;
elseif frequency >= 3000 && frequency < 150e3
        Eref = 87; Href = 5; Bref = 6.25;
elseif frequency >= 150e3 && frequency < 1e6
        Eref = 87; Href = 0.73/(frequency/1e6); Bref = 0.92/(frequency/1e6);
elseif frequency >= 1e6 && frequency < 10e6
        Eref = 87/((frequency/1e6)^0.5); Href = 0.73/(frequency/1e6); Bref = 0.92/(frequency/1e6);
elseif frequency >= 10e6
        Eref = 28; Href = 0.073; Bref = 0.092;               
end

Eref = sqrt(2)*Eref;
Href = sqrt(2)*Href;
Bref = sqrt(2)*Bref*1e-6;   %Bref é expresso em uT na ICNIRP 1998.

%% Determinação da janela de simulação do arquivo aberto
clc;
disp('***************************************************************************************');
disp('Para melhor representação do domínio de simulação, determine os limites manualmente,');
disp('dando foco à região de interesse. Ao utilizar o modo automático, todo o domínio de simulação');
disp('será analisado, resultando em grande aumento no tempo de execução da rotina.');
disp('Além disso, o arquivo .fem a ser analisado deve conter um domínio de simulação com dimensões');
disp('suficientes para caracterizar o sistema.');

auto = input ('Definir região de extração de dados manualmente(0) ou automaticamente(1)?:');

if auto ==1
    
    maxXY = mi_selectnode(10e3,10e3);   %Testa um ponto distante da janela de simulação, no primeiro quadrante.
    Xmax = maxXY(1); Ymax = maxXY(2);
    mi_clearselected;

    minXY = mi_selectnode(-10e3,-10e3);   %Testa um ponto distante da janela de simulação, no terceiro quadrante.
    Xmin = minXY(1); Ymin = minXY(2); 
    mi_clearselected;

elseif auto == 0
    
    Xmax = input('Qual o valor de Xmax (mm)?:');
    Ymax = input('Qual o valor de Ymax (mm)?:');
    Xmin = input('Qual o valor de Xmin (mm)?:');
    Ymin = input('Qual o valor de Ymin (mm)?:');
end

windowArea = (Xmax - Xmin)*(Ymax - Ymin);

%% Definição das matrizes de dados

clc;
disp('***************************************************************************************');
sprintf('A domínio de simulação possui área de %d mm2',windowArea);
disp('Esta rotina pode extrair resultados da densidade de fluxo magnético (B) e da intensidade de');
disp('campo magnético (H) para um grande número de pontos do domínio de simulação. Como B = uH,');
disp('a análise de apenas uma das duas grandezas é suficiente para a maioria das aplicações de projeto.');
disp('O tempo de execução desta rotina aumentará se for necessário extrair dados de B e H simultaneamente.');
disp('Além disso, a resolução para extração de dados pode elevar o tempo de execução. Para uma análise rápida,');
disp('resoluções de 2 a 5 mm podem ser suficientes.');

Bretrieve = input('Digite (1) para obter resultados de B(T) ou (0) para ignorar:');
Hretrieve = input('Digite (1) para obter resultados de H(A/m) ou (0) para ignorar:');
resolution = input('Qual a distância entre pontos a serem salvos nas matrizes de dados (ex.: 1 mm)?:');

for y = Ymin:resolution:Ymax 
    
    if Ymin < 0
        yoffset = y - Ymin;               %Faz offset das coordenadas para evitar indice negativo.
    else
        yoffset = y;
    end
    
    estimatedTime = retrieveTime*(Ymax-y);           %Tempo estimado para fim da obtenção de dados, em minutos.
    clc;
    disp(sprintf('Tempo estimado para fim da obtenção de dados: %4.1f segundos',estimatedTime));
    
    clockY = tic;
    for x = Xmin:resolution:Xmax
        
        if Xmin < 0
           xoffset = x - Xmin;
        else
            xoffset = x;
        end
                
        if Bretrieve ==1
            B = mo_getb(x,y);
            Bx= abs(B(1));      
            By= abs(B(2));
            Bm(xoffset+1,yoffset+1) = sqrt(Bx^2 + By^2); %Soma 1 para evitar indice (0,0).
            
            Baux =  Bm(xoffset+1,yoffset+1);
            if Baux >= Bref
                cloudB = [cloudB; [xoffset+1 yoffset+1]];
            end
            
        end
        
        if Hretrieve == 1
            H = mo_geth(x,y);
            Hx= abs(H(1));      
            Hy= abs(H(2));
            Hm(xoffset+1,yoffset+1) = sqrt(Hx^2 + Hy^2);
            Haux =  Hm(xoffset+1,yoffset+1);
            if Haux >= Href
                cloudH = [cloudH; [xoffset+1 yoffset+1]];
            end
            
        end     
                 
    end
    retrieveTime = toc(clockY);        
end

Bm = Bm';
Hm = Hm';


clc;
disp('*****************************************************');
disp('Somente é possível plotar resultados se eles tiverem sido coletados no passo anterior.');
disp('Se escolhida, a exibição dos resultados pode levar alguns minutos.');


if Bretrieve == 1
        
    showB = input('Digite (0) para plotar resultado simplificado de B(T), (1) para detalhado ou (2) para ambos:');   
    
    if (showB == 1 || showB == 2)
        disp('Aguarde...');
        figure(n)
        contourf(Bm,'LineStyle','none');          %,'LevelStep', 0.00002,'LineStyle','none');
        colormap jet; colorbar;
        xlabel('Simulation Domain Width (mm)'); ylabel('Simulation Domain Length (mm)');
        title('Magnetic Flux Density Map (|B|,T)');
        n = n+1;    %Atualiza indice da figura.
    end
    
    if (showB == 0 || showB == 2)
        figure(n)
        scatter(cloudB(:,1),cloudB(:,2),'*');
        xlabel('Simulation Domain Width (mm)'); ylabel('Simulation Domain Length (mm)');
        title('Magnetic Flux Density Boundary (|B|> Bref,T)');
        n = n+1;    %Atualiza indice da figura.
    end
 
limitx = max(cloudB(:,1));
limity = max(cloudB(:,2));
    
end

if Hretrieve == 1
     
    showH = input('Digite (0) para plotar resultado simplificado de H(A/m), (1) para detalhado ou (2) para ambos:');   
    
    if (showH == 1 || showH == 2)  
        disp('Aguarde...');
        figure(n)
        contourf(Hm,'LineStyle','none'); %,'LevelStep', 0.00002,'LineStyle','none');
        colormap jet; colorbar;
        xlabel('Simulation Domain Width (mm)'); ylabel('Simulation Domain Length (mm)');
        title('Magnetic Field Intensity Map (|H|,A/m)');
        n = n+1;    %Atualiza indice da figura.
    end
    
    if (showH == 0 || showH == 2)  
        
        figure(n)
        scatter(cloudH(:,1),cloudH(:,2),'*');
        xlabel('Simulation Domain Width (mm)'); ylabel('Simulation Domain Length (mm)');
        title('Magnetic Field Intensity Boundary (|H| > Href,A/m)');
        n = n+1;    %Atualiza indice da figura.
    end
    
limitx = max(cloudH(:,1));
limity = max(cloudH(:,2));

end

if (Hretrieve ==1 || Bretrieve ==1)
    clc;
    disp('OBS: Ajuste o número de isolinhas (LevelStep) das imagens geradas, para melhor análise do problema.');
end

    disp('***********************************************************************************************');
    disp('De maneira simplificada, para este conversor os limites de segurança com base na ICNIRP 1998');
    disp(sprintf('são atendidos para x > %d mm e y > %d mm', limitx,limity));

closefemm;