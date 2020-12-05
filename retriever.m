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
disp('***** Ponto de Opera��o *****')
disp('A an�lise de Campo Magn�tico pode ser feita para qualquer ponto de opera��o')
disp('do conversor. Entretanto, � interessante obter resultados para as condi��es')
disp('mais cr�ticas.')
emiCurrent = input('Corrente no emissor (A):');
recCurrent = input('Corrente no receptor (A):');
frequency = input('Freq��ncia de opera��o (kHz):');
frequency = frequency*1e3;

mi_setcurrent('emitter',emiCurrent);
mi_setcurrent('receiver',-recCurrent); 
mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));
mi_analyze(1);
mi_loadsolution;

% Tabela de Limites de Campo El�trico (Eref), Intensidade de campo
% magn�tico (Href) e densidade de fluxo magn�tico (Bref), de acordo com
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
Bref = sqrt(2)*Bref*1e-6;   %Bref � expresso em uT na ICNIRP 1998.

%% Determina��o da janela de simula��o do arquivo aberto
clc;
disp('***************************************************************************************');
disp('Para melhor representa��o do dom�nio de simula��o, determine os limites manualmente,');
disp('dando foco � regi�o de interesse. Ao utilizar o modo autom�tico, todo o dom�nio de simula��o');
disp('ser� analisado, resultando em grande aumento no tempo de execu��o da rotina.');
disp('Al�m disso, o arquivo .fem a ser analisado deve conter um dom�nio de simula��o com dimens�es');
disp('suficientes para caracterizar o sistema.');

auto = input ('Definir regi�o de extra��o de dados manualmente(0) ou automaticamente(1)?:');

if auto ==1
    
    maxXY = mi_selectnode(10e3,10e3);   %Testa um ponto distante da janela de simula��o, no primeiro quadrante.
    Xmax = maxXY(1); Ymax = maxXY(2);
    mi_clearselected;

    minXY = mi_selectnode(-10e3,-10e3);   %Testa um ponto distante da janela de simula��o, no terceiro quadrante.
    Xmin = minXY(1); Ymin = minXY(2); 
    mi_clearselected;

elseif auto == 0
    
    Xmax = input('Qual o valor de Xmax (mm)?:');
    Ymax = input('Qual o valor de Ymax (mm)?:');
    Xmin = input('Qual o valor de Xmin (mm)?:');
    Ymin = input('Qual o valor de Ymin (mm)?:');
end

windowArea = (Xmax - Xmin)*(Ymax - Ymin);

%% Defini��o das matrizes de dados

clc;
disp('***************************************************************************************');
sprintf('A dom�nio de simula��o possui �rea de %d mm2',windowArea);
disp('Esta rotina pode extrair resultados da densidade de fluxo magn�tico (B) e da intensidade de');
disp('campo magn�tico (H) para um grande n�mero de pontos do dom�nio de simula��o. Como B = uH,');
disp('a an�lise de apenas uma das duas grandezas � suficiente para a maioria das aplica��es de projeto.');
disp('O tempo de execu��o desta rotina aumentar� se for necess�rio extrair dados de B e H simultaneamente.');
disp('Al�m disso, a resolu��o para extra��o de dados pode elevar o tempo de execu��o. Para uma an�lise r�pida,');
disp('resolu��es de 2 a 5 mm podem ser suficientes.');

Bretrieve = input('Digite (1) para obter resultados de B(T) ou (0) para ignorar:');
Hretrieve = input('Digite (1) para obter resultados de H(A/m) ou (0) para ignorar:');
resolution = input('Qual a dist�ncia entre pontos a serem salvos nas matrizes de dados (ex.: 1 mm)?:');

for y = Ymin:resolution:Ymax 
    
    if Ymin < 0
        yoffset = y - Ymin;               %Faz offset das coordenadas para evitar indice negativo.
    else
        yoffset = y;
    end
    
    estimatedTime = retrieveTime*(Ymax-y);           %Tempo estimado para fim da obten��o de dados, em minutos.
    clc;
    disp(sprintf('Tempo estimado para fim da obten��o de dados: %4.1f segundos',estimatedTime));
    
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
disp('Somente � poss�vel plotar resultados se eles tiverem sido coletados no passo anterior.');
disp('Se escolhida, a exibi��o dos resultados pode levar alguns minutos.');


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
    disp('OBS: Ajuste o n�mero de isolinhas (LevelStep) das imagens geradas, para melhor an�lise do problema.');
end

    disp('***********************************************************************************************');
    disp('De maneira simplificada, para este conversor os limites de seguran�a com base na ICNIRP 1998');
    disp(sprintf('s�o atendidos para x > %d mm e y > %d mm', limitx,limity));

closefemm;