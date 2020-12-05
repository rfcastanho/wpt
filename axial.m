function axialMis = axial(fileName,logDesired,emiOpenCurrent,frequency,y2_limite_inf_mi,y1_limite_sup,emiCalc,recCalc,usingShield,emiTurn,recTurn)

if nargin <11
    error('Faltam argumentos de entrada')
end
% clear;
% clc;
% emiOpenCurrent = 13;
% frequency = 22e3;
% recCalc = 115e-6;
% emiCalc = 180e-6;
% emiTurn = 45;
% recTurn = 47;
% usingShield = 1;

recPosition = 0;
iterCounter = 0;
testCurrent = 2;

openfemm;
try
    vv=ver;
    opendocument([cd,'/result_dimensions.fem']);
catch
    opendocument('result_dimensions.fem');
end
mi_saveas('temp_axial.fem');

 y2_limite_inf_mi = 41.7625;
 y1_limite_sup = -6.38054;
mi_selectgroup(2);
mi_selectgroup(4);
mi_movetranslate(0,-(y2_limite_inf_mi-y1_limite_sup-0.1));
mi_clearselected;

%% Fase 1

disp('Iniciando Estudo de Desalinhamento Axial ...');

distStep = input('Passo de incremento da distância (mm):');
maxIter = input('Número de incrementos:');

if maxIter == 0
    maxIter = 1;
end

mi_setcurrent('emitter',emiOpenCurrent); 
mi_setcurrent('receiver',0); 

while (iterCounter <= maxIter)

    iterCounter = iterCounter + 1;
    
    if usingShield == 1
        mi_setcurrent('receiver',testCurrent);  %If any kind of magnetic shield is used, self inductances must be calculated after
        mi_setcurrent('emitter',0);             % each variation in spatial position.
        mi_analyze;
        mi_loadsolution;
        rr=mo_getcircuitproperties('receiver');
        recCalc = abs(rr(3)/rr(1));
        recCalcVec(iterCounter) = recCalc;
        
        mi_setcurrent('emitter',emiOpenCurrent); 
        mi_setcurrent('receiver',0);
        mi_analyze;
        mi_loadsolution;
        r=mo_getcircuitproperties('emitter');
        emiCalc = abs(r(3)/r(1));
        emiCalcVec(iterCounter) = emiCalc;
    else
        mi_analyze;
        mi_loadsolution;
    end
    
    u=mo_getcircuitproperties('receiver');
    recOpenVoltage = abs(u(2));
    
    r=mo_getcircuitproperties('emitter');
    mutualCalc = abs(recOpenVoltage/(2*pi*frequency*emiOpenCurrent));
    emiLeak(iterCounter) = emiCalc -(emiTurn/recTurn)*mutualCalc;
    recLeak(iterCounter) = recCalc -(emiTurn/recTurn)*mutualCalc;
    
    couplingCalc(iterCounter) = mutualCalc/(sqrt(recCalc*emiCalc));
    dist(iterCounter) = recPosition; 
    mi_selectgroup(2);
    mi_selectgroup(4);
    mi_movetranslate(0,distStep)
    recPosition = recPosition + distStep;

end

figure;
plot(dist,couplingCalc,'r');
xlabel('Emitter-Receiver Axial Distance, e (mm)'); ylabel('Coupling Coefficient, k');
title('Axial Misalignment Analysis');

figure;
    plot(dist,emiLeak*1e6,'r');
    emiLeg = ('Emitter');
    hold on
    plot(dist,recLeak*1e6,'k');
    recLeg = ('Receiver');
    xlabel('Emitter-Receiver Axial Distance, e (mm)'); ylabel('Leakage Inductances, Lp, Ls (uH)');
    title('Self-Inductance Analysis');
    legend(emiLeg,recLeg);
    hold off

if usingShield == 1
    figure;
    plot(dist,emiCalcVec*1e6,'r');
    emiLeg = ('Emitter');
    hold on
    plot(dist,recCalcVec*1e6,'k');
    recLeg = ('Receiver');
    xlabel('Emitter-Receiver Axial Distance, e (mm)'); ylabel('Self-Inductances, Lp, Ls (uH)');
    title('Self-Inductance Analysis');
    legend(emiLeg,recLeg);
    hold off
end

closefemm;

delete('temp_axial.fem');
delete('temp_axial.ans');

axialMis(1) = 1;

end
