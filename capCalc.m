%Calculations for Sensitivity Analysis
clear; clc;

disp('Compensation types: 1 = Series, 0 = Parallel');
disp(' ');
disp('Enter values for nominal condition:');
emiType = input('Emitter Type: ');
recType = input('Receiver Type:  ');
emiCalc = input('Emitter Inductance: ');
recCalc = input('Receiver Inductance: ');
mutualCalc = input('Mutual Inductance: ');
emiResistance = input('Emitter Resistance: ');
recResistance = input('Receiver Resistance: ');
%emiVoltage = input('Input RMS Voltage: ');
Su = input('Su (VA): ');

emiCalc = emiCalc*1e-6;
recCalc = recCalc*1e-6;
mutualCalc = mutualCalc*1e-6;

coupling = mutualCalc/(sqrt(emiCalc*recCalc));
frequency = 38400;
emiCurrent = 15*sqrt(2);
rint = 11*1e-3;

w = 2*pi*frequency;
Sc = 500;       %Assuming that desired compensated power is 500 VA.
QFactor = Sc/Su;
Load = w*recCalc/QFactor;

recCapacitance = 1/(w^2*recCalc);

% REF1 = WANG, STIELAU, COVIC: DESIGN CONSIDERATIONS FOR A CONTACTLESS ELECTRIC VEHICLE BATTERY CHARGER, 
% IEEE TRANSACTIONS ON INDUSTRIAL ELECTRONICS, VOL. 52, NO. 5, OCTOBER 2005 

% REF2 = Thesis "Wireless Power Transfer do E-Mobility", by Venugopal Prasanth, 2012.

if (emiType == 1) && (recType == 1) %SS Topology
    
    emiCapacitance = recCalc*recCapacitance/emiCalc;  % REF1, table II, SS topology.
    emiCapacitance2 = 1/(w^2*emiCalc);                % REF2, page 27, table II, SS topology.
    
elseif (emiType == 1) && (recType == 0) %SP Topology
    
    emiCapacitance = recCalc^2*recCapacitance/(emiCalc*recCalc - mutualCalc^2); % REF1, table II, SP topology.
    emiCapacitance2 = 1/(w^2*(emiCalc - (mutualCalc^2/recCalc)));               % REF2, page 27, table II, SP topology.
    
elseif (emiType == 0) && (recType == 1) %PS Topology
    
    emiCapacitance = recCalc*recCapacitance/(emiCalc + (mutualCalc^4/(emiCalc*recCalc*recCapacitance*(real(Load))^2))); % REF1, table II, PS topology, but correcting R to R^2.
    emiCapacitance2 = emiCalc/((w^2*mutualCalc^2/real(Load))^2 + w^2*emiCalc^2);                                        % REF2, page 27, table II, PS topology.
    
elseif (emiType == 0) && (recType == 0) %PP Topology
    
    emiCapacitance = (recCalc^2*recCapacitance*(emiCalc*recCalc - mutualCalc^2))/(mutualCalc^4*recCapacitance*(real(Load))^2/recCalc + (emiCalc*recCalc - mutualCalc^2)^2); % REF1, table II, PP topology, but correcting R to R^2.
    emiCapacitance2 = (emiCalc - (mutualCalc^2/recCalc))/(((mutualCalc^2*real(Load))/recCalc^2)^2 + w^2*(emiCalc - (mutualCalc^2/recCalc))^2); % REF2, page 27, table II, PP topology.
end


    if emiType == 1
        c1 = emiResistance + 1i*w*(emiCalc-mutualCalc)+1/(1i*w*emiCapacitance)+1i*w*mutualCalc + rint;
        m = 1;
    elseif emiType == 0
        t = ((rint^2)/(1/(1i*w*emiCapacitance)+ rint) - rint);
        c1 = emiResistance + 1i*w*(emiCalc-mutualCalc)+1i*w*mutualCalc - t;
        m = 1 - (rint/(1/(1i*w*emiCapacitance)+ rint));
    end

    if recType == 1
        c2 = recResistance + 1i*w*(recCalc-mutualCalc)+1/(1i*w*recCapacitance)+Load+1i*w*mutualCalc;
    elseif recType == 0
        c2 = recResistance + 1i*w*(recCalc-mutualCalc)+Load/(1+1i*w*Load*recCapacitance)+1i*w*mutualCalc;
    end

    emiVoltage = emiCurrent*(c1 + (w^2*mutualCalc^2)/c2)/m;           %Emitter current will be calculated only if emitter input voltage is not null.

    load mutualVec.txt;
    A = mutualVec; 
    load emiCalcVec.txt;
    B = emiCalcVec;   
    load recCalcVec.txt;
    C = emiCalcVec;
    
for n = 1:length(mutualVec)
   
    mutualCalc = A(n);
    emiCalc    = B(n);
    recCalc    = C(n);
    
    if emiType == 1
        c1 = emiResistance + 1i*w*(emiCalc-mutualCalc)+1/(1i*w*emiCapacitance)+1i*w*mutualCalc + rint;
        m = 1;
    elseif emiType == 0
        t = ((rint^2)/(1/(1i*w*emiCapacitance)+ rint) - rint);
        c1 = emiResistance + 1i*w*(emiCalc-mutualCalc)+1i*w*mutualCalc - t;
        m = 1 - (rint/(1/(1i*w*emiCapacitance)+ rint));
    end

    if recType == 1
        c2 = recResistance + 1i*w*(recCalc-mutualCalc)+1/(1i*w*recCapacitance)+Load+1i*w*mutualCalc;
    elseif recType == 0
        c2 = recResistance + 1i*w*(recCalc-mutualCalc)+Load/(1+1i*w*Load*recCapacitance)+1i*w*mutualCalc;
    end

    emiCurrent = (m*emiVoltage)/(c1 + (w^2*mutualCalc^2)/c2);           %Emitter current will be calculated only if emitter input voltage is not null.
    recCurrent = 1i*w*mutualCalc*emiCurrent/(c2);

    if recType == 1
        loadCurrent = recCurrent; 
    elseif recType == 0
        loadVoltage = (Load*recCurrent)/(1 + 1i*w*Load*recCapacitance);
        loadCurrent = recCurrent - 1i*w*recCapacitance*loadVoltage;
    end
    
    SL(n) = Load*loadCurrent*conj(loadCurrent);
    
end



