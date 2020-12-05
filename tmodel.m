function modelTransf = tmodel (rint,emiType,recType,emiCalc,mutualCalc,emiCapacitance,w,emiResistance,recResistance,recCalc,recCapacitance,Load,emiVoltage,emiCurrent,mode,recCurrent)

emiCurrentVariation = strcmp(mode,'emiCurrentVar');
recCurrentVariation = strcmp(mode,'recCurrentVar');

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

if emiCurrentVariation == 0
    
    emiCurrent = (m*emiVoltage)/(c1 + (w^2*mutualCalc^2)/c2);           %Emitter current will be calculated only if emitter input voltage is not null.
%     if abs(emiCurrent) > 40                 %Current limitation
%         emiAng = angle(emiCurrent);
%         emiAbs = 40;
%         emiCurrent = emiAbs*cos(emiAng) + 1i*emiAbs*sin(emiAng);
%     end
end

if recCurrentVariation == 0
    recCurrent = 1i*w*mutualCalc*emiCurrent/(c2);
end 

if recType == 1
    loadCurrent = recCurrent; 
elseif recType == 0
    loadVoltage = (Load*recCurrent)/(1 + 1i*w*Load*recCapacitance);
    loadCurrent = recCurrent - 1i*w*recCapacitance*loadVoltage;
end

modelTransf(1) = emiCurrent;
modelTransf(2) = recCurrent;
modelTransf(3) = loadCurrent;

end