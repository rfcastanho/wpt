%Emitter Compensation Capacitance

function capFunction = emicap (emiType,recType,recCalc,emiCalc,mutualCalc,Load,frequency)

w = 2*pi*frequency;
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

capFunction(1) = emiCapacitance;
capFunction(2) = recCapacitance;

end
