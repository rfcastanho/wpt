%Emitter Compensation Capacitance
clear; clc;

type = input('Enter "1" for SS config, "2" for SP, "3" for PS or "4" for PP:');

% emiCalc = 120e-6;
% recCalc = 190e-6;
% ka = 0.3;
% mutualCalc = ka*sqrt(emiCalc*recCalc); %Test configuration 1
% Load = 9.17;
% frequency = 38.4e3;
% w = 2*pi*frequency;

emiCalc = 120e-6;
recCalc = 190e-6;
ka = 0.408;
mutualCalc = ka*sqrt(emiCalc*recCalc);     %Test configuration 2
Load = 15.2807; %Era 3.62 para PP e SP; %Era 0.9 para o SS e PS.
frequency = 38400;
w = 2*pi*frequency;

recCapacitance = 1/(w^2*recCalc);

% REF1 = WANG, STIELAU, COVIC: DESIGN CONSIDERATIONS FOR A CONTACTLESS ELECTRIC VEHICLE BATTERY CHARGER, 
% IEEE TRANSACTIONS ON INDUSTRIAL ELECTRONICS, VOL. 52, NO. 5, OCTOBER 2005 

% REF2 = Thesis "Wireless Power Transfer do E-Mobility", by Venugopal Prasanth, 2012.

if type == 1
    
    emiCapacitance = recCalc*recCapacitance/emiCalc;  % REF1, table II, SS topology.
    emiCapacitance2 = 1/(w^2*emiCalc);                % REF2, page 27, table II, SS topology.
    
elseif type == 2
    
    emiCapacitance = recCalc^2*recCapacitance/(emiCalc*recCalc - mutualCalc^2); % REF1, table II, SP topology.
    emiCapacitance2 = 1/(w^2*(emiCalc - (mutualCalc^2/recCalc)));               % REF2, page 27, table II, SP topology.
    
elseif type == 3
    
    emiCapacitance = recCalc*recCapacitance/(emiCalc + (mutualCalc^4/(emiCalc*recCalc*recCapacitance*(real(Load))^2))); % REF1, table II, PS topology, but correcting R to R^2.
    emiCapacitance2 = emiCalc/((w^2*mutualCalc^2/real(Load))^2 + w^2*emiCalc^2);                                        % REF2, page 27, table II, PS topology.
    
elseif type == 4
    
    emiCapacitance = (recCalc^2*recCapacitance*(emiCalc*recCalc - mutualCalc^2))/(mutualCalc^4*recCapacitance*(real(Load))^2/recCalc + (emiCalc*recCalc - mutualCalc^2)^2); % REF1, table II, PP topology, but correcting R to R^2.
    emiCapacitance2 = (emiCalc - (mutualCalc^2/recCalc))/(((mutualCalc^2*real(Load))/recCalc^2)^2 + w^2*(emiCalc - (mutualCalc^2/recCalc))^2); % REF2, page 27, table II, PP topology.
end

disp('Receiver Capacitance:');
disp(recCapacitance);

disp('Emitter Capacitance:');
disp(emiCapacitance);
disp(emiCapacitance2);