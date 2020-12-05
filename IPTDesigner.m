%##########################################################################
%File: IPTDesigner.m                                                      #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: August 23, 2014                                             #
%Description: Designs DD, DDQ, BP and Circular IPT Systems easily !       #
%##########################################################################
clear; clc;

%% Initial Definitions and Constraints

primConcluded = 0;                                                         %Primary design not concluded;
n = 0;                                                                     %Generic counter;
Np = 1; Ns = 1;                                                            %Primary and secondary initial number of turns.
u0 = 4*pi*1e-7;                                                            %Magnetic permeability [H/m].
urcopper = 1;                                                              %Copper relative magnetic permeability.
ro = 1.72e-8;                                                              %Copper resistivity [Ohm.m].
J = 3.5;                                                                   %Current density [A/mm2].
matLayerDist = 1;                                                          %Distance between different material layers [mm];
iterP = 0; iterS = 0; %Generic iteration counters.

projectName = input('Enter project name:','s');                            %Here the user can define a project name.

targetSu   = input('Enter target peak uncompensated power, Su (VA): ');
targetVoc  = input('Enter target secondary peak open circuit voltage, Voc (V): ');
targetk    = input('Enter target coupling coefficient (enter "0" if any coupling is acceptable), k: ');
disp(' ');
Ip = input('Enter emitter peak current, Ip (A): ');
emiVin = input('Enter emitter First Harmonic peak Voltage, Vin,FHA (V): ');
f0  = input('Enter resonant frequency, f0 (kHz): ');
targetLp   = input('Enter target emitter inductance, Lp (uH): ');
disp(' ');
disp('Possible pad shapes are Circular (C),Bipolar (BP) or Double-D(DD).');
primType   = input('Enter primary type: ', 's');
secType    = input('Enter secondary type: ', 's');
disp(' ');
disp('Below, enter mechanical restrictions according Possible Geometries.pdf file');
disp('Primary configuration:');
primDimension = input('Enter primary dimension, wp (mm): ');
primAlThickness = input('Enter primary Aluminium thickness, Al (mm): ');
primFeThickness = input('Enter primary Ferrite thickness, Fe (mm): ');
secDimension  = input('Enter secondary dimension, ws (mm): ');
disp('Secondary configuration:');
secAlThickness = input('Enter secondary Aluminium thickness, Al (mm): ');
secFeThickness = input('Enter secondary Ferrite thickness, Fe (mm): ');
depth            = input('Enter simulation depth, de (mm): ');
disp(' ');
disp('Now, define primary-secondary relative position.');
horizontalDist   = input('Enter horizontal distance, dc (mm): ');
verticalDist     = input('Enter vertical distance, e (mm): ');

% primType = 'C';
% secType = 'C';
% targetSu = 400; targetVoc = 60; Ip = 21.2;  %Lixo, apagar.
% targetIsc = 6*sqrt(2);
% emiVin = 250; f0 = 38.4; targetLp = 50;
% primAlThickness = 2; primFeThickness = 3; secAlThickness = 2;
% secFeThickness = 3; depth = 200; secDimension = 150;
% primDimension = 200; horizontalDist = 30; verticalDist = 20;  %Lixo, apagar.


f0  = f0*1e3; targetLp = targetLp*1e-6;                       %Adjusting units.
w = 2*pi*f0;
horizontalDist = abs(horizontalDist); verticalDist = abs(verticalDist);
targetIsc  = targetSu/targetVoc;                              %Calculated target secondary peak short-circuit current.

%% Open standard FEM design

basicFile = 'blankDesign';
%basicFile = sprintf('blank_%s-%s',primType,secType);
targetFile = sprintf('%s.fem',basicFile);
targetFile2 = sprintf('/%s.fem',basicFile);

openfemm;                            %Open FEMM
try
    vv=ver;
    opendocument([cd,targetFile2]);  %Open desired file.
catch
    opendocument(targetFile);
end
mi_saveas('temp.fem');        %Save temporary copy.
mi_probdef(f0,'millimeters','planar',1e-8,depth,20,(0));

%% Winding Design

disp(' ');
disp('Do you wish to calculate primary and secondary litz conductors automatically,');
disp('or enter your own conductor diameters ?');
AutoLitz = input('Enter "1" for Automatic Litz or "0" for your own conductors:');

totalClock = tic;

if AutoLitz == 1
    
    delta         = 1e3*sqrt(2*ro/((2*pi*f0)*urcopper*u0));                %Penetretaion depth [mm].
    awgEquivalent = round((log(2*delta/8.2514694))/(-0.115943));           %Equivalent AWG conductor.

    if awgEquivalent < 0
        awgEquivalent = 0;
    end

    awgChosen        = awgEquivalent + 1;                                  %The chosen conductor must be larger than the calculated equivalent.
    primCrossSection = (Ip/sqrt(2))/J;                                     %Primary cross section, calculated with RMS of Ip.
    rcondp           = sqrt(primCrossSection/pi);                          %Radius of necessary primary conductor.
    
    limitedIsc       = targetIsc/2;                                        %Since Isc may not occur in practice, the secondary winding doesn't need to be designed to carry Isc.
                                                                           %Here, Isc/2 is assumed. This can be changed.
    secCrossSection  = (limitedIsc/sqrt(2))/J;                             %Secondary cross section, calculated with RMS of Isc.
    rconds           = sqrt(secCrossSection/pi);                           %Radius of necessary secondary conductor.

    rskin  = (8.2514694*exp(-0.115943*awgChosen))/2;                       %Converting AWG to mm.
                                                               
    nlitzp = ceil((rcondp/rskin)^2);                                       %Primary Litz strands.
    nlitzs = ceil((rconds/rskin)^2);                                       %Secondary Litz strands.
   
    mi_deletematerial('Litz primary');                          
    mi_deletematerial('Litz secondary');
    mi_addmaterial('Litz primary',1,1,0,0,58,0,0,1,5,0,0,nlitzp,2*rskin);  %Create Litz wires with previously calculated specifications.
    mi_addmaterial('Litz secondary',1,1,0,0,58,0,0,1,5,0,0,nlitzs,2*rskin);

    windp = 2*rcondp;  winds = 2*rconds; 
    primWire = 'Litz primary'; secWire = 'Litz secondary';
    
else
    
    primCondDiameter = input('Enter primary conductor diameter (mm):');
    secCondDiameter = input('Enter secondary conductor diameter (mm):');

    mi_deletematerial('Solid primary');                          
    mi_deletematerial('Solid secondary');
    mi_addmaterial('Solid primary',1,1,0,0,58,0,0,1,5,0,0,1,primCondDiameter);  %Create solid wires with previously entered specifications.
    mi_addmaterial('Solid secondary',1,1,0,0,58,0,0,1,5,0,0,1,secCondDiameter); 
    
    windp = primCondDiameter;  winds = secCondDiameter;
    primWire = 'Solid primary'; secWire = 'Solid secondary';
    limitedIsc       = targetIsc/2;  
end


%% Pad Design

%Build primary:
mi_drawrectangle(-primDimension/2,0,primDimension/2,primAlThickness);
mi_drawrectangle(-primDimension/2,primAlThickness,primDimension/2,primAlThickness+matLayerDist);
mi_drawrectangle(-primDimension/2,primAlThickness+matLayerDist,primDimension/2,primAlThickness+matLayerDist+primFeThickness);
mi_drawrectangle(-primDimension/2,primAlThickness+matLayerDist+primFeThickness,primDimension/2,primAlThickness+2*matLayerDist+primFeThickness);
mi_drawrectangle(-primDimension/2,primAlThickness+2*matLayerDist+primFeThickness,-primDimension/2+windp,primAlThickness+2*matLayerDist+primFeThickness+windp);
mi_drawrectangle(primDimension/2,primAlThickness+2*matLayerDist+primFeThickness,primDimension/2-windp,primAlThickness+2*matLayerDist+primFeThickness+windp);

mi_addblocklabel(0,primAlThickness/2); mi_addblocklabel(0,primAlThickness + matLayerDist/2);
mi_addblocklabel(0,primAlThickness + matLayerDist + primFeThickness/2);
mi_addblocklabel(0,primAlThickness + matLayerDist + primFeThickness + matLayerDist/2);
mi_addblocklabel(-primDimension/2 + windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2);
mi_addblocklabel(primDimension/2 - windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2);

mi_addblocklabel(0,-20); mi_selectlabel(0,-20);mi_setblockprop('Air',0,0.5,0,0,0,0); mi_clearselected;

mi_selectlabel(0,primAlThickness/2); mi_setblockprop('Aluminum, 1100',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(0,primAlThickness + matLayerDist/2); mi_setblockprop('Air',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(0,primAlThickness + matLayerDist + primFeThickness/2); mi_setblockprop('Soft magnetic ferrite (Fe-Ni-Zn-V)',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(0,primAlThickness + matLayerDist + primFeThickness + matLayerDist/2); mi_setblockprop('Air',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(-primDimension/2 + windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); mi_setblockprop(primWire,0,0.5,'primary',0,0,-Np); mi_clearselected;
mi_selectlabel(primDimension/2 - windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); mi_setblockprop(primWire,0,0.5,'primary',0,0,Np); mi_clearselected;

mi_addblocklabel(0,-20); mi_selectlabel(0,-20);mi_setblockprop('Air',0,2,0,0,0,0); mi_clearselected;

ax  = primDimension/2-windp;
ay1 = primAlThickness+2*matLayerDist+primFeThickness+windp;
ay2 = primAlThickness+2*matLayerDist+primFeThickness-windp;
mi_selectnode(ax,ay1); mi_selectnode(ax,ay2);mi_setgroup(1);

bx  = -primDimension/2+windp;
by1 = primAlThickness+2*matLayerDist+primFeThickness+windp; primMaxHeight = by1;
by2 = primAlThickness+2*matLayerDist+primFeThickness-windp;
mi_selectnode(bx,by1); mi_selectnode(bx,by2);mi_setgroup(2); mi_clearselected;
sec0 = primMaxHeight + verticalDist;

%Build secondary:
mi_drawrectangle(horizontalDist-secDimension/2,sec0,horizontalDist-secDimension/2+winds,sec0+winds);
mi_drawrectangle(horizontalDist-secDimension/2,sec0+winds,horizontalDist-secDimension/2+secDimension,sec0+winds+matLayerDist);
mi_drawrectangle(horizontalDist-secDimension/2,sec0+winds+matLayerDist,horizontalDist-secDimension/2+secDimension,sec0+winds+matLayerDist+secFeThickness);
mi_drawrectangle(horizontalDist-secDimension/2,sec0+winds+matLayerDist+secFeThickness,horizontalDist-secDimension/2+secDimension,sec0+winds+2*matLayerDist+secFeThickness);
mi_drawrectangle(horizontalDist-secDimension/2,sec0+winds+2*matLayerDist+secFeThickness,horizontalDist+secDimension-secDimension/2,sec0+winds+2*matLayerDist+secFeThickness+secAlThickness);
mi_drawrectangle(horizontalDist-secDimension/2+secDimension-winds,sec0,horizontalDist+secDimension-secDimension/2,sec0+winds);

mi_addblocklabel(horizontalDist,sec0+winds+matLayerDist/2); mi_addblocklabel(horizontalDist,sec0+winds+matLayerDist+primFeThickness/2);
mi_addblocklabel(horizontalDist,sec0+winds+matLayerDist+primFeThickness+matLayerDist/2);
mi_addblocklabel(horizontalDist,sec0+winds+2*matLayerDist+primFeThickness+primAlThickness/2);
mi_addblocklabel(horizontalDist-secDimension/2+winds/2,sec0+winds/2);
mi_addblocklabel(horizontalDist+secDimension/2-winds/2,sec0+winds/2);

mi_selectlabel(horizontalDist,sec0+winds+2*matLayerDist+primFeThickness+primAlThickness/2); mi_setblockprop('Aluminum, 1100',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(horizontalDist,sec0+winds+matLayerDist/2); mi_setblockprop('Air',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(horizontalDist,sec0+winds+matLayerDist+primFeThickness/2); mi_setblockprop('Soft magnetic ferrite (Fe-Ni-Zn-V)',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(horizontalDist,sec0+winds+matLayerDist+primFeThickness+matLayerDist/22); mi_setblockprop('Air',0,0.5,0,0,0,0); mi_clearselected;
mi_selectlabel(horizontalDist-secDimension/2+winds/2,sec0+winds/2); mi_setblockprop(secWire,0,0.5,'secondary',0,0,-Ns); mi_clearselected;
mi_selectlabel(horizontalDist+secDimension/2-winds/2,sec0+winds/2); mi_setblockprop(secWire,0,0.5,'secondary',0,0,Ns); mi_clearselected;

mx  = horizontalDist+secDimension/2-winds;
my1 = sec0;
my2 = sec0+winds;
mi_selectnode(mx,my1); mi_selectnode(mx,my2);mi_setgroup(3);

px  = horizontalDist-secDimension/2+winds;
py1 = sec0;
py2 = sec0+winds;
mi_selectnode(px,py1); mi_selectnode(px,py2);mi_setgroup(4); mi_clearselected;

%% Simulation Domain Definition

apoint = mi_selectnode(10000,10000);    mi_clearselected; ax = apoint(1,1); ay = apoint(1,2);
bpoint = mi_selectnode(10000,-10000);   mi_clearselected; bx = bpoint(1,1); by = bpoint(1,2);
cpoint = mi_selectnode(-10000,-10000);  mi_clearselected; cx = cpoint(1,1); cy = cpoint(1,2);
dpoint = mi_selectnode(-10000,10000);   mi_clearselected; dx = dpoint(1,1); dy = apoint(1,2);

dim = [ax bx cx dx ay by cy dy];
max = abs(max(dim)); 

mi_drawrectangle(-3*max,-1*max,3*max,2*max);

mi_selectsegment(0,10000); mi_selectsegment(0,-10000);
mi_selectsegment(10000,0); mi_selectsegment(-10000,0); 
mi_setsegmentprop('Mixed',0,1,0,0);  mi_clearselected;

%% Parameter Limits

kMin = 0.2;     %Minimum acceptable coupling due to Q = (1/k) <= 5.
kMax = 0.5;     %Maximum acceptable coupling.

calculatedkMin = sqrt(targetSu/(w*targetLp*(Ip/sqrt(2))^2));  %Minimum coupling that satisfies target Su.

if calculatedkMin < kMin
    disp('Target Su can not be achieved !');
else
    kMin = calculatedkMin;
end
    
mi_setcurrent('secondary',0);   
mi_setcurrent('primary',Ip);   
mi_analyze(1);
mi_loadsolution;
r = mo_getcircuitproperties('primary'); LpCalc = abs(r(3)/r(1)); rpCalc = real(r(2)/r(1));
primError = abs((targetLp - LpCalc)/targetLp);

primSidesDist = primDimension-2*Np*windp;
groupA = 1; groupB = 2;
%% Iterations

%Primary design
emiClock = tic;
while (primError > 0.05 && primConcluded == 0)
    
    iterP = iterP + 1;
    lastError = primError;
    Np = Np + 1;
    primSidesDist = primDimension-2*Np*windp;
    
    if primSidesDist <= 0
       disp('Primary was not designed ! Try modifying restrictions.');
       break;
    end
    
    mi_selectgroup(groupA); mi_movetranslate(-windp,0);
    mi_selectgroup(groupB); mi_movetranslate(windp,0);   
    mi_selectlabel(-primDimension/2 + windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); 
    mi_setblockprop(primWire,0,0.5,'primary',0,0,-Np); mi_clearselected;
    mi_selectlabel(primDimension/2 - windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); 
    mi_setblockprop(primWire,0,0.5,'primary',0,0,Np); mi_clearselected;
    
    mi_setcurrent('primary',Ip);   
    mi_setcurrent('secondary',0);    
    mi_analyze(1);
    mi_loadsolution;
    r = mo_getcircuitproperties('primary'); LpCalc = abs(r(3)/r(1)); rpCalc = real(r(2)/r(1)); emiVinCalc = abs(r(2));
    h = mo_getcircuitproperties('secondary');secOpenCalc = abs(h(2)); mutualCalc = abs(secOpenCalc/(w*Ip));

    
    primError = abs((targetLp - LpCalc)/targetLp)
    
    if lastError <= primError
        Np = Np - 1;                            %Undo last iteration.
        primSidesDist = primDimension-2*Np*windp;
        mi_selectgroup(groupA); mi_movetranslate(windp,0);
        mi_selectgroup(groupB); mi_movetranslate(-windp,0);   
        mi_selectlabel(-primDimension/2 + windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); 
        mi_setblockprop(primWire,0,0.5,'primary',0,0,-Np); mi_clearselected;
        mi_selectlabel(primDimension/2 - windp/2,primAlThickness + 2*matLayerDist + primFeThickness + windp/2); 
        mi_setblockprop(primWire,0,0.5,'primary',0,0,Np); mi_clearselected;
        primError = lastError;

        mi_setcurrent('primary',Ip);   
        mi_setcurrent('secondary',0);    
        mi_analyze(1);
        mi_loadsolution;
        r = mo_getcircuitproperties('primary'); LpCalc = abs(r(3)/r(1)); rpCalc = real(r(2)/r(1)); emiVinCalc = abs(r(2));
        h = mo_getcircuitproperties('secondary');secOpenCalc = abs(h(2)); mutualCalc = abs(secOpenCalc/(w*Ip));
        primConcluded = 1;
    end
    
end
emiTime = toc (emiClock); %Total time for emitter calculation.

%Secondary design

secConcluded = 0;
secSidesDist = secDimension-2*Ns*winds;
groupA = 3; groupB = 4;

mi_setcurrent('primary',0);   
mi_setcurrent('secondary',limitedIsc);   
mi_analyze(1);
mi_loadsolution;
r = mo_getcircuitproperties('secondary'); LsCalc = abs(r(3)/r(1)); rsCalc = real(r(2)/r(1));
   
coupling = mutualCalc/(sqrt(LpCalc*LsCalc));
SuCalc = w*(coupling^2)*LpCalc*(Ip/sqrt(2))^2;
VocCalc = secOpenCalc;

recClock = tic;
while (VocCalc < targetVoc || SuCalc < targetSu || IscCalc < targetIsc)    
    
    iterS = iterS + 1;
    lastCoupling = coupling; lastSu = SuCalc; lastVoc = VocCalc;
    Ns = Ns + 1;
    secSidesDist = secDimension-2*Ns*winds;
    
    if secSidesDist <= 0
        disp('Target inductance can not be achieved with specified restrictions.');
        disp('Secondary was not designed ! Try modifying restrictions.');
        break;    
    end
    
    mi_selectgroup(groupA); mi_movetranslate(-winds,0);
    mi_selectgroup(groupB); mi_movetranslate(winds,0);   
    mi_selectlabel(horizontalDist-secDimension/2+winds/2,sec0+winds/2); 
    mi_setblockprop(secWire,0,0.5,'secondary',0,0,-Ns); mi_clearselected;
    mi_selectlabel(horizontalDist+secDimension/2-winds/2,sec0+winds/2); 
    mi_setblockprop(secWire,0,0.5,'secondary',0,0,Ns); mi_clearselected;
        
    mi_setcurrent('primary',0);
    mi_setcurrent('secondary',limitedIsc);   
    mi_analyze(1);
    mi_loadsolution;
    r = mo_getcircuitproperties('secondary'); LsCalc = abs(r(3)/r(1)); rsCalc = real(r(2)/r(1));
    
    mi_setcurrent('primary',Ip);                            %Only emitter is excited.
    mi_setcurrent('secondary',0);
    mi_analyze(1);
    mi_loadsolution;
    r = mo_getcircuitproperties('primary'); LpCalc = abs(r(3)/r(1)); rpCalc = real(r(2)/r(1)); emiVinCalc = abs(r(2));
    h = mo_getcircuitproperties('secondary');secOpenCalc = abs(h(2)); mutualCalc = abs(secOpenCalc/(w*Ip));
    
    coupling = mutualCalc/(sqrt(LpCalc*LsCalc));
    VocCalc = secOpenCalc;
    SuCalc = w*(coupling^2)*LpCalc*(Ip/sqrt(2))^2;
    IscCalc = mutualCalc*(Ip/sqrt(2))/LsCalc;
    
    cV = targetLp*targetVoc^2/(emiVin^2*kMin^2);
    cS = sqrt(targetSu*LsCalc/(w*(Ip/sqrt(2))^2));
    cI = targetIsc*(rsCalc + 1i*w*LsCalc)/(1i*w*Ip);
    
    disp(' ');
    disp(sprintf('Lp: %d',LpCalc));
    disp(sprintf('Ls: %d',LsCalc));
    disp(sprintf('k: %d',coupling));
    disp(sprintf('cV: %d',cV));
    disp(sprintf('VocCalc: %d',VocCalc));
    disp(sprintf('M: %d',mutualCalc));
    disp(sprintf('cS: %d',cS));
    disp(sprintf('Su (VA): %d',SuCalc));
    disp(sprintf('cI: %d',cI));
    disp(sprintf('Isc (ARMS): %d',IscCalc));
    
    if (lastCoupling > coupling || lastSu > SuCalc || lastVoc > VocCalc)
        disp('Target inductance can not be achieved with specified restrictions.');
        break
    end
    
end
recTime = toc (recClock); %Total time for receiver calculation.

totalTime = toc (totalClock); %Total algorithm run time.

j = 1;


