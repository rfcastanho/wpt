%##########################################################################
%File: mag_eqv.m                                                          #
%Created by: Rodolfo Castanho Fernandes  (rodolfocf@ieee.org)             #
%Last Update: March 29, 2014                                              #
%Description: Main routine.                                               #
%##########################################################################

clearvars -except Mprt   %Mprt is a matrix that stores results of mutiple simulations.
clc                      %If necessary, clear Mprt using "clear all" in command line.
%% Version Control
ver = '0.1';
lastUpdate = '29 Mar 2014';

%% Initial Config
global Mprt
tst = exist('Mprt');    %Verify if Mprt already exist, If not, set it as an empty matrix.

if tst ==0
    Mprt = [];
end

planarPossible = 0;     %Here, initial control flags are defined.
proceedMutual = 0;
proceedMisalign = 0;
logDesired = 0;
str2 = '0';
nameCompare = 1;
%% Project Name

disp('*****************');
disp('This step requires basic project specifications.');
disp(' ');
disp('If a project name is not entered below, some results will not be made available. Hence, it is strongly recommended that');
disp('a project name is used.');
projectName = input('Enter project name (type "0" if a report is not desired): ', 's'); %Here the user can define a project name.
nameCompare = strcmp(projectName,str2);                 %nameCompare = 0 if user ented a valid name in the line above.
                                                        %If projectName = "0", the report is not created.    
if nameCompare == 0
    fileName = sprintf('results_%s.txt',projectName);   %Create a txt file with chosen project name.
    datafile = fopen(fileName,'wt');                    %Prepare txt file to be edited.
    fprintf(datafile,'##########################################################################');      %This is the basic header.
    fprintf(datafile,'\n# Iterative Method for Inductive Power Transfer Coils Representation     #');    %This is the basic header.
    fprintf(datafile,'\n# Algorithm Version: %s                                                 #',ver); %This is the basic header.
    fprintf(datafile,'\n# Last Update: %s                                               #',lastUpdate);  %This is the basic header.
    fprintf(datafile,'\n##########################################################################\n');  %This is the basic header.
    fclose(datafile);
    logDesired = 1;  
else
    fileName = 0;
    logDesired = 0;
end

%% Input Information
% emiInduct = 180;          %While developing the algorithms, this commented lines can be used to enter a set of test parameters.
% recInduct = 115;
% recOpenVoltage = 100;
% emiOpenCurrent = 15;
% recShortCurrent = 3;
% frequency = 30;

emiInduct =    input('Target emitter self-inductance(uH): ');                           %Emitter Self-Inductance.
recInduct =    input('Target receiver self-inductance(uH): ');                          %Receiver Self-Inductance.
recOpenVoltage =   input('Receiver open circuit voltage(peak): ');                      %Receiver open circuit voltage, peak value.
emiOpenCurrent =   input('Emitter current with open receiver circuit(peak):');          %Emmiter current when receiver circuit is open, peak value.
recShortCurrent =   input('Maximum expected receiver current(peak):');                  %Maximum current measured at the receiver circuit, peak value. It can be the short circuit current.
frequency =     input('Operating frequency (kHz):');                                    %Operating frequency. Usually this is the resonant circuit frequency.

if logDesired == 1                      %If report was requested, the above inputs are saved in .txt file.
    datafile = fopen(fileName,'at'); 
    fprintf(datafile,'\n ****** Input Parameters ******');
    fprintf(datafile,'\n Target emitter self-inductance: %d uH',emiInduct);
    fprintf(datafile,'\n Target receiver self-inductance: %d uH',recInduct);
    fprintf(datafile,'\n Receiver open circuit voltage: %d V peak',recOpenVoltage);
    fprintf(datafile,'\n Emitter current with open receiver circuit: %d A peak',emiOpenCurrent);
    fprintf(datafile,'\n Maximum expected receiver current: %d A peak',recShortCurrent);
    fprintf(datafile,'\n Operating frequency: %d kHz',frequency);
    fclose(datafile);
end


%% Main Program

mainTimer = tic;                                            %Internal timer, not used currently.
emiInduct = emiInduct*1e-6; recInduct = recInduct*1e-6;     %Convert input inductances to uH.
frequency = frequency*1000;                                 %Convert input frequency to kHz.

clc;
disp('*****************');
disp('This step will design electrical conductors based on frequency and emitter and receiver currents.');
disp(' ');

wind = winding(emiOpenCurrent,recShortCurrent,frequency);   %Call routine to design electrical conductors, 
                                                            %based on input frequency and maximum expected currents.
x1 = wind(1);
x2 = wind(2);
y1 = wind(3);                                                                   
y2 = wind(4);
y1_limite_inf = wind(5);
y2_limite_inf = wind(6);
xfill2 = wind(7);
yfill2 = wind(8);
xfill1 = wind(9);
yfill1 = wind(10);

y2_limite_inf_mi = y2_limite_inf;   %Save a copy of y2_limite_inf, to use in axial.m if coupling.m is previously executed.

clc;
disp('*****************');
disp('Electrical conductors for both emitter and receiver were determined.');
disp('In this step, coil form (shape) and mechanical restrictions must be defined.');
disp(' ');
form = input('Type "0" for concentrated circular windings or "1" for spiral winding:');          %Here, user decides if windings are concentrated (circular) or spiral.
restriction = input('Type "0" for max. number of turns or "1" for max. diameter restriction:');  %Then, user decides what kind of mechanical restriction will be used (max number of turns
                                                                                                %or max diameter).

%      Table of Results:    
%  Form          Restriction
%   0               0           = Circular with fixed number of turns (diameter will be adjusted)
%   0               1           = Circular with fixed diameter (number of turns will be adjusted)
%   1               0           = Spiral with fixed number of turns(diameter will be adjusted)
%   1               1           = Spiral with fixed diameter (number of turns will be adjusted)

choiceVector = [form,restriction];  %This vector represents the choice.

if choiceVector == [0,0]
      
    planarPossible = 1;         %If "Circular with fixed number of turns" was chosen, user must inform number of turns for emitter and receiver.
    emiTurn =    input('Emitter number of turns: ');
    recTurn =    input('Receiver number of turns: ');
    conc_nfix = confix(fileName,logDesired,frequency,emiInduct,recInduct,emiTurn,recTurn,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent);                             
    
    y1_limite_sup = conc_nfix(1); %Outputs of the above function.
    emiCalc = conc_nfix(2);
    recCalc = conc_nfix(3);
    emiRadius = conc_nfix(4);
    recRadius = conc_nfix(5);
    usingShield = conc_nfix(6);

elseif choiceVector == [0,1]
        
    planarPossible = 0;         %If "Circular with fixed diameter" was chosen, user must inform max diameters for emitter and receiver.
    emiMaxDiameter =    input('Emitter maximum diameter (Dp), in mm: ');            
    recMaxDiameter =    input('Receiver maximum diameter (Ds), in mm: ');          
    conc_dfix = codfix(fileName,logDesired,frequency,emiMaxDiameter,recMaxDiameter,emiInduct,recInduct,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent);

    y1_limite_sup = conc_dfix(1);   %Outputs of the above function.
    emiCalc = conc_dfix(2);
    recCalc = conc_dfix(3);
    emiRadius = conc_dfix(4);
    recRadius = conc_dfix(5); 
    emiTurn = conc_dfix(6);
    recTurn = conc_dfix(7);
    usingShield = conc_dfix(8);
        
elseif choiceVector == [1,0]
    
    planarPossible = 0;         %If "Spiral with fixed number of turns" was chosen, user must inform number of turns for emitter and receiver.
    emiTurn =    input('Emitter number of turns: ');
    recTurn =    input('Receiver number of turns: ');
    panc_nfix = panfix(fileName,logDesired,frequency,emiInduct,recInduct,emiTurn,recTurn,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent);    
        
    y1_limite_sup = panc_nfix(1);   %Outputs of the above function.
    emiCalc = panc_nfix(2);
    recCalc = panc_nfix(3);
    % emiRadius = panc_nfix(4);
    % recRadius = panc_nfix(5);
    usingShield = panc_nfix(6);

elseif choiceVector == [1,1]
        
    planarPossible = 0;         %If "Spiral with fixed diameter" was chosen, user must inform max diameters for emitter and receiver.
    emiMaxDiameter =    input('Emitter maximum diameter (Dp), in mm: ');            
    recMaxDiameter =    input('Receiver maximum diameter (Ds), in mm: ');    
    panc_dfix = padfix(fileName,logDesired,frequency,emiMaxDiameter,recMaxDiameter,emiInduct,recInduct,x1,x2,y1,y2,y1_limite_inf,y2_limite_inf,xfill2,yfill2,xfill1,yfill1,emiOpenCurrent,recShortCurrent);

    y1_limite_sup = panc_dfix(1);   %Outputs of the above function.
    emiCalc = panc_dfix(2);
    recCalc = panc_dfix(3);
    emiRadius = panc_dfix(4);
    recRadius = panc_dfix(5);
    usingShield = panc_dfix(6);
    emiTurn = panc_dfix(7);
    recTurn = panc_dfix(8);
end

disp(' ');
disp('Equivalent geometry obtained. Type any key to continue.');
pause;

clc;        %When routine reaches this point, the input parameters are already represented according to the desired coil form.
disp('*****************');
disp('Based on input parameters and desired coil form, a equivalent geometric representation');
disp('is already available in the file "result_dimensions.fem". This geometry will be analyzed in the next steps.');
disp('Note 1: Resulting graphs will be displayed, but text results are shown in the text report only.');
disp('Note 2: Save graphs for further reference, if necessary. Graphs are not saved automatically');
disp(' ');
disp('Initially, the equivalent axial separation distance (e) can be determined based on input parameters.');
disp('Parameter (e) is the distance between emitter and receiver (the gap) for which the Receiver Open Voltage occurs when the');
disp('Emitter Current is applied at the chosen Frequency. In other words, (e) is the gap that corresponds to the desired coupling factor');
disp('or mutual inductance.');

disp(' ');
mutualDesired = input ('Type "1" to perform coupling analysis or "0" to skip this step:');

if mutualDesired == 1
    mutual = coupling(logDesired,fileName,recOpenVoltage,emiOpenCurrent,frequency,y2_limite_inf,y1_limite_sup,yfill2,emiCalc,recCalc,usingShield);               %Obtenção do fator de acoplamento correto  
    proceedMutual = mutual(1);   %Outputs of the above function.
    yfill2 = mutual(2);
    dist = mutual(3);
    y2_limite_inf = mutual(4);
end

clc;
disp('*****************');
disp('This step is able to perform axial misalignment analysis for the obtained coil configuration.');
axialDesired = input ('Type "1" to perform axial misalignment analysis or "0" to skip this step:');

if axialDesired == 1
    axialMis = axial(fileName,logDesired,emiOpenCurrent,frequency,y2_limite_inf_mi,y1_limite_sup,emiCalc,recCalc,usingShield,emiTurn,recTurn);
    axxb = axialMis(1);         %Outputs of the above function.
end

clc;

if (proceedMutual == 1) && (planarPossible == 1) && (emiTurn == recTurn) && (emiInduct == recInduct)
    
    disp('*****************');
    disp('This is a special coil where emitter and receiver are equal (geometrically and electrically).');
    disp('Therefore, a planar representation is possible. If this is step is performed, it will be possible to');
    disp('analyze emitter and receiver in terms of lateral and angular misalignment.');
    disp('Note 1: This is an aproximation when compared to the axissymetric simulation. In most cases, number of turns');
    disp('and self-inductances will slightly differ from those obtained in previous calculations. User discretion must be'); 
    disp('employed to decide if the resulting aproximation error is not acceptable or not.');
    disp(' ');
    planarDesired = input('Type "1" to perform planar representation or "0" to skip this step:');
    
    if planarDesired == 1  
        planar = rep_planar(fileName,logDesired,frequency,y1,y2_limite_inf,emiRadius,emiTurn,emiInduct,dist); %Chama função de representação planar.
        proceedMisalign = planar(1);     %Outputs of the above function.
        emiCalc = planar(2);
    end
    
    if proceedMisalign == 1
        
       clc; 
       disp('*****************');
       disp('Planar representation was performed.');
       disp(' ');
       misalignDesired = input('Type (1) to perform misalignment analyses or "0" to skip this step:');
       
       if misalignDesired == 1
           misAl = misalign(frequency,emiCalc,emiOpenCurrent,yfill2,y2_limite_inf,dist,emiRadius);
           xxa = misAl(1);       %Outputs of the above function.
       end
    end
%       
% else
%     disp('Nota: Representação planar não é possível.');

end

clc;
disp ('***** End of geometric representation *****');

executionTime = toc(mainTimer);

% This is the end of the routine.
