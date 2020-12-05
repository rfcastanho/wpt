% FEM-Integrated Sensitivity Analysis
clear; clc;

clockSens = tic;
detail = 5; 
dev = 1/100; %Acceptable deviation for obtained inductances in respect to nominal (reference) values.
resultsAgree = 0;
updatePosition = 0;
satDetect = 0;
loadCondition = 1;
n = 0;
contfig = 1;

geomName = input('Enter file name (without extension): ', 's');      %Here, user defines a previous simulated geometry to be opened.
targetFile = sprintf('%s.fem',geomName);
targetFile2 = sprintf('/%s.fem',geomName);

openfemm;                            %Open FEMM
try
    vv=ver;
    opendocument([cd,targetFile2]);  %Open desired file.
catch
    opendocument(targetFile);
end
mi_saveas('temp_target.fem');        %Save temporary copy.

mi_analyze(1);
mi_loadsolution;
v = mo_getprobleminfo;
simType = v(1);
depth = v(3);

if simType == 0
    simType = 'planar';
elseif simType == 1
    simType = 'axi';
end

%% Define Electric Circuit Compensation Topology

disp('From the desired file to be analyzed, and considering that the loaded file represents');
disp('the magnetic stage in its NOMINAL OPERATING CONDITIONS, define:');

% emiCurrentRef = input('Emitter nominal peak current (A):'); %These are the expected nominal peak currents.
% recCurrentRef = input('Receiver nominal peak current (A):');
% frequencyRef =  input('System nominal frequency (kHz):');

emiCurrentRef = 10; recCurrentRef = 20; frequencyRef = 30;  %deletar isto !!

frequencyRef = frequencyRef*1e3; w = 2*pi*frequencyRef;
emiCurrent = emiCurrentRef; recCurrent = recCurrentRef; frequency = frequencyRef;

mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));  %Confirm nominal excitation.
mi_setcurrent('emitter',emiCurrent);                            %Only emitter is excited.
mi_setcurrent('receiver',0);
mi_analyze(1);
mi_loadsolution;
r = mo_getcircuitproperties('emitter'); emiCalcRef = abs(r(3)/r(1)); emiResistanceRef = real(r(2)/r(1));           %Here, emiCalc is obtained "as is", considering 
h = mo_getcircuitproperties('receiver');emiOpenCalc = abs(h(2)); mutualCalcRef = abs(emiOpenCalc/(w*emiCurrent));  %interactions with of nearby elements or and magnetic saturation.

mi_setcurrent('emitter',0);
mi_setcurrent('receiver',-recCurrent);                          %Now, only receiver is excited.
mi_analyze(1);
mi_loadsolution;
v = mo_getcircuitproperties('receiver'); recCalcRef = abs(v(3)/v(1)); recResistanceRef = real(v(2)/v(1));  %Also, recCalc is obtained "as is". 
couplingRef = mutualCalcRef/(sqrt(emiCalcRef*recCalcRef));

% disp('From the desired electric circuit to be considered, define:');
% 
% emiType = input('Type "1" if emitter is series compensated or "0" if it is parallel compensated:');
% 
% recType = input('Type "1" if receiver is series compensated or "0" if it is parallel compensated:');
% 
% emiVoltage = input('Emitter peak input voltage (First Harmonic Aproximation) (V):');
% emiVoltage = emiVoltage + 1i*0; emiVoltageRef = emiVoltage;
%
% rint = input('Source internal resistance (Ohms):'); 
% 
% loadReal = input('Real part of the load:');
% loadImag = input('Imaginary part of the load:'); 
% Load = loadReal + 1i*loadImag; loadRef = Load;

%##########################################################################
emiType = 1; recType = 1; %deletar isto !!
emiVoltage = 200; emiVoltageRef = emiVoltage; rint = 4.3; loadReal = 5;   %deletar isto !!
loadImag = 0; Load = loadReal + 1i*loadImag; loadRef = Load;              %deletar isto !!
%##########################################################################

if emiType == 1;
    emiComp = 'Series';
elseif emiType == 0;
    emiComp = 'Parallel';
end

if recType == 1;
    recComp = 'Series';
elseif recType == 0;
    recComp = 'Parallel';
end

folder = sprintf('%s %s %s Results',geomName,emiComp,recComp);
mkdir(folder);      %Create folder to save results.

%% Calculate/Define Compensation Capacitances
clc;
disp('Compensation capacitances can be automatically calculated to compensate reactances of the loaded geometry.')
autoCap = input('Enter "1" to calculate capacitances automatically or "0" to enter values manually:');

if autoCap == 1
    chosenCompCalc = 'automatically';
    recCalc = recCalcRef; emiCalc = emiCalcRef; mutualCalc = mutualCalcRef; %Capacitance is calculated based on nominal inductances
    capFunction = emicap (emiType,recType,recCalc,emiCalc,mutualCalc,Load,frequency);
    emiCapacitance = capFunction(1);
    recCapacitance = capFunction(2);
else
    emiCapacitance = input('Enter emitter capacitance (uF):');
    recCapacitance = input('Enter receiver capacitance (uF):');
    emiCapacitance = emiCapacitance*1e-6; recCapacitance = recCapacitance*1e-6;
    chosenCompCalc = 'manually';
end
    
%% B-Field Analysis
clc;
disp('B-Field analysis will excite emitter and receiver simultaneously with respective calculated currents.');
disp('This is used to calculate B (magnetic flux density) at a chosen reference point.');
disp('If B-Field Analysis is performed, sensibility analysis will take longer to complete.');
BAnalysis = input('Enter "1" to perform B-Field analysis or "0" to ignore:');

if BAnalysis == 1
    disp('For the reference point where B will be calculated, enter:');
    coordX = input('x coordinate for the reference point:');
    coordY = input('y coordinate for the reference point:');
    BDecision = sprintf('performed for %d (x coordinate) and %d (y coordinate).',coordX,coordY);
%     mi_addnode(coordX,coordY);
else
    BDecision = 'NOT performed';
end

%% Save inputs to log file

fileName = sprintf('Results log for %s.txt',geomName);   %Create a txt file with chosen project name.
datafile = fopen(fileName,'wt');                    %Prepare txt file to be edited.
fprintf(datafile,'##########################################################################');      %This is the basic header.
fprintf(datafile,'\n#                     Sensibility Analysis                               #');    %This is the basic header.
fprintf(datafile,'\n                                                                         #'); %This is the basic header.
fprintf(datafile,'\n ****** Emiter Input Parameters ******');
fprintf(datafile,'\n Nominal emitter self-inductance: %d H',emiCalcRef);
fprintf(datafile,'\n Nominal emitter resistance: %d Ohm',emiResistanceRef);
fprintf(datafile,'\n Emitter compensation: %s ',emiComp);
fprintf(datafile,'\n Emitter compensation capacitance: %d F',emiCapacitance);
fprintf(datafile,'\n');
fprintf(datafile,'\n ****** Receiver Input Parameters ******');
fprintf(datafile,'\n Nominal receiver self-inductance: %d H',recCalcRef);
fprintf(datafile,'\n Nominal receiver resistance: %d Ohm',recResistanceRef);
fprintf(datafile,'\n Receiver compensation: %s ',recComp);
fprintf(datafile,'\n Receiver compensation capacitance: %d F',recCapacitance);
fprintf(datafile,'\n');
fprintf(datafile,'\n *** NOTE: Compensation capacitances were defined %s',chosenCompCalc);
fprintf(datafile,'\n');
fprintf(datafile,'\n Nominal mutual inductance: %d H',mutualCalcRef);
fprintf(datafile,'\n Nominal coupling coefficient: %d ',couplingRef);
fprintf(datafile,'\n');
fprintf(datafile,'\n ****** Nominal Excitations and Load Parameters ******');
fprintf(datafile,'\n Emitter nominal current: %d A peak',emiCurrentRef);
fprintf(datafile,'\n Emitter input voltage: %d V peak',emiVoltageRef);
fprintf(datafile,'\n Emitter input voltage source internal resistance: %d Ohms',rint);
fprintf(datafile,'\n Receiver nominal current: %d A peak',recCurrentRef);
fprintf(datafile,'\n Nominal frequency: %d Hz',frequencyRef);
fprintf(datafile,'\n Real part of nominal load: %d Ohms',loadReal);
fprintf(datafile,'\n Imaginary part of nominal load: %d Ohms',loadImag);

fprintf(datafile,'\n Loaded analysis was %s',BDecision);
fclose(datafile);

%% Frequency Sensibility Analysis
clc;
freqAnalysis = input('Enter "1" to perform frequency sensibility analysis or "0" to ignore it:');

if freqAnalysis == 1
     
    execFrequency = 'Frequency sensibility analysis performed.';
    freqVar= input('Enter frequency variation (eg., 10, for +/- 10% of nominal value):');
    freqMinimum = frequencyRef*(1 - freqVar/100);
    freqMaximum = frequencyRef*(1 + freqVar/100);

    freqRange = freqMaximum - freqMinimum;
    freqStep = (freqRange/2)/(ceil(detail/2));

    disp('During frequency sensibility analysis it is possible to perform freqyency sweep for different load conditions.');
    disp('This is particularly useful to investigate occurance of bifurcation phenomena.');
    disp('Below, enter minimum and maximum expected loads to perform analysis or enter both "0"');
    disp('to perform frequency analysis for nominal load only.');
    
    realLoadMin = input('Enter real part of minimum load condition (ohm):');
    imagLoadMin = input('Enter imaginary part of minimum load condition (ohm):');
    loadMin = realLoadMin + 1i*imagLoadMin;
    realLoadMax = input('Enter real part of maximum load condition (ohm):');
    imagLoadMax = input('Enter imaginary part of maximum load condition (ohm):');
    loadMax = realLoadMax + 1i*imagLoadMax;
    
    if abs(loadMax) == 0
        loadCondition = 2;
    end
    
    if abs(loadMin) == 0
        loadCondition = 3;
    end
    
    while loadCondition < 4
    
    emiCalc = emiCalcRef; recCalc = recCalcRef; mutualCalc = mutualCalcRef; emiResistance = emiResistanceRef; recResistance = recResistanceRef;
    frequency = frequencyRef;
        
    if loadCondition == 3
        selectLoad = 'Nominal Load';
        Load = loadRef;
    elseif loadCondition == 2
        selectLoad = 'Minimum Load';
        Load = loadMin;
    elseif loadCondition == 1;
        selectLoad = 'Maximum Load';
        Load = loadMax;
    end
    
    subfolder = sprintf('Frequency Analysis for %s', selectLoad);
    mkdir(subfolder);
    
        for frequency = freqMinimum:freqStep:freqMaximum            
        
            n = n + 1;
            w = 2*pi*frequency;
            freq(n) = frequency;
            mode = 'normal';
            resultsAgree = 0;
            emiCalcPrev = emiCalcRef; recCalcPrev = recCalcRef; mutualCalcPrev = mutualCalcRef; emiResistancePrev = emiResistanceRef; recResistancePrev = recResistanceRef;
            emiCalc = emiCalcPrev; recCalc = recCalcPrev; mutualCalc = mutualCalcPrev; emiResistance = emiResistancePrev; recResistance = recResistancePrev;
    
            while resultsAgree == 0
        
                modelTransf = tmodel (rint,emiType,recType,emiCalc,mutualCalc,emiCapacitance,w,emiResistance,recResistance,recCalc,recCapacitance,Load,emiVoltage,emiCurrent,mode,recCurrent); %Call T-Model function.
                emiCurrent = modelTransf(1);    %Outputs of the T-Model function are emitter current, receiver current and load current;
                recCurrent = modelTransf(2);
                loadCurrent = modelTransf(3);
    
                mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));  %Update excitation.
                mi_setcurrent('emitter',abs(emiCurrent));                       %Only emitter is excited.
                mi_setcurrent('receiver',0);
                mi_analyze(1);
                mi_loadsolution;
                r = mo_getcircuitproperties('emitter'); emiCalc = abs(r(3)/r(1)); emiResistance = real(r(2)/r(1));              %Here, emiCalc is obtained "as is", considering 
                h = mo_getcircuitproperties('receiver');emiOpenCalc = abs(h(2)); mutualCalc = abs(emiOpenCalc/(w*emiCurrent));  %interactions with of nearby elements or and magnetic saturation.
    
                mi_setcurrent('emitter',0);
                mi_setcurrent('receiver',abs(recCurrent));                       %Now, only receiver is excited.
                mi_analyze(1);
                mi_loadsolution;
                v = mo_getcircuitproperties('receiver'); recCalc = abs(v(3)/v(1)); recResistance = real(v(2)/v(1));  %Also, recCalc is obtained "as is". 

                emiDelta = abs(emiCalc - emiCalcPrev)/emiCalcPrev; 
                recDelta = abs(recCalc - recCalcPrev)/recCalcPrev; 
                mutualDelta = abs(mutualCalc - mutualCalcPrev)/mutualCalcPrev;
                emiResistanceDelta = abs(emiResistance - emiResistancePrev)/emiResistancePrev; recResistanceDelta = abs(recResistance - recResistancePrev)/recResistancePrev;
        
                emiCalcPrev = emiCalc; recCalcPrev = recCalc; mutualCalcPrev = mutualCalc; emiResistancePrev = emiResistance; recResistancePrev = recResistance;
        
                if (emiDelta < dev) &&  (recDelta < dev) && (mutualDelta < dev) && (emiResistanceDelta < dev) && (recResistanceDelta < dev)
                    resultsAgree = 1;    
                end
            end
        
            loadCur(n) = loadCurrent; 
            emiCur(n) = emiCurrent;
            recCur(n) = recCurrent;
        
            inputsCalc = inputPar(emiType,emiVoltage,emiCurrent,rint,emiCapacitance,w);
            sourceCurrentAngle(n) = inputsCalc(1);
        
            SL(n) = Load*loadCurrent*conj(loadCurrent);
            emiQFactor(n) = w*emiCalc/emiResistance;
            recQFactor(n) = w*recCalc/recResistance;
        
            if BAnalysis == 1
               Bfunction = bcalc(emiCurrent,recCurrent,coordX,coordY); 
               rmsB(n) = Bfunction(1);
               
%                if frequency == frequencyRef
%                 rmsBRef = rmsB(n);                              %Save reference value for normalization purposes.
%                end
               
            end
            
            if frequency == frequencyRef                      
                loadCurRef = loadCurrent;                       %Save reference values for normalization purposes.
                emiQRef = emiQFactor(n);
                recQRef = recQFactor(n);
                SLRef = SL(n);
            end
            
        end
        
        if BAnalysis == 1
            figure(contfig);
            plot(freq/frequencyRef,rmsB/1e-6);        %Plot active load power vs. frequency.
            xlabel('Frequency (kHz)'); ylabel('RMS B-Field (uT)');
            title('Frequency Sensibility Analysis');
            saveas(figure(contfig),'rmsB_Freq.fig')
            close(gcf);
            contfig = contfig + 1;
        end
    
        figure(contfig);
        plot(freq/frequencyRef,sourceCurrentAngle);        %Plot active load power vs. frequency.
        xlabel('Frequency (kHz)'); ylabel('Source Current Angle (deg)');
        title('Frequency Sensibility Analysis');
        saveas(figure(contfig),'sourceCurrent_Freq.fig')
        close(gcf);
        contfig = contfig + 1;
    
        figure(contfig);
        plot(freq/frequencyRef,real(SL)/(real(SLRef)));        %Plot active load power vs. frequency.
        xlabel('Frequency (kHz)'); ylabel('Active Power (W)');
        title('Frequency Sensibility Analysis');
        saveas(figure(contfig),'SL_Freq.fig')
        close(gcf);
        contfig = contfig + 1;
    
        figure(contfig);
        plot(freq/frequencyRef,emiQFactor/emiQRef);        %Plot emitter unloaded quality factor vs. frequency.
        xlabel('Frequency (kHz)'); ylabel('Emitter Quality Factor (unloaded)');
        title('Frequency Sensibility Analysis');
        saveas(figure(contfig),'emiQFactor_Freq.fig');
        close(gcf);
        contfig = contfig + 1;
    
        figure(contfig);
        plot(freq/frequencyRef,recQFactor/recQRef);        %Plot emitter unloaded quality factor vs. frequency.
        xlabel('Frequency (kHz)'); ylabel('Receiver Quality Factor (unloaded)');
        title('Frequency Sensibility Analysis');
        saveas(figure(contfig),'recQFactor_Freq.fig');
        close(gcf);
        contfig = contfig + 1;
    
        save('SL.mat','SL');save('freq.mat','freq'); save('emiQFactor.mat','emiQFactor'); save('loadCur.mat','loadCur'); save('emiCur.mat','emiCur');
        save('recCur.mat','recCur'); save('recQFactor.mat','recQFactor');  save('sourceCurrentAngle.mat','sourceCurrentAngle');
        movefile('SL.mat',subfolder); movefile('freq.mat',subfolder); movefile('emiQFactor.mat',subfolder); movefile('recQFactor.mat',subfolder);
        movefile('SL_Freq.fig',subfolder); movefile('emiQFactor_Freq.fig',subfolder); movefile('recQFactor_Freq.fig',subfolder);
        movefile('sourceCurrentAngle.mat',subfolder); movefile('sourceCurrent_Freq.fig',subfolder); movefile('loadCur.mat',subfolder);
        movefile('emiCur.mat',subfolder);movefile('recCur.mat',subfolder);
    
        if BAnalysis == 1
            save('rmsB.mat','rmsB'); movefile('rmsB_Freq.fig',subfolder);  movefile('rmsB.mat',subfolder);
        end
    
        movefile(subfolder,folder);
    
        n = 0;
        SL = 0; emiCalc = 0; recCalc = 0; rmsB = 0; emiQFactor = 0; recQFactor = 0;
        loadCur = 0; emiCur = 0; recCur = 0; sourceCurrentAngle = 0;
        loadCondition = loadCondition + 1;
    
    end
else
    execFrequency = 'Frequency sensibility analysis NOT performed.';
end

if freqAnalysis == 1
closefemm;
delete('temp_target.fem');
delete('temp_target.ans');
end

%% Emitter Current Sensibility Analysis
clc;
emiCurrentAnalysis = input('Enter "1" to perform emitter current sensibility analysis or "0" to ignore it:');

if emiCurrentAnalysis == 1
    
    subfolder = 'Emitter Current Analysis';
    mkdir(subfolder);
    
    execEmiCurrent = 'Emitter current sensibility analysis performed.';
    openfemm;                            %Open FEMM
    
    try
        vv=ver;
        opendocument([cd,targetFile2]);  %Open desired file.
    catch
        opendocument(targetFile);
    end
    mi_saveas('temp_target.fem');        %Save temporary copy.
    
    emiCurVar= input('Enter emitter current variation (eg., 10, for +/- 10% of nominal value):');
    emiCurrentMinimum = emiCurrentRef*(1 - emiCurVar/100);
    emiCurrentMaximum = emiCurrentRef*(1 + emiCurVar/100);

    emiCurrentRange = emiCurrentMaximum - emiCurrentMinimum;
    emiCurrentStep = (emiCurrentRange/2)/(ceil(detail/2));

    emiCalc = emiCalcRef; recCalc = recCalcRef; mutualCalc = mutualCalcRef; emiResistance = emiResistanceRef; recResistance = recResistanceRef;
    frequency = frequencyRef; w = 2*pi*frequencyRef; emiVoltage = emiVoltageRef;

    for emiCurrent = emiCurrentMinimum:emiCurrentStep:emiCurrentMaximum

        n = n + 1;
        emiCur(n) = emiCurrent;
        mode = 'emiCurrentVar';         %Force to enter Emitter Current Variation mode in TModel function.
        resultsAgree = 0;
        emiCalcPrev = emiCalcRef; mutualCalcPrev = mutualCalcRef; emiResistancePrev = emiResistanceRef;
        emiCalc = emiCalcPrev; mutualCalc = mutualCalcPrev; emiResistance = emiResistancePrev;
        
        while resultsAgree == 0
                  
            modelTransf = tmodel (rint,emiType,recType,emiCalc,mutualCalc,emiCapacitance,w,emiResistance,recResistance,recCalc,recCapacitance,Load,emiVoltage,emiCurrent,mode,recCurrent); %Call T-Model function.   
            recCurrent = modelTransf(2); %Outputs of the T-Model function are emitter current (already known), receiver current and load current;
            loadCurrent = modelTransf(3);
            
            mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));  %Update excitation.
            mi_setcurrent('emitter',abs(emiCurrent));                       %Only emitter is excited.
            mi_setcurrent('receiver',0);
            mi_analyze(1);
            mi_loadsolution;
            r = mo_getcircuitproperties('emitter'); emiCalc = abs(r(3)/r(1)); emiResistance = real(r(2)/r(1));              %Here, emiCalc is obtained "as is", considering 
            h = mo_getcircuitproperties('receiver');emiOpenCalc = abs(h(2)); mutualCalc = abs(emiOpenCalc/(w*emiCurrent));  %interactions with of nearby elements or and magnetic saturation.

            %NOTE: At this point, receiver parameter would be re-calculated. However, it is assumed that receiver inductance
            %and resistance are fixed. This means that receiver is not subject to magnetic saturation.
            
            emiDelta = abs(emiCalc - emiCalcPrev)/emiCalcPrev;
            mutualDelta = abs(mutualCalc - mutualCalcPrev)/mutualCalcPrev;
            emiResistanceDelta = abs(emiResistance - emiResistancePrev)/emiResistancePrev;
         
            emiCalcPrev = emiCalc; mutualCalcPrev = mutualCalc; emiResistancePrev = emiResistance;
            
            if ((abs(emiCalc - emiCalcRef)/emiCalcRef) > 0.2)
                clc;
                disp('Possible magnetic saturation detected in emitter !');
            end

            if (emiDelta < dev) && (mutualDelta < dev) && (emiResistanceDelta < dev)
                resultsAgree = 1;
            end
        end
        
        loadCur(n) = loadCurrent; 
        recCur(n) = recCurrent;
        
        SL(n) = Load*loadCurrent*conj(loadCurrent);
        emiCalc(n) = emiCalc;
        
        if BAnalysis == 1
           Bfunction = bcalc(emiCurrent,recCurrent,coordX,coordY); 
           rmsB(n) = Bfunction(1);
        end
        
    end
    
    if BAnalysis == 1
        figure(contfig);
        plot(emiCur,rmsB/1e-6);        %Plot active load power vs. frequency.
        xlabel('Emitter Peak Current (A)'); ylabel('RMS B-Field (uT)');
        title('Emitter Current Sensibility Analysis');
        saveas(figure(contfig),'rmsB_emiCur.fig')
        close(gcf);
        contfig = contfig + 1;
    end
    
    figure(contfig);
    plot(emiCur,real(SL));        %Plot active load power vs. Emitter peak current.
    xlabel('Emitter Peak Current (A)'); ylabel('Active Power (W)');
    title('Emitter Current Sensibility Analysis');
    saveas(figure(contfig),'SL_emiCur.fig')
    close(gcf);
    contfig = contfig + 1;
    
    figure(contfig);
    plot(emiCur,emiCalc/1e-6);        %Plot emitter self-inductance vs. emitter peak current.
    xlabel('Emitter Peak Current (A)'); ylabel('Emitter Inductance (uH)');
    title('Emitter Current Sensibility Analysis');
    saveas(figure(contfig),'emiCalc_emiCur.fig');
    close(gcf);
    contfig = contfig + 1;
      
    save('emiCur.mat','emiCur');save('SL.mat','SL'); save('emiCalc.mat','emiCalc'); movefile('emiCur.mat',subfolder);
    save('loadCur.mat','loadCur');save('recCur.mat','recCur');
    movefile('SL.mat',subfolder); movefile('emiCalc.mat',subfolder); movefile('SL_emiCur.fig',subfolder); movefile('emiCalc_emiCur.fig',subfolder);
    movefile('loadCur.mat',subfolder); movefile('recCur.mat',subfolder);
         
    if BAnalysis == 1
        save('rmsB.mat','rmsB'); movefile('rmsB_emiCur.fig',subfolder);  movefile('rmsB.mat',subfolder);
    end
    
    movefile(subfolder,folder);
    
    n = 0;
    SL = 0; emiCalc = 0; recCalc = 0; rmsB = 0;
    loadCur = 0; emiCur = 0; recCur = 0; sourceCurrentAngle = 0;
    
    closefemm;
    delete('temp_target.fem');
    delete('temp_target.ans');

else
    execEmiCurrent = 'Emitter current sensibility analysis NOT performed.';
end

%% Receiver Current Sensibility Analysis
clc;
recCurrentAnalysis = input('Enter "1" to perform receiver current sensibility analysis or "0" to ignore it:');

if recCurrentAnalysis == 1
    
    subfolder = 'Receiver Current Analysis';
    mkdir(subfolder);
    
    execRecCurrent = 'Receiver current sensibility analysis performed.';
    
    openfemm;                            %Open FEMM
    try
        vv=ver;
        opendocument([cd,targetFile2]);  %Open desired file.
    catch
        opendocument(targetFile);
    end
    mi_saveas('temp_target.fem');        %Save temporary copy.
    
    recCurVar= input('Enter receiver current variation (eg., 10, for +/- 10% of nominal value):');
    recCurrentMinimum = recCurrentRef*(1 - recCurVar/100);
    recCurrentMaximum = recCurrentRef*(1 + recCurVar/100);

    recCurrentRange = recCurrentMaximum - recCurrentMinimum;
    recCurrentStep = (recCurrentRange/2)/(ceil(detail/2));

    emiCalc = emiCalcRef; recCalc = recCalcRef; mutualCalc = mutualCalcRef; emiResistance = emiResistanceRef; recResistance = recResistanceRef;
    frequency = frequencyRef; w = 2*pi*frequencyRef; emiVoltage = emiVoltageRef;

    for recCurrent = recCurrentMinimum:recCurrentStep:recCurrentMaximum

        n = n + 1;
        recCur(n) = recCurrent;
        mode = 'recCurrentVar';
        resultsAgree = 0;
        recCalcPrev = recCalcRef; mutualCalcPrev = mutualCalcRef; recResistancePrev = recResistanceRef;
        mutualCalc = mutualCalcPrev; recResistance = recResistancePrev; recCalc = recCalcPrev;
        
        while resultsAgree == 0
                  
            modelTransf = tmodel (rint,emiType,recType,emiCalc,mutualCalc,emiCapacitance,w,emiResistance,recResistance,recCalc,recCapacitance,Load,emiVoltage,emiCurrent,mode,recCurrent); %Call T-Model function.   
            loadCurrent = modelTransf(3);
            
            mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));  %Update excitation.
            mi_setcurrent('receiver',abs(recCurrent));                       %Only emitter is excited.
            mi_setcurrent('emitter',0);
            mi_analyze(1);
            mi_loadsolution;
            v = mo_getcircuitproperties('receiver'); recCalc = abs(v(3)/v(1)); recResistance = real(v(2)/v(1));
            
            %NOTE: At this point, emitter parameters would be re-calculated. However, it is assumed that emitter inductance,
            % resistance and the coupling coefficient are fixed. This means that emitter is not subject to magnetic saturation.
            
            recDelta = abs(recCalc - recCalcPrev)/recCalcPrev;
            recResistanceDelta = abs(recResistance - recResistancePrev)/recResistancePrev;
            recCalcPrev = recCalc; recResistancePrev = recResistance;
            
            if ((abs(recCalc - recCalcRef)/recCalcRef) > 0.2)
                clc;
                disp('Possible magnetic saturation detected in emitter !');
            end

            if (recDelta < dev) && (recResistanceDelta < dev)
                resultsAgree = 1;
            end
        end
        
        loadCur(n) = loadCurrent; 
        emiCur(n) = emiCurrent;
        
        SL(n) = Load*loadCurrent*conj(loadCurrent);
        recCalc(n) = recCalc;
        
        if BAnalysis == 1
           Bfunction = bcalc(emiCurrent,recCurrent,coordX,coordY); 
           rmsB(n) = Bfunction(1);
        end
        
    end
   
    if BAnalysis == 1
        figure(contfig);
        plot(recCur,rmsB/1e-6);        %Plot active load power vs. frequency.
        xlabel('Receiver Peak Current (A)'); ylabel('RMS B-Field (uT)');
        title('Receiver Current Sensibility Analysis');
        saveas(figure(contfig),'rmsB_recCur.fig')
        close(gcf);
        contfig = contfig + 1;
    end
    
    figure(contfig);
    plot(recCur,real(SL));        %Plot active load power vs. receiver peak current.
    xlabel('Receiver Peak Current (A)'); ylabel('Active Power (W)');
    title('Receiver Current Sensibility Analysis');
    saveas(figure(contfig),'SL_recCur.fig')
    close(gcf);
    contfig = contfig + 1;
    
    figure(contfig);
    plot(recCur,recCalc/1e-6);        %Plot receiver self-inductance vs. receiver peak current.
    xlabel('Receiver Peak Current (A)'); ylabel('Receiver Inductance (uH)');
    title('Receiver Current Sensibility Analysis');
    saveas(figure(contfig),'recCalc_recCur.fig');
    close(gcf);
    contfig = contfig + 1;
      
    save('recCur.mat','recCur');save('SL.mat','SL'); 
    save('recCalc.mat','recCalc'); save('loadCur.mat','loadCur'); save('emiCur.mat','emiCur');
    movefile('recCur.mat',subfolder); movefile('SL.mat',subfolder); movefile('recCalc.mat',subfolder);
    movefile('SL_recCur.fig',subfolder); movefile('recCalc_recCur.fig',subfolder);
    movefile('loadCur.mat',subfolder);movefile('emiCur.mat',subfolder);
    
    
    if BAnalysis == 1
        save('rmsB.mat','rmsB'); movefile('rmsB_recCur.fig',subfolder);  movefile('rmsB.mat',subfolder);
    end
    
    movefile(subfolder,folder);
    
    n = 0;
    SL = 0; emiCalc = 0; recCalc = 0; rmsB = 0;
    loadCur = 0; emiCur = 0; recCur = 0; sourceCurrentAngle = 0;
    
    closefemm;
    delete('temp_target.fem');
    delete('temp_target.ans');

else
    execRecCurrent = 'Receiver current sensibility analysis NOT performed.';
end


%% Coupling Coefficient Sensibility Analysis
clc;
couplingAnalysis = input('Enter "1" to perform coupling coefficient sensibility analysis or "0" to ignore it:');

if couplingAnalysis == 1

    subfolder = 'Coupling Coefficient Analysis';
    mkdir(subfolder);
    
    execCoupling = 'Coupling coefficient sensibility analysis performed.';
    openfemm;                            %Open FEMM
    try
        vv=ver;
        opendocument([cd,targetFile2]);  %Open desired file.
    catch
        opendocument(targetFile);
    end
    mi_saveas('temp_target.fem');        %Save temporary copy.
    
    distanceIncrease = input('Define increase in nominal emitter-receiver distance (mm):');
    distanceDecrease = input('Define decrease in nominal emitter-receiver distance (mm):');

    distanceStep = (distanceIncrease - distanceDecrease)/detail;
    distanceDecrease = - distanceDecrease; %Negative, because it is below the nominal distance.
    
    emiCalc = emiCalcRef; recCalc = recCalcRef; mutualCalc = mutualCalcRef; emiResistance = emiResistanceRef; recResistance = recResistanceRef;
    emiCalcPrev = emiCalcRef; recCalcPrev = recCalcRef; mutualCalcPrev = mutualCalcRef; emiResistancePrev = emiResistanceRef; recResistancePrev = recResistanceRef;
    frequency = frequencyRef; w = 2*pi*frequencyRef; emiVoltage = emiVoltageRef;
    
    %First, positive emitter-receiver distance variation.
    for distance = distanceDecrease:distanceStep:distanceIncrease           %Assume that the loaded file representes receiver in nominal position.

        if n == 0
            mi_selectgroup(2);
            mi_movetranslate(0,distanceDecrease);
            mi_clearselected;
        end
        
        n = n + 1;
        dist(n) = distance;
        mode = 'normal';
        resultsAgree = 0;
        %emiCalcPrev = emiCalcRef; recCalcPrev = recCalcRef; mutualCalcPrev = mutualCalcRef; emiResistancePrev = emiResistanceRef; recResistancePrev = recResistanceRef;
        emiCalc = emiCalcPrev; recCalc = recCalcPrev; mutualCalc = mutualCalcPrev; emiResistance = emiResistancePrev; recResistance = recResistancePrev; emiVoltage = emiVoltageRef;

        if updatePosition == 1
            updatePosition = 0;
            mi_selectgroup(2);
            mi_movetranslate(0,distanceStep);
            mi_clearselected;
        end
        
        while resultsAgree == 0
        
            modelTransf = tmodel (rint,emiType,recType,emiCalc,mutualCalc,emiCapacitance,w,emiResistance,recResistance,recCalc,recCapacitance,Load,emiVoltage,emiCurrent,mode,recCurrent); %Call T-Model function.   
            emiCurrent = modelTransf(1);%Outputs of the T-Model function are emitter current, receiver current and load current;
            recCurrent = modelTransf(2); 
            loadCurrent = modelTransf(3);
    
            mi_probdef(frequency,'millimeters',simType,1e-8,depth,20,(0));  %Update excitation.
            mi_setcurrent('emitter',abs(emiCurrent));                       %Only emitter is excited.
            mi_setcurrent('receiver',0);
            mi_analyze(1);
            mi_loadsolution;
            r = mo_getcircuitproperties('emitter'); emiCalc = abs(r(3)/r(1)); emiResistance = real(r(2)/r(1));              %Here, emiCalc is obtained "as is", considering 
            h = mo_getcircuitproperties('receiver');emiOpenCalc = abs(h(2)); mutualCalc = abs(emiOpenCalc/(w*emiCurrent));  %interactions with of nearby elements or and magnetic saturation.
            
            mi_setcurrent('emitter',0);
            mi_setcurrent('receiver',abs(recCurrent));                       %Now, only receiver is excited.
            mi_analyze(1);
            mi_loadsolution;
            v = mo_getcircuitproperties('receiver'); recCalc = abs(v(3)/v(1)); recResistance = real(v(2)/v(1));  %Also, recCalc is obtained "as is". 

            emiDelta = abs(emiCalc - emiCalcPrev)/emiCalcPrev; 
            recDelta = abs(recCalc - recCalcPrev)/recCalcPrev; 
            mutualDelta = abs(mutualCalc - mutualCalcPrev)/mutualCalcPrev;
            emiResistanceDelta = abs(emiResistance - emiResistancePrev)/emiResistancePrev; recResistanceDelta = abs(recResistance - recResistancePrev)/recResistancePrev;
        
            emiCalcPrev = emiCalc; recCalcPrev = recCalc; mutualCalcPrev = mutualCalc; emiResistancePrev = emiResistance; recResistancePrev = recResistance;
        
            if (emiDelta < dev) &&  (recDelta < dev) && (mutualDelta < dev) && (emiResistanceDelta < dev) && (recResistanceDelta < dev)
                resultsAgree = 1;
                updatePosition = 1;
            end
        end
        
        loadCur(n) = loadCurrent; 
        emiCur(n) = emiCurrent;
        recCur(n) = recCurrent;
        
        inputsCalc = inputPar(emiType,emiVoltage,emiCurrent,rint,emiCapacitance,w);
        sourceCurrentAngle(n) = inputsCalc(1);
        
        SL(n) = Load*loadCurrent*conj(loadCurrent);
        emiCalcInd(n) = emiCalc;
        recCalcInd(n) = recCalc;
        coupling(n) = mutualCalc/(sqrt(emiCalc*recCalc));
        
         if BAnalysis == 1
           Bfunction = bcalc(emiCurrent,recCurrent,coordX,coordY); 
           rmsB(n) = Bfunction(1);
        end
        
    end
    
    if BAnalysis == 1
        figure(contfig);
        plot(dist,rmsB/1e-6);        %Plot active load power vs. frequency.
        xlabel('Emitter-Receiver Distance (mm)'); ylabel('RMS B-Field (uT)');
        title('Coupling Sensibility Analysis');
        saveas(figure(contfig),'rmsB_coupling.fig')
        close(gcf);
        contfig = contfig + 1;
    end
    
    figure(contfig);
    plot(dist,sourceCurrentAngle);        %Plot active load power vs. frequency.
    xlabel('Emitter-Receiver Distance (mm)'); ylabel('Source Current Angle (deg)');
    title('Coupling Sensibility Analysis');
    saveas(figure(contfig),'sourceCurrent_coupling.fig')
    close(gcf);
    contfig = contfig + 1;
    
    figure(contfig);
    plot(dist,real(SL));                  %Plot active load power vs. emitter-receiver distance.
    xlabel('Emitter-Receiver Distance (mm)'); ylabel('Active Power (W)');
    title('Coupling Sensibility Analysis');
    saveas(figure(contfig),'SL_dist.fig')
    close(gcf);
    contfig = contfig + 1;
    
    figure(contfig);
    plot(dist,emiCalcInd/1e-6);        %Plot emitter self-inductance vs. emitter-receiver distance.
    xlabel('Emitter-Receiver Distance (mm)'); ylabel('Emitter Inductance (uH)');
    title('Coupling Coefficient Sensibility Analysis');
    saveas(figure(contfig),'emiCalc_dist.fig');
    close(gcf);
    contfig = contfig + 1;
    
    figure(contfig);
    plot(dist,recCalcInd/1e-6);        %Plot emitter self-inductance vs. emitter-receiver distance.
    xlabel('Emitter-Receiver Distance (mm)'); ylabel('Receiver Inductance (uH)');
    title('Coupling Coefficient Sensibility Analysis');
    saveas(figure(contfig),'recCalc_dist.fig');
    close(gcf);
    contfig = contfig + 1;
      
    figure(contfig);
    plot(dist,coupling);        %Plot emitter self-inductance vs. emitter-receiver distance.
    xlabel('Emitter-Receiver Distance (mm)'); ylabel('Coupling Coefficient');
    title('Coupling Coefficient Sensibility Analysis');
    saveas(figure(contfig),'coupling_dist.fig');
    close(gcf);
    contfig = contfig + 1;
    
    save('dist.mat','dist');save('SL.mat','SL'); save('coupling.mat','coupling');
    save('emiCalcInd.mat','emiCalcInd');save('recCalcInd.mat','recCalcInd');
    save('sourceCurrentAngle.mat','sourceCurrentAngle'); movefile('sourceCurrentAngle.mat',subfolder);
    movefile('dist.mat',subfolder); movefile('SL.mat',subfolder); movefile('emiCalcInd.mat',subfolder);
    movefile('recCalcInd.mat',subfolder); movefile('SL_dist.fig',subfolder); movefile('emiCalc_dist.fig',subfolder);
    movefile('recCalc_dist.fig',subfolder); movefile('coupling_dist.fig',subfolder); movefile('coupling.mat',subfolder);
    movefile('sourceCurrent_coupling.fig',subfolder);
    save('loadCur.mat','loadCur'); save('emiCur.mat','emiCur');save('recCur.mat','recCur');
    movefile('loadCur.mat',subfolder);movefile('emiCur.mat',subfolder);movefile('recCur.mat',subfolder);
    
    if BAnalysis == 1
        save('rmsB.mat','rmsB'); movefile('rmsB_coupling.fig',subfolder);  movefile('rmsB.mat',subfolder);
    end
    
    movefile(subfolder,folder);
    
    n = 0;
    SL = 0; emiCalc = 0; recCalc = 0; rmsB = 0;
    loadCur = 0; emiCur = 0; recCur = 0;
    
    closefemm;
    delete('temp_target.fem');
    delete('temp_target.ans');
    
else
execCoupling = 'Coupling coefficient sensibility analysis NOT performed.';
end

sensTime = toc(clockSens);  

datafile = fopen(fileName,'at'); 
fprintf(datafile,'\n ');
fprintf(datafile,'\n ****** Analysis Performed ******');
fprintf(datafile,'\n %s',execFrequency);
fprintf(datafile,'\n %s',execEmiCurrent);
fprintf(datafile,'\n %s',execRecCurrent);
fprintf(datafile,'\n %s',execCoupling);
fprintf(datafile,'\n ');
fprintf(datafile,'\n Elapsed time for Sensibility Analysis: %d minutes',sensTime/60);
fclose(datafile);
movefile(fileName,folder);
copyfile(targetFile,folder);
clc;
disp('*************************************');
disp('Sensitivity Analysis concluded.');
disp('*************************************');
                                               