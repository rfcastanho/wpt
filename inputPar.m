% Input Parameters calculation

function inputsCalc = inputPar(emiType,emiVoltage,emiCurrent,rint,emiCapacitance,w)

if (emiType == 1) %SS and SP Topologies
    
    sourceCurrent = emiCurrent; %Emitter current is the source output current;
    
elseif (emiType == 0) %PP and PS Topologies
    
    emiCapCurrent = (emiVoltage - rint*emiCurrent)/(1/(1i*w*emiCapacitance)+rint);
    sourceCurrent = emiCapCurrent + emiCurrent;
    
end

[ang,mod] = cart2pol(real(sourceCurrent),imag(sourceCurrent));
angle = radtodeg(ang);

%Here, could calculate input power with emiVoltage and sourceCurrent

inputsCalc(1) = angle;

end