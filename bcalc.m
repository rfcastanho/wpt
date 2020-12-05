%B-Field Calculation

function Bfunction = bcalc(emiCurrent,recCurrent,coordX,coordY)

mi_setcurrent('emitter',emiCurrent);  %Emitter and receiver excited.
mi_setcurrent('receiver',recCurrent);
mi_analyze(1);
mi_loadsolution;

bb = mo_getb(coordX,coordY);
Bx = abs(bb(1)); By = abs(bb(2)); 
modB = sqrt(Bx^2 + By^2); 
rmsB = modB/(sqrt(2)); %RMS value of B in desired reference point.
                                                                                   %RMS is chosen because ICNIRP 1998 uses it for reference.
 
Bfunction(1) = rmsB;
                                                                                   
end