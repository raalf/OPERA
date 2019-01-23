% clc
clear

strSPANSPACING = 'HALFCOSINE';
strCHORDSPACING = 'HALFCOSINE';
valALPHA = 4;
AR = 7;

valMAXTIME = 20;
valDELTIME = 1;

vecM = 5;
vecN = 8;

xtcr = 0.25; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr)
OPERA_MAIN