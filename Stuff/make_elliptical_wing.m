% clc
clear

strSPANSPACING = 'HALFCOSINE';
% strSPANSPACING = 'COSINE';
strCHORDSPACING = 'HALFCOSINE';
valALPHA = 4;
AR = 7;

valMAXTIME = 60;
valDELTIME = 0.1;

vecM = 6;
vecN = 12;

xtcr = 0; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr, 'YES')
run('OPERA_MAIN');