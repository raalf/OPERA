% clc
clear

strSPANSPACING = 'HALFCOSINE';
% strSPANSPACING = 'COSINE';
strCHORDSPACING = 'HALFCOSINE';
valALPHA = 20;
AR = 7;

valMAXTIME = 60;
valDELTIME = 0.1;

vecM = 1;
vecN = 7;

xtcr = 1; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr, 'YES')
run('OPERA_MAIN');