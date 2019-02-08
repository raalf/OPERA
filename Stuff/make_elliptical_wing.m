% clc
clear

% strSPANSPACING = 'HALFCOSINE';
strSPANSPACING = 'COSINE';
strCHORDSPACING = 'NORMAL';
valALPHA = 4;
AR = 7;

valMAXTIME = 20;
valDELTIME = 1;

vecM = 6;
vecN = 11;

xtcr = 0; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr)
run('OPERA_MAIN');