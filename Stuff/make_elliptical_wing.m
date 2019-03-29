% clc
clear

strSPANSPACING = 'HALFCOSINE';
% strSPANSPACING = 'COSINE';
strCHORDSPACING = 'HALFCOSINE';
valALPHA = 20;
AR = 7;

valMAXTIME = 60;
valDELTIME = 0.1;

vecM = 4;
vecN = 6;

xtcr = 1; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr, 'NO')
run('OPERA_MAIN');