% clc
clear

% strSPANSPACING = 'NORMAL';
strSPANSPACING = 'HALFCOSINE';
strCHORDSPACING = 'NORMAL';
valALPHA = 4;
AR = 7;

valMAXTIME = 60;
valDELTIME = 0.2;

vecM = 2;
vecN = 5;

xtcr = 1; % 1 is straight TE, 0 is straight LE

elliptical_wing_o_matic(vecN, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr, 'NO')
run('OPERA_MAIN');