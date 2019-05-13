function [matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, valCOLL, valRPM, valJ] = fcnOPREAD(strFILE)

fp = fopen(strFILE);

%% Analysis Parameters
% Wing or rotor?
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
strATYPE = fscanf(fp,'%s',1);

% Thin or thick surface?
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
strA2TYPE = fscanf(fp,'%s',1);

% Wing or rotor?
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
strA3TYPE = fscanf(fp,'%s',1);
strATYPE = {strATYPE, strA2TYPE, strA3TYPE};

% Chordwise element spacing
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
strSPACING = fscanf(fp,'%s',1);

% Maximum timesteps
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valMAXTIME = fscanf(fp,'%lf');

% Timestep size
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valDELTIME = fscanf(fp,'%lf');

% Fixed or relax?
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
strRELAX = fscanf(fp,'%s',1);

flagRELAX = 0;
if strcmp(strRELAX, 'RELAXED')
    flagRELAX = 1;
end

% Reading sequence of alphas
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
seqALPHA = fscanf(fp,'%lf');

% Reading sequence of sideslip angles
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
seqBETA = fscanf(fp,'%lf');

% Reading density
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valDENSITY = fscanf(fp,'%lf');

if strcmpi(strATYPE{1}, 'ROTOR')
    valAREA = nan;
    valSPAN = nan;
    
    % Diameter
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valDIAM = fscanf(fp,'%lf');
    
    % RPM
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valRPM = fscanf(fp,'%lf'); 
    
    % Collective
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valCOLL = fscanf(fp,'%lf'); 
    
    % Advance ratio
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valJ = fscanf(fp,'%lf'); 
else
    valDIAM = nan;
    valRPM = nan;
    valCOLL = nan;
    valJ = nan;
    
    % Reading reference area
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valAREA = fscanf(fp,'%lf');

    % Reading reference span
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    valSPAN = fscanf(fp,'%lf');    
end

%% Reading panel/wing/lifting line information
% Reading No. of panels
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valPANELS = fscanf(fp,'%lf');

%% Reading panel information and geometry

vecN = zeros(valPANELS,1);
vecM = zeros(valPANELS,1);
strAIRFOIL = [];
vecSYM = zeros(valPANELS,1);

for i = 1:valPANELS
    % Reading spanwise 'n'
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecN(i) = fscanf(fp,'%lf',1);
    
    % Reading chordwise 'm'
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecM(i) = fscanf(fp,'%lf',1);
    
    % Reading symmetry information
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    strSYM = fscanf(fp,'%s',1);
    
    if strcmp(strSYM, 'YES')
        vecSYM(i) = 1;
    else
        vecSYM(i) = 0;
    end
    
    % Reading trailing edge information
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    strTE = fscanf(fp,'%s',1);
    
    if strcmp(strTE, 'YES')
        vecPANELTE(i) = 1;
    else
        vecPANELTE(i) = 0;
    end
    
    % Reading leading edge information
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    strLE = fscanf(fp,'%s',1);
    
    if strcmp(strLE, 'YES')
        vecPANELLE(i) = 1;
    else
        vecPANELLE(i) = 0;
    end
    
    % Reading spanwise DVE spacing for this panel
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    strPSPACE{i} = fscanf(fp,'%s',1);    
    
    % Skipping geometry column header
    fgets(fp);
    fgets(fp);
    
    % Reading geometry
    % Explanation below:
    %{
        info_geometry(x,y,z)
            x is for the left or right point
                1 left
                2 right
            y is for the values
                1 x
                2 y
                3 z
                4 chord
                5 epsilon
                6 boundary condition
            z is panel number
    %}

    matGEOM(1,:,i) = fscanf(fp,'%lf',5);
    try
        strAIRFOIL{i,1} = strtrim(fgetl(fp));
    end
        matGEOM(2,:,i) = fscanf(fp,'%lf',5);
    try
        strAIRFOIL{i,2} = strtrim(fgetl(fp));
    end
    
end

fclose(fp);
[matPOINTS, matTEPOINTS, matLEPOINTS, vecULS, vecDVESYM] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM, vecPANELTE, vecPANELLE, strATYPE, strAIRFOIL, strSPACING, strPSPACE);


end

