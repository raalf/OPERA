function [flgRELAX, flgSTEADY, valMAXTIME, valDELTIME, valDENSITY, valUINF, valALPHA, valBETA, ...
    valROLL, valFPA, valTRACK, valAREA, valSPAN, matPOINTS, matTEPOINTS, matLEPOINTS, ...
    vecDVESURFACE, vecDVEFLIP, valROTORS, valWINGS, vecROTORRPM, vecROTORDIAM, ...
    vecROTORBLADES, matROTORHUB, matROTORAXIS, vecDVEWING, vecDVEROTOR] = fcnXMLREAD(filename)

inp = fcnXML2STRUCT(filename);
VAP = inp.VAP;

%% Settings
if strcmpi(VAP.settings.relax.Text, 'true') flgRELAX = true; else flgRELAX = false; end
if strcmpi(VAP.settings.steady.Text, 'true') flgSTEADY = true; else flgSTEADY = false; end
valMAXTIME = floor(str2double(VAP.settings.maxtime.Text));
valDELTIME = str2double(VAP.settings.delta_time.Text);

%% Conditions
valDENSITY = str2double(VAP.conditions.density.Text);
valUINF = str2double(VAP.conditions.speed.Text);

valALPHA = str2double(VAP.conditions.alpha.Text);
valBETA = str2double(VAP.conditions.beta.Text);
valROLL = str2double(VAP.conditions.roll.Text);
valFPA = str2double(VAP.conditions.fpa.Text);
valTRACK = str2double(VAP.conditions.track.Text);

valAREA = str2double(VAP.conditions.ref_area.Text);
valSPAN = str2double(VAP.conditions.ref_span.Text);


%%
vecROTORRPM = [];
vecROTORDIAM = [];
vecROTORBLADES = [];
matROTORHUB = []; 
matROTORAXIS = [];

%% Loading in wings and rotors

try valWINGS = max(size(VAP.wing)); catch; valWINGS = 0; end
try valROTORS = max(size(VAP.rotor)); catch; valROTORS = 0; end

k = 1;
kk = 1;
kkk = 1;
o = 1;
p = 1;

% Loading Wings
for j = 1:valWINGS
    
    try win = VAP.wing{1,j}; catch; win = VAP.wing; end
    
    tmpWINGINCID(k,1) = str2double(win.incidence.Text);
    
    try tmpWINGORIG(k,:) = [str2double(win.vehicle_x.Text) str2double(win.vehicle_y.Text) str2double(win.vehicle_z.Text)];
    catch; tmpWINGORIG(k,:) = [0 0 0]; end
    
    vecPANELS(k,1) = max(size(win.panel));
    for m = 1:vecPANELS(k,1)
        
        try pan = win.panel{1,m}; catch; pan = win.panel; end
        
        vecNtemp(kk,1) = floor(str2double(pan.spanwise_elements.Text));
        vecMtemp(kk,1) = str2double(win.chordwise_elements.Text);
        
        vecSECTIONS(kk,1) = max(size(pan.section));
        for n = 1:vecSECTIONS(kk,1)
            sec = pan.section{1,n};
            
            matSECTIONS(kkk,:) = [str2double(sec.wing_x.Text) + tmpWINGORIG(k,1) str2double(sec.wing_y.Text) + tmpWINGORIG(k,2) str2double(sec.wing_z.Text) + tmpWINGORIG(k,3) str2double(sec.chord.Text) tmpWINGINCID(k)+str2double(sec.twist.Text)];
            try tmpAIRFOIL{kkk,:} = sec.airfoil.Text; catch tmpAIRFOIL{kkk,:} = 'NACA 0015'; end
            
            vecSECTIONPANEL(kkk,1) = kk;
            
            kkk = kkk + 1;
        end
        
        panel_wing(kk,1) = o;
        panel_rotors(kk,1) = 0;
        
        kk = kk + 1;
    end
    
    o = o + 1;
    k = k + 1;
end


%% Loading Rotors
for j = 1:valROTORS
    
    try rot = VAP.rotor{1,j}; catch; rot = VAP.rotor; end
    
    vecROTORRPM(p,1) = str2double(rot.rpm.Text);
    vecROTORDIAM(p,1) = str2double(rot.ref_diam.Text);
    
    matROTORHUB(p,:) = [str2double(rot.veh_x_hub.Text) str2double(rot.veh_y_hub.Text) str2double(rot.veh_z_hub.Text)];
    matROTORAXIS(p,:) = [str2double(rot.veh_x_axis.Text) str2double(rot.veh_y_axis.Text) str2double(rot.veh_z_axis.Text)];
    
    vecROTORBLADES(p,:) = floor(str2double(rot.blades.Text));
    
    vecROTORM(p,1) = floor(str2double(rot.chordwise_elements.Text));
    
    try vecCOLLECTIVE(p,1) = str2double(rot.collective.Text); catch; vecCOLLECTIVE(p,1) = 0; end;
    
    vecPANELS(k,1) = max(size(rot.panel));
    
    flip = 1;
    try if strcmpi(rot.rotation_direction.Text,'CW'); flip = -1; end; end
    vecROTORRPM(p,1) = vecROTORRPM(p,1)*flip;
    
    for m = 1:vecPANELS(k,1)
        
        try pan = rot.panel{1,m}; catch; pan = rot.panel; end
        
        vecNtemp(kk,1) = floor(str2double(pan.spanwise_elements.Text));
        vecMtemp(kk,1) = str2double(rot.chordwise_elements.Text);
        
        vecSECTIONS(kk,1) = max(size(pan.section));
        
        for n = 1:vecSECTIONS(kk,1)
            sec = pan.section{1,n};
            
            matSECTIONS(kkk,:) = [str2double(sec.rotor_x.Text) + matROTORHUB(p,1) flip*str2double(sec.rotor_y.Text) + matROTORHUB(p,2) str2double(sec.rotor_z.Text) + matROTORHUB(p,3) str2double(sec.chord.Text) str2double(sec.twist.Text) + vecCOLLECTIVE(p,1)];
            try tmpAIRFOIL{kkk,:} = sec.airfoil.Text; catch tmpAIRFOIL{kkk,:} = 'NACA 0015'; end
            vecSECTIONPANEL(kkk,1) = kk;
            
            kkk = kkk + 1;
        end
        
        panel_wing(kk,1) = 0;
        panel_rotors(kk,1) = p;
        
        kk = kk + 1;
    end
    
    p = p + 1;
    k = k + 1;
end
valPANELS = sum(vecPANELS);

%% Reorganizing data for output
k = 1;
for i = 1:valPANELS
    sections = matSECTIONS(vecSECTIONPANEL == i,:);
    len = size(sections,1);
    foils = {tmpAIRFOIL{vecSECTIONPANEL == i}}';
    
    for j = 1:len - 1
        matGEOM(1,:,k) = sections(j,:);
        matGEOM(2,:,k) = sections(j+1,:);
        cellAIRFOIL{k,1} = foils{j};
        cellAIRFOIL{k,2} = foils{j+1};
        vecPANELWING(k,1) = panel_wing(i);
        vecPANELROTOR(k,1) = panel_rotors(i);
        
        vecN(k,1) = vecNtemp(i);
        vecM(k,1) = vecMtemp(i);
        k = k + 1;
    end
end

valPANELS = size(matGEOM,3);


%% Generate DVEs
[matPOINTS, matTEPOINTS, matLEPOINTS, vecDVESURFACE, vecDVEFLIP, vecDVEWING, vecDVEROTOR] = fcnGENERATEDVES(valPANELS, matGEOM, vecN, vecM, cellAIRFOIL, vecPANELWING, vecPANELROTOR);

end








