clc
clear

load('LS.mat')

clf;

type = 'ROT';
% type = 'WING';

vecTRAILINGDVE = 1;
vecTRAILINGEDGE = [3];

if strcmp(type, 'ROT') == 1
    % moving the wing so we can rotate around z-axis
    translation = [-0.5 2 0];
    VLST = VLST + translation;
    CENTER = CENTER + translation;
    
    [hFig1] = fcnPLOTBODY(0, DVE, NELE, VLST, ELST, DVECT, CENTER);
    
    ROTORIG = [0 0 0];
    DELROT = pi/10;
    ROTRPM = 1;
    
    % UINF for vertex list
    radius = VLST - ROTORIG;
    radius = sum(abs(sqrt(radius(:,1:2).^2)),2);
    
    VUINF = 2.*pi.*ROTRPM.*radius;
    
    % UINF for control point locations
    radius = CENTER - ROTORIG;
    radius = sum(abs(sqrt(radius(:,1:2).^2)),2);
    
    CUINF = 2.*pi.*ROTRPM.*radius;
    
    % Rotatio
    % translate space so that the rotation axis passes through the origin
    % perform the desired rotation by theta about the z axis
    % apply the inverse of step (1)
    
    tempVLST = VLST - ROTORIG;
    tempCENTER = CENTER - ROTORIG;
    
    ROT = [cos(DELROT) -sin(DELROT) 0; sin(DELROT) cos(DELROT) 0; 0 0 1];
    
    VLST2 = (ROT*tempVLST')' + ROTORIG;
    CENTER2 = (ROT*tempCENTER')' + ROTORIG;
    
    % Spit out wake
    
    % Old trailing edge of wing
    old_te = VLST(ELST(vecTRAILINGEDGE,:),:);
    
    % New trailing edge of wing
    new_te = VLST2(ELST(vecTRAILINGEDGE,:),:);
    
    temp(:,:,1) = [new_te(1:end/2,:); old_te(1:end/2,:)];
    temp(:,:,2) = [new_te((end/2)+1:end,:); old_te((end/2)+1:end,:)];
    temp(:,:,3) = [old_te(1:end/2,:); new_te((end/2)+1:end,:)];
    
    % Creating wake HDVEs
    [WTR, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnTRIANG(temp);
   
    VLST = VLST2;
    CENTER = CENTER2;
    
    hold on
    [hFig1] = fcnPLOTBODY(0, DVE, NELE, VLST, ELST, DVECT, CENTER);
    hold off
    
else
    UINF = 1;
    ALPHA = deg2rad(10);
    BETA = deg2rad(10);
    
    [hFig1] = fcnPLOTBODY(0, DVE, NELE, VLST, ELST, DVECT, CENTER);
    
    UINF = [UINF*cos(ALPHA)*cos(BETA) UINF*sin(BETA) UINF*sin(ALPHA)*cos(BETA)];
    
    DELTIME = 2;
    
    translation = DELTIME.*UINF;
    
    VLST = VLST + translation;
    CENTER = CENTER + translation;
    
    hold on
    
    [hFig1] = fcnPLOTBODY(0, DVE, NELE, VLST, ELST, DVECT, CENTER);
    hold off
    
end