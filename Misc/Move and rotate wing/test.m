clc
clear

load('LS.mat')

clf;

% type = 'ROT';
type = 'WING';

vecTE = [3];

if strcmp(type, 'ROT') == 1
    % moving the wing so we can rotate around z-axis
    translation = [-0.5 2 0];
    VLST = VLST + repmat(translation, length(VLST(:,1)),1);
    CENTER = CENTER + repmat(translation, length(CENTER(:,1)),1);
    
    [hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER);
    
    ROTORIG = [0 0 0];
    DELROT = pi/10;
    ROTRPM = 1;
    
    [VLST, CENTER, VUINF, CUINF, matNEWWAKE] = fcnMOVEROTOR(ROTORIG, DELROT, ROTRPM, VLST, CENTER, ELST, vecTE);
    
    hold on
    [hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER);
    hold off
    
else
    UINF = 1;
    ALPHA = deg2rad(10);
    BETA = deg2rad(10);
    DELTIME = 2;
    
    [hFig1] = fcnPLOTBODY(0, DVE, NELE, VLST, ELST, DVECT, CENTER);
    
    [VLST, CENTER, VUINF, CUINF, matNEWWAKE] = fcnMOVEWING(ALPHA, BETA, DELTIME, VLST, CENTER, ELST, vecTE);
    
    hold on
    [hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER);
    hold off
    
end

[WTR, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnTRIANG(matNEWWAKE);

[hFig1] = fcnPLOTWAKE(1, WDVE, WNELE, WVLST, WELST, WDVECT, WCENTER);
