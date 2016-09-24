clear
clc
clf
tic

% Analysis Type and Geometry File

ATYPE = 'LS'; % Lifting Surface
STL = 'CAD Geom/simple_liftingsurface.stl';
% STL = 'CAD Geom/quad.stl';

% STL = 'Cad Geom/lifting_split.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

% A2TYPE = 'ROT';
A2TYPE = 'WING';
valMAXTIME = 5;
valDELTIME = 0.3;
vecTE = [3];

seqALPHA = deg2rad(5);
seqBETA = 0;

%% Triangulating Geometry

[TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, ...
    PLEX, DVECT, ALIGN, VATT, VNORM, CENTER] = fcnIMPORTGEOM(STL, ATYPE);

%% D-Matrix Creation

[D, R] = fcnDWING(EATT, PLEX, NELE, ELOC, ALIGN, VLST, VNORM, CENTER, DVE, DVECT);


%% Alpha Loop
for ai = 1:length(seqALPHA)
    valALPHA = seqALPHA(ai);
    for bi = 1:length(seqBETA)
        valBETA = seqBETA(bi);
        
        matWAKEGEOM = [];
        for valTIMESTEP = 1:valMAXTIME
            %% Timestep to solution
            %   Move wing
            %   Generate new wake elements
            %   Create W-Matrix and W-Resultant
            %   Solve W-Matrix
            %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
            %   Calculate surface normal forces
            %   Calculate DVE normal forces
            %   Calculate induced drag
            %   Calculate cn, cl, cy, cdi
            %   Calculate viscous effects
            
            % Moving the wing or rotor
            if strcmp(A2TYPE,'WING') == 1
                [VLST, CENTER, VUINF, CUINF, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, VLST, CENTER, ELST, vecTE);
            elseif strcmp(A2TYPE,'ROT') == 1
                [VLST, CENTER, VUINF, CUINF, matNEWWAKE] = fcnMOVEROTOR(ROTORIG, DELROT, ROTRPM, VLST, CENTER, ELST, vecTE);
            end
            
            % Generating new wake elements
            [matWAKEGEOM, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM);
            
            
            
        end
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER);
[hFig1] = fcnPLOTWAKE(0, WDVE, WNELE, WVLST, WELST, WDVECT, WCENTER);
%% End


toc

% whos
