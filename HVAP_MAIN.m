clear
clc
% clf
tic

% Analysis Type and Geometry File

ATYPE = 'LS'; % Lifting Surface
% STL = 'CAD Geom/simple_liftingsurface.stl';
% STL = 'CAD Geom/quad.stl';
% STL = 'CAD Geom/2quad.stl';
STL = 'CAD Geom/pyramid.stl';

% STL = 'Cad Geom/lifting_split.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

% A2TYPE = 'ROT';
A2TYPE = 'WING';
valMAXTIME = 1;
valDELTIME = 0.3;
vecTE = []';
vecSYM = []';

seqALPHA = deg2rad(5);
seqBETA = 0;

%% Triangulating Geometry

[TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, ...
    PLEX, DVECT, ALIGN, VATT, VNORM, CENTER] = fcnIMPORTGEOM(STL, ATYPE);

%% D-Matrix Creation

matD = fcnDWING2(ATYPE, EATT, PLEX, NELE, ELOC, ELST, ALIGN, VLST, CENTER, DVE, DVECT, vecTE, vecSYM);
valDLEN = length(matD);

%% Alpha Loop
for ai = 1:length(seqALPHA)
    valALPHA = seqALPHA(ai);
    for bi = 1:length(seqBETA)
        valBETA = seqBETA(bi);
        
        if strcmp(A2TYPE,'WING') == 1
            [matVUINF, matCUINF] = fcnUINFWING(valALPHA, valBETA, VLST, CENTER);
        elseif strcmp(A2TYPE,'ROT') == 1
            [matVUINF, matCUINF] = fcnUINFROTOR(vecROTORIG, valROTRPM, matVLST, matCENTER);
        end
        
        % Building wing resultant
        vecR = fcnRWING(ATYPE, valDLEN, 0, VLST, ELST, EATT, CENTER, DVECT, matVUINF, matCUINF, vecTE);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, NELE);
        
        
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
%                 [VLST, CENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, VLST, CENTER, ELST, vecTE);
            elseif strcmp(A2TYPE,'ROT') == 1
%                 [VLST, CENTER, matNEWWAKE] = fcnMOVEROTOR(ROTORIG, DELROT, ROTRPM, VLST, CENTER, ELST, vecTE);
            end
            
            % Generating new wake elements
%             [matWAKEGEOM, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM);
            
            
            
        end
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER, PLEX, matCOEFF);
% [hFig1] = fcnPLOTWAKE(0, WDVE, WNELE, WVLST, WELST, WDVECT, WCENTER);
%% End


toc

% whos
