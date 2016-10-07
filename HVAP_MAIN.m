clear
clc
% clf
tic

% Analysis Type and Geometry File

strATYPE = 'LS'; % Lifting Surface
% STL = 'CAD Geom/simple_liftingsurface.stl';
% strSTL = 'CAD Geom/quad.stl';
% strSTL = 'CAD Geom/2quad.stl';
strSTL = 'CAD Geom/pyramid.stl';

% STL = 'Cad Geom/lifting_split.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

strA2TYPE = 'WING';
valMAXTIME = 1;
valDELTIME = 0.3;
vecTE = []';
vecSYM = []';

seqALPHA = deg2rad(5);
seqBETA = 0;

%% Triangulating Geometry

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
    matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER] = fcnIMPORTGEOM(strSTL, strATYPE);

%% D-Matrix Creation

matD = fcnDWING3(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecSYM, matVATT);
valDLEN = length(matD);

%% Alpha Loop
for ai = 1:length(seqALPHA)
    valALPHA = deg2rad(seqALPHA(ai));
    for bi = 1:length(seqBETA)
        valBETA = deg2rad(seqBETA(bi));
        
        vecUINF = fcnUINFWING(valALPHA, valBETA);  
        
        % Building wing resultant
        vecR = fcnRWING(strATYPE, valDLEN, 0, matEATT, matCENTER, matDVECT, vecUINF, vecTE);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
        matWAKEGEOM = [];
%         for valTIMESTEP = 1:valMAXTIME
%             %% Timestep to solution
%             %   Move wing
%             %   Generate new wake elements
%             %   Create W-Matrix and W-Resultant
%             %   Solve W-Matrix
%             %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
%             %   Calculate surface normal forces
%             %   Calculate DVE normal forces
%             %   Calculate induced drag
%             %   Calculate cn, cl, cy, cdi
%             %   Calculate viscous effects
%             
%             % Moving the wing
%             [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matELST, vecTE);
%             
%             % Generating new wake elements
% %             [matWAKEGEOM, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM);
%             
%             
%             
%         end
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
% [hFig1] = fcnPLOTWAKE(0, WDVE, WNELE, WVLST, WELST, WDVECT, WCENTER);
%% End


toc

% whos
