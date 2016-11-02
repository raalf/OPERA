clear
clc

tic

%% Header

disp('=====================================================================================');
disp('+---------------+                                      ');
disp('| RYERSON       |  HVAP V1.0   []                         []  ');
disp('| APPLIED       |              ||   ___     ___     ___   ||  ');
disp('| AERODYNAMICS  |              ||  /   \   /| |\   /   \  ||  ');
disp('| LABORATORY OF |              || |  O  |__|| ||__|  O  | ||  ');
disp('| FLIGHT        |              ||  \___/--/^^^^^\--\___/  ||  ');
disp('+---------------+              ||________|       |________||  __');
disp('        .-----------------/  \-++--------|   .   |--------++-/  \-----------------. ');
disp('       /.---------________|  |___________\__(*)__/___________|  |________---------.\');
disp('                 |    |   ''$$''   |                       |   ''$$''   |    |       ');
disp('                (o)  (o)        (o)                     (o)        (o)  (o)         ');
disp('=====================================================================================');

%% Preamble
% Analysis Type and Geometry File

strATYPE = 'LS'; % Lifting Surface
strSTL = 'CAD Geom/simple_liftingsurface.stl';
% strSTL = 'CAD Geom/quad.stl';
% strSTL = 'CAD Geom/2quad.stl';
% strSTL = 'CAD Geom/pyramid.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

strA2TYPE = 'WING';
valMAXTIME = 6;
valDELTIME = 0.3;
vecTE = [3]';
vecSYM = []';

seqALPHA = 30;
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
        
        matWADJE = [];
        matWELST = [];
        matWDVE = [];
        valWNELE = [];
        matWEATT = [];
        matWEIDX = [];
        matWELOC = [];
        matWPLEX = [];
        matWDVECT = [];
        matWALIGN = [];
        matWVATT = [];
        matWVNORM = [];
        matWCENTER = [];
        matWAKEGEOM = [];
        
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
        % Building wing resultant
        vecR = fcnRWING(strATYPE, valDLEN, 0, matEATT, matCENTER, matDVECT, vecUINF, vecTE);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
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
            
            %             % Moving the wing
            [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matELST, vecTE);
            
            % Generating new wake elements
            if any(vecTE)
                [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
                    matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM);
            end
            
            
        end
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
if any(vecTE)
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
end

%% End


toc

% whos
