clear
clc

tic

% warning off

%% Header

disp('=====================================================================================');
disp('+---------------+                                      ');
disp('| RYERSON       |  OPERA V1.0  []                         []  ');
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
% strSTL = 'CAD Geom/simple_liftingsurface.stl';

% strSTL = 'Cad Geom/wing_simple_short.stl';

% strSTL = 'Cad Geom/wing_simple.stl';

% strSTL = 'Cad Geom/quad.stl';
% strSTL = 'CAD Geom/quad-mix.stl';
strSTL = 'CAD Geom/quad-mix2.stl';
% strSTL = 'Cad Geom/quad-align.stl';
% strSTL = 'Cad Geom/quad-align-wing.stl';
% strSTL = 'Cad Geom/quad-align-wing-stretch.stl';

% strSTL = 'CAD Geom/2quad.stl';
% strSTL = 'CAD Geom/pyramid.stl';
%
% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

strA2TYPE = 'WING';
valMAXTIME = 0;
valDELTIME = 0.3;
flagRELAX = 0;
vecLE = []';
vecSYM = []';

vecTE = [];
vecTE = [26 19 12 3];
% vecTE = [52 45 38 31 24 17 3 5 61 68 75 82 89 96]';

seqALPHA = 30;
seqBETA = 0;

%% Triangulating Geometry

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
    matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnIMPORTGEOM(strSTL, strATYPE);

% trimesh(TR)
%% D-Matrix Creation

matD = fcnDWING7(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT);
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
        valWSIZE = [];
        matWCOEFF = [];
        matWVLST = [];
        
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
        % Building wing resultant
        vecR = fcnRWING(strATYPE, valDLEN, 0, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
        for valTIMESTEP = 1:valMAXTIME
            %% Timestep to solution
            %   Move wing
            %   Generate new wake elements
            %   Create and solve WD-Matrix for new elements
            %   Solve wing D-Matrix with wake-induced velocities
            %   Solve entire WD-Matrix
            %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
            %   Calculate surface normal forces
            %   Calculate DVE normal forces
            %   Calculate induced drag
            
            valWSIZE = length(vecTE);
            
            % Moving the wing
            [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matELST, vecTE);
            
            % Generating new wake elements
            if any(vecTE)
                [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
                    matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER, matWCOEFF] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, matWCOEFF, vecTE, matEATT);
                
                % Rebuild wing resultant
                vecR = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX);
                
                matCOEFF = fcnSOLVED(matD, vecR, valNELE);
                
                if flagRELAX == 1
                    [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER] ...
                        = fcnRELAX(valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX);
                end
            end
        end
    end
end

%% Plot
%
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
if any(vecTE)
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
end

%% End

% % matCOEFF = [0 1 0 1 1; 0 1 0 1 1];
% 
% i = 1;
% 
% point = [0 0 0]
% 
% l_point = fcnTOLOC(i, point, matDVE, matDVECT, matVLST)
% 
% 
% 
% eta = l_point(1);
% xsi = l_point(2);
% 
% circ = matCOEFF(i,1).*(eta.^2) + matCOEFF(i,2).*eta + matCOEFF(i,3).*(xsi.^2)  + matCOEFF(i,4).*xsi + matCOEFF(i,5)
% % gcirc = fcnTOGLOB(i, [eta xsi circ], matDVE, matDVECT, matVLST)
%     
% i = 2;
% 
% point = [0 0 0]
% 
% l_point = fcnTOLOC(i, point, matDVE, matDVECT, matVLST)
% 
% eta = l_point(1);
% xsi = l_point(2);
% 
% circ = matCOEFF(i,1).*(eta.^2) + matCOEFF(i,2).*eta + matCOEFF(i,3).*(xsi.^2) + matCOEFF(i,4).*xsi + matCOEFF(i,5)
% % gcirc = fcnTOGLOB(i, [eta xsi circ], matDVE, matDVECT, matVLST)

