clear
clc

tic

% warning off

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
% strSTL = 'CAD Geom/simple_liftingsurface.stl';

% strSTL = 'Cad Geom/wing5.stl';

% strSTL = 'Cad Geom/quad.stl';
% strSTL = 'CAD Geom/quad-mix.stl';
% strSTL = 'Cad Geom/quad-align.stl';
strSTL = 'Cad Geom/quad-align-wing.stl';
% strSTL = 'Cad Geom/quad-align-wing-stretch.stl';

% strSTL = 'CAD Geom/2quad.stl';
% strSTL = 'CAD Geom/pyramid.stl';
% 
% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

strA2TYPE = 'WING';
valMAXTIME = 0;
valDELTIME = 0.3;
vecTE = [22 21 23]';
vecLE = [16 10 2]';
vecSYM = []';

seqALPHA = 50;
seqBETA = 0;

%% Triangulating Geometry

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
    matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER] = fcnIMPORTGEOM(strSTL, strATYPE);
matALIGN(11,:,:) = 1;

trimesh(TR)
%% D-Matrix Creation

matD = fcnDWING6(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT);
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
            end
        end
    end
end

%% Plot
% 
[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF);
if any(vecTE)
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
end

%% End

dve = 1;
vert = 1;

lambdas = [1 0 0; 0 1 0; 0 0 1];

lam = lambdas(find(matDVE(dve,:,1) == vert),:);

eta = (lam(1)*matPLEX(1,1,dve) + lam(2)*matPLEX(2,1,dve) + lam(3)*matPLEX(3,1,dve));
xi = (lam(1)*matPLEX(1,2,dve) + lam(2)*matPLEX(2,2,dve) + lam(3)*matPLEX(3,2,dve));

gamma = matCOEFF(dve,1)*(eta^2) + matCOEFF(dve,2)*eta + matCOEFF(dve,3)*(xi^2) + matCOEFF(dve,4)*xi + matCOEFF(dve,5)
vort_eta = matCOEFF(dve,1)*eta + matCOEFF(dve,2)
vort_xi = matCOEFF(dve,3)*xi + matCOEFF(dve,4)
toc

gamma_globe = fcnTOGLOB(dve, [eta xi gamma], matDVE, matDVECT, matVLST)

% whos
