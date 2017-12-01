clear
% clc

%% Preamble
%
% strFILE = 'inputs/simple_wing.dat'
% strFILE = 'inputs/nonplanar.dat'
% strFILE = 'inputs/2dve.dat';
strFILE = 'inputs/4dve.dat';
% strFILE = 'inputs/Stock_Test1.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG(strATYPE, matPOINTS);

valTIMESTEP = 0;
flagRELAX = 0;

valMAXTIME = 0;
valALPHA = 20

vecLE = [];
vecTE = [];
try
[vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST);
end
%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = [];

vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,1), matROTANG(vecTEDVE,2), matROTANG(vecTEDVE,3)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM, matVSCOMB);

valDLEN = length(matD);


vecUINF = fcnUINFWING(valALPHA, 0);
matUINF = repmat(vecUINF, valNELE, 1);

% Building wing resultant
vecR = fcnRWING(strATYPE, valDLEN, 0, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, 0, [], [], [], [], [], [], [], matVNORM, matVLST);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

%%

% % Solving for wing coefficients
% [matCOEFF] = [...
%     0 0 0 -1 -1 0;
%     0 0 0 -1 -1 0; ...
%     ];

% matCOEFF = [matCOEFF(:,4:6) matCOEFF(:,1:3)]

% matCOEFF = [0 1 0 0 0 0; 0 1 0 0 0 0];
% matCOEFF = repmat(matCOEFF(1,:),2,1);


%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
view([-30 17])

% granularity = .2;
% x = -0.6:granularity:1.2;
% y = -0.2:granularity:1.2;
% % y = ones(size(x)) - 0.5;
% z = -0.2:granularity:0.2;
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% 
% % fpg = [-0.5 0.5 0.15];
% 
% [s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);
% 
% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% q_ind = s_ind;
% hold on
% % quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off

view([0 90])







