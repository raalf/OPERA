clear
% clc

%% Preamble
%
% strFILE = 'inputs/simple_wing.dat'
% strFILE = 'inputs/nonplanar.dat'
strFILE = 'inputs/2dve.dat';
% strFILE = 'inputs/4dve.dat';
% strFILE = 'inputs/4dve_nosym.dat'
% strFILE = 'inputs/Stock_Test1.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);
strATYPE = 'THIN';

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

matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM);

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

granularity = .05;
x = -0.6:granularity:1.2;
% y = -0.2:granularity:1.2;
y = zeros(size(x)) + 0.33;
z = -0.2:granularity:0.2;
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [-0.5 0.5 0.15];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);

% q_ind = s_ind + repmat(vecUINF,size(s_ind,1),1);
q_ind = s_ind;
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off

hold on
q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF;
quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 0.25, 'g')
hold off

view([0 90])

%%
v_row = 1:size(matVLST,1);

ii = 1;
iii = 1;
for jj = 1:length(v_row)
    
    [~, b] = ismembertol(v_row(jj), matDVE,'OutputAllIndices',true);

    [elem(:,1) elem(:,2)] = ind2sub(size(matDVE), cell2mat(b));

    for i = 1:size(elem,1)
        point = matPLEX(elem(i,2),:,elem(i,1));
        vort = [2.*matCOEFF(elem(i,1),3).*point(:,1) + matCOEFF(elem(i,1),4) + matCOEFF(elem(i,1),5).*point(:,2), 2.*matCOEFF(elem(i,1),1).*point(:,2) + matCOEFF(elem(i,1),2) + matCOEFF(elem(i,1),5).*point(:,1), point(:,2).*0];
        vort_glob(ii,:) = fcnSTARGLOB(vort, matROTANG(elem(i,1),1), matROTANG(elem(i,1),2), matROTANG(elem(i,1),3));
        point_glob(ii,:) = matVLST(v_row(jj),:);
        v_num(ii,:) = v_row(jj);
        ii = ii + 1;
    end
    
    vort_avg(iii,:) = mean(vort_glob(ii-size(elem,1):ii-1,:),1);
    point_avg(iii,:) = point_glob(ii-1,:);
    
    iii = iii + 1;
    elem = [];
    point = [];
    vort = [];
end

hFig28 = figure(28);
clf(28);

subplot(1,2,1)
scatter(v_num, vort_glob(:,1));
ylabel('Gamma_x','FontSize',15);
xlabel('Vertex Number','FontSize',15);

grid minor
box on
axis tight

subplot(1,2,2)
scatter(v_num, vort_glob(:,2));
ylabel('Gamma_y','FontSize',15);
xlabel('Vertex Number','FontSize',15);

grid minor
box on
axis tight








