function [matNEWWAKECOEFF] = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, vecSPANDIR, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT)

% This function finds the coefficients of the post-trailing edge elements, by setting the spanwise component equal to the opposite of the circulation along the trailing edge of the wing (spanwise)
% There is no chordwise change in wake coefficients now (steady only)

lamb_circ = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

lamb_vort = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1; ...
    ];

len = length(vecWLEDVE);

%% Circulation at TE of wing/LE of wake Prelim
% Evaluated at the mid-point of each edge which splits two HDVEs

%(x,y) of all three vertices of HDVEs at trailing edge of the wing and the corresponding leading edge of the post te wake HDVE
x1 = [reshape(matPLEX(1,1,vecTEDVE),len,1) reshape(matWPLEX(1,1,vecWLEDVE),len,1)];
x2 = [reshape(matPLEX(2,1,vecTEDVE),len,1) reshape(matWPLEX(2,1,vecWLEDVE),len,1)];
x3 = [reshape(matPLEX(3,1,vecTEDVE),len,1) reshape(matWPLEX(3,1,vecWLEDVE),len,1)];
y1 = [reshape(matPLEX(1,2,vecTEDVE),len,1) reshape(matWPLEX(1,2,vecWLEDVE),len,1)];
y2 = [reshape(matPLEX(2,2,vecTEDVE),len,1) reshape(matWPLEX(2,2,vecWLEDVE),len,1)];
y3 = [reshape(matPLEX(3,2,vecTEDVE),len,1) reshape(matWPLEX(3,2,vecWLEDVE),len,1)];


%% Midpoint of edge to set circulation of wake opposite to wing
lmb1 = reshape(lamb_circ([nonzeros(matELOC(vecTE,:)) ones(len,1)],1),len,2); % Ones(le,1) because wake LE is always local edge 1
lmb2 = reshape(lamb_circ([nonzeros(matELOC(vecTE,:)) ones(len,1)],2),len,2);
lmb3 = reshape(lamb_circ([nonzeros(matELOC(vecTE,:)) ones(len,1)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
resultant1 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma1 = [b1(:,2), b2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*3)-2],1,3) + repmat([0:2],len,1)]',[],1);
circ1 = zeros(len, valWNELE*3);
circ1(sub2ind(size(circ1),rows,col1)) = reshape(gamma1',[],1);

%% Vertices to set wake circulation opposite to wing (along TE)

[ta,tb,~] = find(matDVE(vecTEDVE,:,1) == repmat(matELST(vecTE,1),1,3));
[~,rb,] = sort(ta);
vnuma(:,1) = tb(rb);

[ta,tb,~] = find(matDVE(vecTEDVE,:,1) == repmat(matELST(vecTE,2),1,3));
[~,rb,] = sort(ta);
vnuma(:,2) = tb(rb);

[ta,tb,~] = find(matWDVE(vecWLEDVE,:,1) == repmat(matWELST(vecWLE,1),1,3));
[~,rb,] = sort(ta);
vnumb(:,1) = tb(rb);

[ta,tb,~] = find(matWDVE(vecWLEDVE,:,1) == repmat(matWELST(vecWLE,2),1,3));
[~,rb,] = sort(ta);
vnumb(:,2) = tb(rb);

%% First vertex
lmb1 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],1),len,2);
lmb2 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],2),len,2);
lmb3 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
resultant2 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma2 = [b1(:,2), b2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*3)-2],1,3) + repmat([0:2],len,1)]',[],1);
circ2 = zeros(len, valWNELE*3);
circ2(sub2ind(size(circ2),rows,col1)) = reshape(gamma2',[],1);

%% Second Vertex
lmb1 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],1),len,2);
lmb2 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],2),len,2);
lmb3 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
resultant3 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma3 = [b1(:,2), b2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*3)-2],1,3) + repmat([0:2],len,1)]',[],1);
circ3 = zeros(len, valWNELE*3);
circ3(sub2ind(size(circ3),rows,col1)) = reshape(gamma3',[],1);

%% Now for the second (aft) post TE wake element
% Evaluated at the mid-point of each edge which splits two HDVEs

idx = all(matWEATT,2); % All edges that split 2 DVEs
nedg = length(matWEATT(idx,1));

%(x,y) of all three vertices of HDVEs in local coordinates
y1 = reshape(matWPLEX(1,2,matWEATT(idx,:)),nedg,2);
y2 = reshape(matWPLEX(2,2,matWEATT(idx,:)),nedg,2);
y3 = reshape(matWPLEX(3,2,matWEATT(idx,:)),nedg,2);

lmb1 = reshape(lamb_circ(matWELOC(idx,:),1),nedg,2);
lmb2 = reshape(lamb_circ(matWELOC(idx,:),2),nedg,2);
lmb3 = reshape(lamb_circ(matWELOC(idx,:),3),nedg,2);

b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(nedg,2);

resultant4 = zeros(nedg,1);

gamma1 = [b1(:,1),b2(:,1),c3(:,1)];
gamma2 = [b1(:,2),b2(:,2),c3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:nedg]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matWEATT(idx,1).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);
col2 = reshape([repmat([(matWEATT(idx,2).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);

circ4 = zeros(nedg, valWNELE*3);

circ4(sub2ind(size(circ4),rows,col1)) = reshape(gamma1',[],1);
circ4(sub2ind(size(circ4),rows,col2)) = reshape(gamma2',[],1);

%% At vertices
clear vnuma vnumb

[ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnuma(:,1) = tb(rb); % Sorting it according to row, cause find returns them jumbled up

[ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnuma(:,2) = tb(rb);

[ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnumb(:,1) = tb(rb);

[ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnumb(:,2) = tb(rb);


%% Vertex 1

% Evaluated at the mid-point of each edge which splits two HDVEs

lmb1 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],1),nedg,2);
lmb2 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],2),nedg,2);
lmb3 = reshape(lamb_vort([vnuma(:,1) vnumb(:,1)],3),nedg,2);

b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(nedg,2);

resultant5 = zeros(nedg,1);

gamma1 = [b1(:,1),b2(:,1),c3(:,1)];
gamma2 = [b1(:,2),b2(:,2),c3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:nedg]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matWEATT(idx,1).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);
col2 = reshape([repmat([(matWEATT(idx,2).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);

circ5 = zeros(nedg, valWNELE*3);

circ5(sub2ind(size(circ5),rows,col1)) = reshape(gamma1',[],1);
circ5(sub2ind(size(circ5),rows,col2)) = reshape(gamma2',[],1);

%% Vertex 2

% Evaluated at the mid-point of each edge which splits two HDVEs

lmb1 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],1),nedg,2);
lmb2 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],2),nedg,2);
lmb3 = reshape(lamb_vort([vnuma(:,2) vnumb(:,2)],3),nedg,2);

b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(nedg,2);

resultant6 = zeros(nedg,1);

gamma1 = [b1(:,1),b2(:,1),c3(:,1)];
gamma2 = [b1(:,2),b2(:,2),c3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:nedg]',1,3)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matWEATT(idx,1).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);
col2 = reshape([repmat([(matWEATT(idx,2).*3)-2],1,3) + repmat([0:2],nedg,1)]',[],1);

circ6 = zeros(nedg, valWNELE*3);

circ6(sub2ind(size(circ6),rows,col1)) = reshape(gamma1',[],1);
circ6(sub2ind(size(circ6),rows,col2)) = reshape(gamma2',[],1);

%% 

d_wake = [circ1; circ2; circ3; circ4; circ5; circ6];
res_wake = [resultant1; resultant2; resultant3; resultant4; resultant5; resultant6];

matNEWWAKECOEFF = d_wake\res_wake;

matNEWWAKECOEFF = reshape(matNEWWAKECOEFF,3,valWNELE,1)';


end

