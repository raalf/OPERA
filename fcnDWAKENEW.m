function [matNEWWAKECOEFF] = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWEIDX)
% This function finds the coefficients of the post-trailing edge elements, by setting the spanwise component equal to the opposite of the circulation along the trailing edge of the wing (spanwise)

lambda_mid = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

lambda_vert = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1; ...
    ];

len = length(vecWLEDVE);

aftdves = [valWNELE - valWSIZE + 1:valWNELE]'; % aft row of post-TE wake HDVEs 

%% Circulation at TE of wing/LE of wake Prelim
% (x,y) of all three vertices of HDVEs at trailing edge of the wing and the corresponding leading edge of the post te wake HDVE
x1 = [reshape(matPLEX(1,1,vecTEDVE),len,1) reshape(matWPLEX(1,1,vecWLEDVE),len,1)];
x2 = [reshape(matPLEX(2,1,vecTEDVE),len,1) reshape(matWPLEX(2,1,vecWLEDVE),len,1)];
x3 = [reshape(matPLEX(3,1,vecTEDVE),len,1) reshape(matWPLEX(3,1,vecWLEDVE),len,1)];
y1 = [reshape(matPLEX(1,2,vecTEDVE),len,1) reshape(matWPLEX(1,2,vecWLEDVE),len,1)];
y2 = [reshape(matPLEX(2,2,vecTEDVE),len,1) reshape(matWPLEX(2,2,vecWLEDVE),len,1)];
y3 = [reshape(matPLEX(3,2,vecTEDVE),len,1) reshape(matWPLEX(3,2,vecWLEDVE),len,1)];

%% Midpoint of edge to set circulation of wake opposite to wing
lmb1 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],1),len,2); % Ones(le,1) because wake LE is always local edge 1
lmb2 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],2),len,2);
lmb3 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],3),len,2);

a2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
a1 = 0.5.*(a2.^2);
b2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
b1 = 0.5.*(b2.^2);
c2 = a2.*b2;
c3 = ones(len,2);

% res1 = [];
% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res1 = sum([a1(:,1), a2(:,1), b1(:,1), b2(:,1), c2(:,1), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma1 = -[a1(:,2), a2(:,2), b1(:,2), b2(:,2), c2(:,2), c3(:,2)];

% circ1 = [];
circ1 = fcnCREATEDSECT(sparse(len, valWNELE*6), len, 6, vecWLEDVE, [], gamma1, []);

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
lmb1 = reshape(lambda_vert([vnuma(:,1) vnumb(:,1)],1),len,2);
lmb2 = reshape(lambda_vert([vnuma(:,1) vnumb(:,1)],2),len,2);
lmb3 = reshape(lambda_vert([vnuma(:,1) vnumb(:,1)],3),len,2);

a2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
a1 = 0.5.*(a2.^2);
b2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
b1 = 0.5.*(b2.^2);
c2 = a2.*b2;
c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res2 = sum([a1(:,1), a2(:,1), b1(:,1), b2(:,1), c2(:,1), c3(:,1)].*matCOEFF(vecTEDVE,:),2);
% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma2 = -[a1(:,2), a2(:,2), b1(:,2), b2(:,2), c2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

circ2 = fcnCREATEDSECT(sparse(len, valWNELE*6), len, 6, vecWLEDVE, [], gamma2, []);

%% Second Vertex
lmb1 = reshape(lambda_vert([vnuma(:,2) vnumb(:,2)],1),len,2);
lmb2 = reshape(lambda_vert([vnuma(:,2) vnumb(:,2)],2),len,2);
lmb3 = reshape(lambda_vert([vnuma(:,2) vnumb(:,2)],3),len,2);

a2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
a1 = 0.5.*(a2.^2);
b2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
b1 = 0.5.*(b2.^2);
c2 = a2.*b2;
c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res3 = sum([a1(:,1), a2(:,1), b1(:,1), b2(:,1), c2(:,1), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma3 = -[a1(:,2), a2(:,2), b1(:,2), b2(:,2), c2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

circ3 = fcnCREATEDSECT(sparse(len, valWNELE*6), len, 6, vecWLEDVE, [], gamma3, []);

%% Now for the second (aft) post TE wake element
% Evaluated at the mid-point of each edge which splits two HDVEs
idx = sort(ismember(matWEATT, vecWLEDVE), 2, 'descend');
idx = idx(:,1);
idx_postte = zeros(size(matWEATT(:,1)));
idx_postte(unique([reshape(unique(matWEIDX(vecWLEDVE,:)),[],1); reshape(unique(matWEIDX(aftdves,2:3)),[],1)])) = 1; % Here, we are going to ignore edges that aren't between the two post-te HDVEs, as we want to ignore the rest of the wake
idx = idx & all(matWEATT,2) & idx_postte; % All post-te row HDVE edges which we must set circulation constant
nedg = length(nonzeros(idx));

vnuma = [];
vnumb = [];
% Typically found at the two vertices that split two HDVEs
% vnuma is local vertices 1 and 2 (columns) for HDVE 1
% vnumb is local vertices 1 and 2 for HDVE 2
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

circ4 = [];
circ5 = [];
circ6 = [];
circ4 = fcnDCIRC(idx, nedg, lambda_mid, valWNELE, matWPLEX, matWEATT, matWELOC);
circ5 = fcnDCIRC2(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, matWELOC);
circ6 = fcnDCIRC2(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, matWELOC);

res4 = zeros(nedg, 1);
res5 = zeros(nedg, 1);
res6 = zeros(nedg, 1);

%% Vorticity
vort1 = [];
vort2 = [];
vort3 = [];
vort4 = [];

temp = reshape(permute(matWPLEX,[1 3 2]),[],3,1);

e1vec = temp(matWEATT(idx,1).*3 + vnuma(:,2) - 3,:) - temp(matWEATT(idx,1).*3 + vnuma(:,1) - 3,:);
e2vec = temp(matWEATT(idx,2).*3 + vnumb(:,2) - 3,:) - temp(matWEATT(idx,2).*3 + vnumb(:,1) - 3,:);
e1vec = e1vec(:,1:2)./repmat(sqrt(e1vec(:,1).^2 + e1vec(:,2).^2),1,2);
e2vec = e2vec(:,1:2)./repmat(sqrt(e2vec(:,1).^2 + e2vec(:,2).^2),1,2);

vort1 = fcnDVORTEDGE(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, e1vec, e2vec);
vort2 = fcnDVORTEDGE(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, e1vec, e2vec);

vort3 = fcnDVORTEDGE(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, [-e1vec(:,2) e1vec(:,1)], [-e2vec(:,2) e2vec(:,1)]);
vort4 = fcnDVORTEDGE(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valWNELE, matWPLEX, matWEATT, [-e1vec(:,2) e1vec(:,1)], [-e2vec(:,2) e2vec(:,1)]);

resv1 = zeros(nedg, 1);
resv2 = zeros(nedg, 1);
resv3 = zeros(nedg, 1);
resv4 = zeros(nedg, 1);

%%
d_wake = [circ1; circ2; circ3; circ4; circ5; circ6; vort1; vort2; vort3; vort4];
res_wake = [res1; res2; res3; res4; res5; res6; resv1; resv2; resv3; resv4];

matNEWWAKECOEFF = d_wake\res_wake;
matNEWWAKECOEFF = reshape(matNEWWAKECOEFF,6,[],1)';


end

