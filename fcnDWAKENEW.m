function [matNEWWAKECOEFF] = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, vecSPANDIR, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWCOEFF, matWALIGN, matWEIDX)

% This function finds the coefficients of the post-trailing edge elements, by setting the spanwise component equal to the opposite of the circulation along the trailing edge of the wing (spanwise)
% There is no chordwise change in wake coefficients now (steady only)

lamb_mid = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

lamb_vert = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1; ...
    ];

len = length(vecWLEDVE);

aftdves = [valWNELE - valWSIZE + 1:valWNELE]'; % aft row of post-TE wake HDVEs 

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
lmb1 = reshape(lamb_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],1),len,2); % Ones(le,1) because wake LE is always local edge 1
lmb2 = reshape(lamb_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],2),len,2);
lmb3 = reshape(lamb_mid([nonzeros(matELOC(vecTE,:)) ones(len,1)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res1 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma1 = [a1(:,2), a2(:,2), b1(:,2), b2(:,2), c3(:,2)];

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*5)-4],1,5) + repmat([0:4],len,1)]',[],1);
circ1 = zeros(len, valWNELE*5);
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
lmb1 = reshape(lamb_vert([vnuma(:,1) vnumb(:,1)],1),len,2);
lmb2 = reshape(lamb_vert([vnuma(:,1) vnumb(:,1)],2),len,2);
lmb3 = reshape(lamb_vert([vnuma(:,1) vnumb(:,1)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res2 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma2 = [a1(:,2), a2(:,2), b1(:,2), b2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*5)-4],1,5) + repmat([0:4],len,1)]',[],1);
circ2 = zeros(len, valWNELE*5);
circ2(sub2ind(size(circ2),rows,col1)) = reshape(gamma2',[],1);

%% Second Vertex
lmb1 = reshape(lamb_vert([vnuma(:,2) vnumb(:,2)],1),len,2);
lmb2 = reshape(lamb_vert([vnuma(:,2) vnumb(:,2)],2),len,2);
lmb3 = reshape(lamb_vert([vnuma(:,2) vnumb(:,2)],3),len,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(len,2);

% The resultant is the circulation at the trailing edge of the wing, found using the coefficients A1 A2 B1 B2 C3 we solved for in DWING/Resultant
res3 = -sum([a1(:,1).*vecSPANDIR(:,1), a2(:,1).*vecSPANDIR(:,1), b1(:,1).*vecSPANDIR(:,2), b2(:,1).*vecSPANDIR(:,2), c3(:,1)].*matCOEFF(vecTEDVE,:),2);

% These are the wake coefficients (spanwise) that we need to solve for, B1, B2 and C3
gamma3 = [a1(:,2), a2(:,2), b1(:,2), b2(:,2), c3(:,2)]; % Wake only changes in the spanwise direction, which is always local xsi

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:len]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(vecWLEDVE.*5)-4],1,5) + repmat([0:4],len,1)]',[],1);
circ3 = zeros(len, valWNELE*5);
circ3(sub2ind(size(circ3),rows,col1)) = reshape(gamma3',[],1);

%% Now for the second (aft) post TE wake element
% Evaluated at the mid-point of each edge which splits two HDVEs

idx = sort(ismember(matWEATT, vecWLEDVE), 2, 'descend');
idx = idx(:,1);
idx_postte = zeros(size(matWEATT(:,1)));
idx_postte(unique([reshape(unique(matWEIDX(vecWLEDVE,:)),[],1); reshape(unique(matWEIDX(aftdves,2:3)),[],1)])) = 1; % Here, we are going to ignore edges that aren't between the two post-te HDVEs, as we want to ignore the rest of the wake
idx = idx & all(matWEATT,2) & idx_postte; % All post-te row HDVE edges which we must set circulation constant

nedg = length(nonzeros(idx));

circ4 = fcnDCIRC(idx, nedg, lamb_mid, valWNELE, matWPLEX, matWEATT, matWELOC);
res4 = zeros(nedg,1);

%% Vorticity

% along that edge equal for both HDVES
% vnuma is local vertices 1 and 2 (columns) for HDVE 1
% vnumb is local vertices 1 and 2 for HDVE 2
[ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnuma2(:,1) = tb(rb); % Sorting it according to row, cause find returns them jumbled up

[ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnuma2(:,2) = tb(rb);

[ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnumb2(:,1) = tb(rb);

[ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnumb2(:,2) = tb(rb);

% First vertex
[vort_e1, vort_x1] = fcnDVORT(idx, vnuma2(:,1), vnumb2(:,1), nedg, lamb_vert, valWNELE, matWPLEX, matWEATT, matWELOC, matWALIGN);

% Second vertex
[vort_e2, vort_x2] = fcnDVORT(idx, vnuma2(:,2), vnumb2(:,2), nedg, lamb_vert, valWNELE, matWPLEX, matWEATT, matWELOC, matWALIGN);

vort5 = [vort_e1; vort_x1; vort_e2; vort_x2];
res5 = zeros(nedg*4,1);

d_wake = [circ1; circ2; circ3; circ4; vort5];
res_wake = [res1; res2; res3; res4; res5];

%% Back half of post-TE wake DVEs

%{
% The wake HDVEs in question:

% If this is the first row of wake elements (timestep = 1) then we set the aft row of wake elements (2nd triangles in the post-te row) to constant circ
if valWSIZE*2 == valWNELE
      
    cols = reshape([repmat([(aftdves.*5)-4],1,2)+repmat([0:1],length(aftdves),1)]',[],1);
    
    rows_rep = reshape([repmat([1:valWSIZE]',1,length(aftdves)*2)]',[],1); 
    cols_rep = repmat(cols,valWSIZE,1);
    
    d_wake(sub2ind(size(d_wake),rows_rep,cols_rep)) = 0;
    
    rows = reshape([repmat([1:length(aftdves)*2]',1,1)]',[],1); 
    cols1 = (aftdves.*5)-4;
    cols2 = cols1 + 1;
    
    circ_chw = zeros(length(aftdves)*2, valWNELE*5);
    circ_chw(sub2ind(size(circ_chw),rows,[cols1; cols2])) = 1;
    
    res_chw = zeros(length(aftdves)*2,1);
    
    d_wake = [d_wake; circ_chw];
    res_wake = [res_wake; res_chw];
       
% Otherwise we set it equal to the leading edge of the row of wake elements behind it   
else
    aftedges = matWEIDX(aftdves,1); % It will always be edge 1 of these HDVEs, due to the way the wake is generated
    
    idx = logical(zeros(length(matWEATT(:,1)),1));
    idx(aftedges) = 1;

    nedg = length(nonzeros(idx));
         
    % Circulation
%     circ_chw2 = fcnDCIRC(idx, nedg, lamb_mid, valWNELE, matWPLEX, matWEATT, matWELOC);
%     res_chw2 = zeros(nedg,1);
    
    % Vorticity
    % along that edge equal for both HDVES
    % vnuma is local vertices 1 and 2 (columns) for HDVE 1
    % vnumb is local vertices 1 and 2 for HDVE 2
    [ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,1),1,3));
    [~,rb,] = sort(ta);
    vnuma3(:,1) = tb(rb); % Sorting it according to row, cause find returns them jumbled up

    [ta,tb,~] = find(matWDVE(matWEATT(idx,1),:,1) == repmat(matWELST(idx,2),1,3));
    [~,rb,] = sort(ta);
    vnuma3(:,2) = tb(rb);

    [ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,1),1,3));
    [~,rb,] = sort(ta);
    vnumb3(:,1) = tb(rb);

    [ta,tb,~] = find(matWDVE(matWEATT(idx,2),:,1) == repmat(matWELST(idx,2),1,3));
    [~,rb,] = sort(ta);
    vnumb3(:,2) = tb(rb);
    
%     % First vertex
%     [vort_echw1, vort_xchw1] = fcnDVORT(idx, vnuma3(:,1), vnumb3(:,1), nedg, lamb_vert, valWNELE, matWPLEX, matWEATT, matWELOC, matWALIGN);
% 
%     % Second vertex
%     [vort_echw2, vort_xchw2] = fcnDVORT(idx, vnuma3(:,2), vnumb3(:,2), nedg, lamb_vert, valWNELE, matWPLEX, matWEATT, matWELOC, matWALIGN);
%     
%     vort_chw = [vort_echw1; vort_xchw1; vort_echw2; vort_xchw2];
%     res_vchw = zeros(length(vort_chw(:,1)),1);
    
%     d_wake = [d_wake; circ_chw2];
%     res_wake = [res_wake; res_chw2];
%     d_wake(:,1:valWSIZE*2*5) = [];
    d_wake = d_wake(:, end - (valWSIZE*2*5) + 1:end);
    
%     % Adding in values of older (2nd post-te row) of wake HDVEs, where the value of circulation is known. Adding this to resultant
%     oldwakedves = setdiff(unique(matWEATT(aftedges,:)), aftdves); % the older wake HDVEs that border with aftdves
%     
%     rows = reshape([repmat([1:length(oldwakedves)*5]',1,1)]',[],1); 
%     col1 = reshape([repmat([(oldwakedves.*5)-4],1,5) + repmat([0:4],nedg,1)]',[],1);
%     
%     circ_chw3 = zeros(length(oldwakedves)*5, valWNELE*5);
%     circ_chw3(sub2ind(size(circ_chw3),rows,col1)) = 1;
%     
%     res_chw3 = reshape(matWCOEFF(oldwakedves,:)',[],1);
%     
%     
%     res_wake = [res_wake; res_chw2; res_vchw; res_chw3];
    
end
%}

cols = reshape([repmat([(aftdves.*5)-4],1,2)+repmat([0:1],length(aftdves),1)]',[],1);

rows_rep = reshape([repmat([1:valWSIZE]',1,length(aftdves)*2)]',[],1); 
cols_rep = repmat(cols,valWSIZE,1);

d_wake(sub2ind(size(d_wake),rows_rep,cols_rep)) = 0;

rows = reshape([repmat([1:length(aftdves)*2]',1,1)]',[],1); 
cols1 = (aftdves.*5)-4;
cols2 = cols1 + 1;

circ_chw = zeros(length(aftdves)*2, valWNELE*5);
circ_chw(sub2ind(size(circ_chw),rows,[cols1; cols2])) = 1;

res_chw = zeros(length(aftdves)*2,1);

d_wake = [d_wake; circ_chw];
res_wake = [res_wake; res_chw];

d_wake = d_wake(:, end - (valWSIZE*2*5) + 1:end);

matNEWWAKECOEFF = d_wake\res_wake;
matNEWWAKECOEFF = reshape(matNEWWAKECOEFF,5,valWSIZE*2,1)';


end

