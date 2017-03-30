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

%% Vorticity

idx = sort(ismember(matWEATT, vecWLEDVE), 2, 'descend');
idx = idx(:,1);
idx_postte = zeros(size(matWEATT(:,1)));
idx_postte(unique([reshape(unique(matWEIDX(vecWLEDVE,:)),[],1); reshape(unique(matWEIDX(aftdves,2:3)),[],1)])) = 1; % Here, we are going to ignore edges that aren't between the two post-te HDVEs, as we want to ignore the rest of the wake
idx = idx & all(matWEATT,2) & idx_postte; % All post-te row HDVE edges which we must set circulation constant

nedg = length(nonzeros(idx));

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

d_wake = vort5;
res_wake = res5;

%% Back half of post-TE wake DVEs

matNEWWAKECOEFF = d_wake\res_wake;
matNEWWAKECOEFF = reshape(matNEWWAKECOEFF,4,valWSIZE*2,1)';


end

