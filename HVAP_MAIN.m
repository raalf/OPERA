clear
clc
tic

% Analysis Type and Geometry File

ATYPE = 'LS'; % Lifting Surface
% STL = 'CAD Geom/simple_liftingsurface.stl';
STL = 'CAD Geom/quad.stl';

% STL = 'Cad Geom/lifting_split.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

%% Triangulating Geometry

[TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, DNORM, PLEX, DVECT, ALIGN] = fcnTRIANG(STL, ATYPE);


%% D-Matrix Creation

[~,~,depth] = size(PLEX);

lamb_circ = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

% lamb_vort(:,:,1) = [ ...
%     1 0 0; ...
%     0 1 0; ...
%     0 0 1; ...
%     ];
% lamb_vort(:,:,2) = [ ...
%     0 1 0; ...
%     0 0 1; ...
%     1 0 0; ...
%     ];

% D = sparse(N*6, N*6);
D = zeros(NELE*6, NELE*6);

idx = all(EATT,2); % All edges that split 2 DVEs
nedg = length(EATT(idx,1));

x1 = reshape(PLEX(1,1,EATT(idx,:)),nedg,2);
x2 = reshape(PLEX(2,1,EATT(idx,:)),nedg,2);
x3 = reshape(PLEX(3,1,EATT(idx,:)),nedg,2);
y1 = reshape(PLEX(1,2,EATT(idx,:)),nedg,2);
y2 = reshape(PLEX(2,2,EATT(idx,:)),nedg,2);
y3 = reshape(PLEX(3,2,EATT(idx,:)),nedg,2);

% Circulation 
lmb1 = reshape(lamb_circ(ELOC(idx,:),1),nedg,2);
lmb2 = reshape(lamb_circ(ELOC(idx,:),2),nedg,2);
lmb3 = reshape(lamb_circ(ELOC(idx,:),3),nedg,2);

a3 = ones(nedg,2);
a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b3 = ones(nedg,2);
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

gamma1 = [a1(:,1),a2(:,1),a3(:,1),b1(:,1),b2(:,1),b3(:,1)];
gamma2 = [a1(:,2),a2(:,2),a3(:,2),b1(:,2),b2(:,2),b3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = repmat([1:nedg]',1,6);
rows = reshape(rows',[],1);

% Column indices for the first half of each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = repmat([(EATT(idx,1).*6)-5],1,6)+repmat([0:5],nedg,1);
col1 = reshape(col1',[],1);

% Assigning the values to D using linear indices
idx1 =  sub2ind(size(D),rows,col1);
D(idx1) = reshape(gamma1',[],1);

% Repeat of above, for the 2nd half of the circulation equations
col2 = repmat([(EATT(idx,2).*6)-5],1,6)+repmat([0:5],nedg,1);
col2 = reshape(col2',[],1);
idx2 =  sub2ind(size(D),rows,col2);
D(idx2) = reshape(gamma2',[],1);

% % This part is tricky
% % Vorticity needs to be constant at vertices, which means for every edge we need to make it constant at the two vertices of that edge
% % This results in 2 vorticity equations per edge
% % Vorticity has two directions as well, eta and xi, and the eta of one DVE will be a combo of the adjacent DVE's eta and xi
% % So we need to go through DVEs and for every DVE, set eta and xi equal at the vertices, with that being equal to a combo of eta and xi
% % vorticity from the neighbouring DVE using the ADJE matrix
% % For reference, I will call 'mm' the current DVE and 'nn' the adjacent DVE
% 
% % Vertices of mm and nn
% lmb11 = reshape(lamb_vort(ELOC(idx,:),1,1),nedg,2); % one vertex from side 1
% lmb12 = reshape(lamb_vort(ELOC(idx,:),1,2),nedg,2); % second vertex from side 1
% lmb21 = reshape(lamb_vort(ELOC(idx,:),2,1),nedg,2); % one vertex from side 2
% lmb22 = reshape(lamb_vort(ELOC(idx,:),2,2),nedg,2); % second vertex from side 2
% lmb31 = reshape(lamb_vort(ELOC(idx,:),3,1),nedg,2);
% lmb32 = reshape(lamb_vort(ELOC(idx,:),3,2),nedg,2);
% 
% % Eta direction of mm and nn, an(:,1) is DVE mm and an(:,2) is DVE nn
% a3 = zeros(nedg,2);
% a2 = ones(nedg,2);
% a1 = (lmb11.*x1+lmb21.*x2+lmb31.*x3);
% 
% % Xi direction of nn
% b3 = zeros(nedg,2);
% b3(:,1) = zeros(nedg,1);
% b2 = ones(nedg,2);
% b2(:,1) = zeros(nedg,1);
% b1 = (lmb11.*y1+lmb21.*y2+lmb31.*y3);
% b1(:,1) = zeros(nedg,1);
% 
% dgamma1 = [a1(:,1),a2(:,1),a3(:,1)]; % DVE mm in the eta direction
% % DVE nn in a combination of the eta and xi directions, multiplied by ALIGN matrix to get the proportions right
% dgamma2 = [a1(:,2).*ALIGN(idx,1,1),a2(:,2).*ALIGN(idx,1,1),a3(:,2).*ALIGN(idx,1,1),b1(:,2).*ALIGN(idx,2,1),b2(:,2).*ALIGN(idx,2,1),b3(:,2).*ALIGN(idx,2,1)];
% 
% % % Xi direction of mm
% % b3 = zeros(nedg,2);
% % b2 = ones(nedg,2);
% % b1 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);


%% Plot

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT);


%% End

clearvars -except TR ADJE ELST VLST DVE NELE EATT EIDX ELOC DNORM PLEX DVECT ALIGN D

toc

whos
