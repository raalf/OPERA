function [D, R] = fcnDWING(EATT, PLEX, NELE, ELOC, ALIGN, VATT)

[~,~,depth] = size(PLEX);

lamb_circ = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

lamb_vort(:,:,1) = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1; ...
    ];
lamb_vort(:,:,2) = [ ...
    0 1 0; ...
    0 0 1; ...
    1 0 0; ...
    ];

% D = sparse(N*5, N*5);
D = zeros(NELE*5, NELE*5);

idx = all(EATT,2); % All edges that split 2 DVEs
nedg = length(EATT(idx,1));

x1 = reshape(PLEX(1,1,EATT(idx,:)),nedg,2);
x2 = reshape(PLEX(2,1,EATT(idx,:)),nedg,2);
x3 = reshape(PLEX(3,1,EATT(idx,:)),nedg,2);
y1 = reshape(PLEX(1,2,EATT(idx,:)),nedg,2);
y2 = reshape(PLEX(2,2,EATT(idx,:)),nedg,2);
y3 = reshape(PLEX(3,2,EATT(idx,:)),nedg,2);

%% Circulation equations
% Evaluated at the mid-point of each edge which splits two HDVEs

lmb1 = reshape(lamb_circ(ELOC(idx,:),1),nedg,2);
lmb2 = reshape(lamb_circ(ELOC(idx,:),2),nedg,2);
lmb3 = reshape(lamb_circ(ELOC(idx,:),3),nedg,2);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(nedg,2);

gamma1 = [a1(:,1),a2(:,1),b1(:,1),b2(:,1),c3(:,1)];
gamma2 = [a1(:,2),a2(:,2),b1(:,2),b2(:,2),c3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = repmat([1:nedg]',1,5);
rows = reshape(rows',[],1);

% Column indices for the first half of each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = repmat([(EATT(idx,1).*5)-4],1,5)+repmat([0:4],nedg,1);
col1 = reshape(col1',[],1);

% Assigning the values to D using linear indices
idx1 =  sub2ind(size(D),rows,col1);
D(idx1) = reshape(gamma1',[],1);

% Repeat of above, for the 2nd half of the circulation equations
col2 = repmat([(EATT(idx,2).*5)-4],1,5)+repmat([0:4],nedg,1);
col2 = reshape(col2',[],1);
idx2 =  sub2ind(size(D),rows,col2);
D(idx2) = reshape(gamma2',[],1);

%% Vorticity at one vertex
% For reference, I will call 'mm' the current DVE and 'nn' the adjacent DVE 

% Vertices of mm and nn
lmb11 = reshape(lamb_vort(ELOC(idx,:),1,1),nedg,2); % one vertex from side 1
lmb12 = reshape(lamb_vort(ELOC(idx,:),1,2),nedg,2); % second vertex from side 1
lmb21 = reshape(lamb_vort(ELOC(idx,:),2,1),nedg,2); % one vertex from side 2
lmb22 = reshape(lamb_vort(ELOC(idx,:),2,2),nedg,2); % second vertex from side 2
lmb31 = reshape(lamb_vort(ELOC(idx,:),3,1),nedg,2);
lmb32 = reshape(lamb_vort(ELOC(idx,:),3,2),nedg,2);

% Eta direction of mm and nn, an(:,1) is DVE mm and an(:,2) is DVE nn
a2 = ones(nedg,2);
a1 = (lmb11.*x1+lmb21.*x2+lmb31.*x3);

% Xi direction of nn
b2 = ones(nedg,2);
b2(:,1) = zeros(nedg,1);
b1 = (lmb11.*y1+lmb21.*y2+lmb31.*y3);
b1(:,1) = zeros(nedg,1);

c3 = zeros(nedg,2);

dgamma1 = [a1(:,1), a2(:,1), zeros(nedg,1) zeros(nedg,1) zeros(nedg,1)]; % DVE mm in the eta direction

% DVE nn in a combination of the eta and xi directions, multiplied by ALIGN matrix to get the proportions right
% c3 is not multiplied by ALIGN, as it's a constant which is applied to both directions
dgamma2 = [a1(:,2).*ALIGN(idx,1,1),a2(:,2).*ALIGN(idx,1,1),b1(:,2).*ALIGN(idx,2,1),b2(:,2).*ALIGN(idx,2,1),c3(:,2)].*-1;

% Row indices of the rows where vorticity equations will go
rows = repmat([1:nedg]',1,5)+nedg;
rows = reshape(rows',[],1);

% Column indices for the first half of each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = repmat([(EATT(idx,1).*5)-4],1,5)+repmat([0:4],nedg,1);
col1 = reshape(col1',[],1);

% Column indices for the second half of each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col2 = repmat([(EATT(idx,2).*5)-4],1,5)+repmat([0:4],nedg,1);
col2 = reshape(col2',[],1);

% Assigning the values to D using linear indices
idx4 =  sub2ind(size(D),rows,col1);
D(idx4) = reshape(dgamma1',[],1);

idx5 = sub2ind(size(D),rows,col2);
D(idx5) = reshape(dgamma2',[],1);

%% Vorticity at second vertex

% Eta direction of mm and nn, an(:,1) is DVE mm and an(:,2) is DVE nn
a2 = ones(nedg,2);
a1 = (lmb12.*x1+lmb22.*x2+lmb32.*x3);

% Xi direction of nn
b2 = ones(nedg,2);
b2(:,1) = zeros(nedg,1);
b1 = (lmb12.*y1+lmb22.*y2+lmb32.*y3);
b1(:,1) = zeros(nedg,1);

c3 = zeros(nedg,2);

dgamma1 = [a1(:,1), a2(:,1), zeros(nedg,1) zeros(nedg,1) zeros(nedg,1)]; % DVE mm in the eta direction

% DVE nn in a combination of the eta and xi directions, multiplied by ALIGN matrix to get the proportions right
% c3 is not multiplied by ALIGN, as it's a constant which is applied to both directions
dgamma2 = [a1(:,2).*ALIGN(idx,1,1),a2(:,2).*ALIGN(idx,1,1),b1(:,2).*ALIGN(idx,2,1),b2(:,2).*ALIGN(idx,2,1),c3(:,2)].*-1;

% Row indices of the rows where vorticity equations will go
rows = repmat([1:nedg]',1,5)+nedg*2;
rows = reshape(rows',[],1);

% Column indices for the first half of each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = repmat([(EATT(idx,1).*5)-4],1,5)+repmat([0:4],nedg,1);
col1 = reshape(col1',[],1);

% Column indices for the second half of each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col2 = repmat([(EATT(idx,2).*5)-4],1,5)+repmat([0:4],nedg,1);
col2 = reshape(col2',[],1);

% Assigning the values to D using linear indices
idx6 =  sub2ind(size(D),rows,col1);
D(idx6) = reshape(dgamma1',[],1);

idx7 = sub2ind(size(D),rows,col2);
D(idx7) = reshape(dgamma2',[],1);

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all verticies, as well as control points
% Where a vertex is shared by multiple elements, the average of the normals is taken



%% Resultant

R = zeros(NELE*5,1);


end

