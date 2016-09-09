clc
clear

%%
LE_1 = [5 4 0];
LE_2 = [5 3 0];
C1 = 2;
C2 = 2;
e1 = 0;
e2 = 0;

fpg = [1 4 5]; % field point in global coords
k = 10;
%% ==============================================================================

X1 = LE_1(1);
X2 = LE_2(1);

Y1 = LE_1(2);
Y2 = LE_2(2);

Z1 = LE_1(3);
Z2 = LE_2(3);

n = 1;
m = 1;

% Calculate twist at leading edge with dihedral correctly modelled
paneldata.rLE = [X1,Y1,Z1];
paneldata.tLE = [X2,Y2,Z2];
paneldata.rchord = C1;
paneldata.tchord = C2;

paneldata.repsilon = deg2rad(e1);
paneldata.tepsilon = deg2rad(e2);

% Read panel corners
panel = reshape(paneldesign(paneldata),3,4)';
panelX = [panel([1;4],1),panel([2;3],1)];
panelY = [panel([1;4],2),panel([2;3],2)];
panelZ = [panel([1;4],3),panel([2;3],3)];

panelTE = [panel(end,:); panel(end-1,:)];


% NChordwise and NSpanwise
% Generate extra points to find control point
N = 1*2;
M = 1*2;


% split Root Chord to chordwise elements
chordX(:,1) = linspace(panelX(1,1),panelX(2,1),M+1)';
chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
chordZ(:,1) = linspace(panelZ(1,1),panelZ(2,1),M+1)';
% split Tip Chord to chordwise elements
chordX(:,2) = linspace(panelX(1,2),panelX(2,2),M+1)';
chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
chordZ(:,2) = linspace(panelZ(1,2),panelZ(2,2),M+1)';

% linspace spanwise elemenets at each chordwise station

% Preallocate the memories for point matrices
PX2 = zeros(M+1,N+1);
PY2 = zeros(M+1,N+1);
PZ2 = zeros(M+1,N+1);

for i = 1:M+1
    if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
        spanwise = linspace(chordY(1,1),chordY(1,2),N+1);
        chordbase = chordY(1,:);
    else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
        spanwise = linspace(chordZ(1,1),chordZ(1,2),N+1);
        chordbase = chordZ(1,:);
    end
    PX2(i,:) = interp1(chordbase,chordX(i,:),spanwise); % X coordinates of all points
    PY2(i,:) = interp1(chordbase,chordY(i,:),spanwise); % Y coordinates of all points
    PZ2(i,:) = interp1(chordbase,chordZ(i,:),spanwise); % Z coordinates of all points
end

% DVE Parameters Calculation
% Calculate Control Points, stored in 3D matrix
CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],m,n,3);
CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],m,n,3);
LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],m,n,3);
TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],m,n,3);
TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],m,n,3);

% Filter Leading Edge Point
LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],m,n,3);
LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],m,n,3);

% Create eta vector for full leading edge
% Non-normalized
LE_vec = LE_Right - LE_Left;

% Create half chord xsi vector
% Non-normalized
xsi_vec = LE_Mid - CP;

DVE_norm = normalize3D(cross(LE_vec, xsi_vec, 3));

% Roll in Degrees -arctan ( Y component / Z component of DVC normal vector)
% atan2d is used here
% roll(nu) right wing up positive
nu = -atan2d(DVE_norm(:,:,2),DVE_norm(:,:,3));

% Pitch in Degrees
% arcsin ( X component of DVE normal vector )
epsilon = asind(DVE_norm(:,:,1));

% Yaw in Degrees
% xsi in local with roll picth, yaw set to zero.. but WHY?
xsi_local = glob_star_3D( xsi_vec,nu,epsilon,zeros(m,n) );
% Magnitude of half chord vector
xsi = (xsi_local(:,:,1).^2+xsi_local(:,:,2).^2+xsi_local(:,:,3).^2).^0.5;
psi = atand(xsi_local(:,:,2)./xsi_local(:,:,1));

% Find eta. bring non-normalized LE_vec to local and half the Y component
LE_vec_local = glob_star_3D( LE_vec,nu,epsilon,psi);
eta = LE_vec_local(:,:,2)./2;

% Find Leading Edge Sweep
% arctan(LE X local component/ LE Y local component)
phi_LE = atand(LE_vec_local(:,:,1)./LE_vec_local(:,:,2));

% Find Trailing Edge Sweep
% Project TE Points onto DVE plane
% (TE_Left / TE_Right) (CP)                   (DVE_norm)
% q(x,y,z) TE point | p(a,b,c) Control Point | n(d,e,f) DVE normal
% q_proj = q - dot(q-p,n)*n
TE_Left_proj = TE_Left-repmat(dot(TE_Left-CP,DVE_norm,3),1,1,3).*DVE_norm;
TE_Right_proj = TE_Right-repmat(dot(TE_Right-CP,DVE_norm,3),1,1,3).*DVE_norm;
TE_vec_proj = TE_Right_proj - TE_Left_proj;

% Rotate the Projected TE on DVE to local reference frame
% arctan(Projected TE local X component/Projected TE local Y component)
TE_vec_proj_local = glob_star_3D( TE_vec_proj,nu,epsilon,psi );
phi_TE = atand(TE_vec_proj_local(:,:,1)./TE_vec_proj_local(:,:,2));

count = 1;

eta = reshape(eta',count,1);%eta(:);
xsi = reshape(xsi',count,1);%xsi(:);
roll = reshape(nu',count,1);%nu(:);
pitch = reshape(epsilon',count,1);%epsilon(:);
yaw = reshape(psi',count,1);%psi(:);
xo = reshape(permute(CP, [2 1 3]),count,3);%reshape(CP(:),count,3);
phiLE = reshape(phi_LE',count,1);%phi_LE(:);
phiTE = reshape(phi_TE',count,1);%phi_TE(:);
Dnorm = reshape(permute(DVE_norm, [2 1 3]),count,3);%reshape(DVE_norm(:),count,3);
TECoordL = reshape(permute(TE_Left_proj, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
TECoordR = reshape(permute(TE_Right_proj, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);
LECoordL = reshape(permute(LE_Left, [2 1 3]),count,3);%reshape(TE_Left_proj(:),count,3);
LECoordR = reshape(permute(LE_Right, [2 1 3]),count,3);%reshape(TE_Right_proj(:),count,3);

DVE_type = 3; % semi-infinite sheet
Temp.DBL_EPS = 1e-14;

[a3x, b3x, c3x] = fcnDVEInduction(Temp, fpg, xo, roll, pitch, phiLE, phiTE, yaw, eta, xsi, DVE_type, k)


%%
% endpoints(:,:,1) = [-0.5 -0.5 0];
% endpoints(:,:,2) = [-0.5 0.5 0];

VLST = [LE_1; LE_2; (TECoordL)];
DVE = [1 2 3];
DNORM = cross(VLST(2,:) - VLST(3,:), VLST(2,:) - VLST(1,:));
DNORM = -DNORM./norm(DNORM);

[PLEX, DVECT] = fcnTRITOLEX(VLST, DNORM);

fpl = fcnTOLOC(1, fpg, DVE, DVECT, VLST, DNORM);

endpoints(:,:,1) = PLEX(1,:);
endpoints(:,:,2) = PLEX(2,:);

% left to right sheet psi = 90 degrees??

[al, bl, cl] = fcnVSIND(endpoints, deg2rad(phiLE), deg2rad(90), fpl, k);

% alpha_r = acos(dot(DVECT(:,:,1),[1 0 0])./(sqrt(sum(abs(DVECT(:,:,1)).^2,2)).*sqrt(sum(abs([1 0 0]).^2,2))));
% beta_r = acos(dot(DVECT(:,:,2),[0 1 0])./(sqrt(sum(abs(DVECT(:,:,2)).^2,2)).*sqrt(sum(abs([0 1 0]).^2,2))));
% gamma_r = acos(dot(DVECT(:,:,3),[0 0 1])./(sqrt(sum(abs(DVECT(:,:,3)).^2,2)).*sqrt(sum(abs([0 0 1]).^2,2))));
% 
% R_alpha = [1 0 0; 0 cos(alpha_r) sin(alpha_r); 0 -sin(alpha_r) cos(alpha_r)];
% R_beta = [cos(beta_r) 0 -sin(beta_r); 0 1 0; sin(beta_r) 0 cos(beta_r)];
% R_gamma = [cos(gamma_r) sin(gamma_r) 0; -sin(gamma_r) cos(gamma_r) 0; 0 0 1];
% Rot = R_alpha.*R_beta.*R_gamma
% 
% b = bl*Rot
% c = cl*Rot

bc = fcnROTVECT([1 1], [bl; cl], DVECT)
% bc = fcnROTVECT(1, bl, DVECT)


% b = fcnTOGLOB(1, bl, DVE, DVECT, VLST)
% fpg2 = fcnTOGLOB(1, fpl, DVE, DVECT, VLST)

%%
% % Plotting global and local to visualize
test_num = 1;
hFig2 = figure(2);
clf(2);
subplot(2,1,1)
patch(PLEX(:,1,test_num),PLEX(:,2,test_num),'b')
hold on
scatter3(fpl(1), fpl(2), fpl(3),100,'ok','filled')
hold off
alpha(0.5);
xlabel('eta-direction','FontSize',15);
ylabel('xi-direction','FontSize',15);
axis equal
grid on
box on
subplot(2,1,2)
patch(VLST(DVE(test_num,:,1),1),VLST(DVE(test_num,:,1),2),VLST(DVE(test_num,:,1),3),'r')
hold on
scatter3(fpg(1), fpg(2), fpg(3),100,'ok','filled')
verbose = 0;
if verbose == 1
    for ii = 1:NELE
        str = sprintf('%d',ii);
        text(DVE(ii,1,3),DVE(ii,2,3),DVE(ii,3,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(VLST(:,1))
        str = sprintf('%d',ii);
        text(VLST(ii,1),VLST(ii,2),VLST(ii,3),str,'Color','g','FontSize',20);
    end
    
    edge1 = VLST(ELST(:,1),:);
    edge2 = VLST(ELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end
    
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,1), DVECT(:,2,1), DVECT(:,3,1), 0.25, 'b') % eta
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,2), DVECT(:,2,2), DVECT(:,3,2), 0.25, 'k') % xi
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,3), DVECT(:,2,3), DVECT(:,3,3), 0.25,'m') % zeta (normal)
    
end

hold off
alpha(0.5);
xlabel('X-direction','FontSize',15);
ylabel('Y-direction','FontSize',15);
zlabel('Z-direction','FontSize',15);
axis equal
grid on
box on





