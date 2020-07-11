function [ CP, LE_Left, LE_Right, TE_Left, TE_Right, vecULS] = fcnPANEL2DVE(panel4corners, i, vecN, vecM, airfoil)
%fcnPANEL2DVE Summary of this function goes here
%   fcnPANEL2DVE takes four corners of a panel and outputs vertices of non-planer DVEs

panelX = [panel4corners([1;4],1),panel4corners([2;3],1)];
panelY = [panel4corners([1;4],2),panel4corners([2;3],2)];
panelZ = [panel4corners([1;4],3),panel4corners([2;3],3)];
%     panelTE = [panel(end,:); panel(end-1,:)];

% NChordwise and NSpanwise
% Generate extra points to find control point
N = vecN(i)*2;
M = vecM(i)*2;

% Chordwise spacing
xspacing = linspace(0,1,M+1);


% split Root Chord to chordwise elements
chordX(:,1) = (panelX(2,1) - panelX(1,1)).*xspacing + panelX(1,1);
chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
chordZ(:,1) = (panelZ(2,1) - panelZ(1,1)).*xspacing + panelZ(1,1);

if ~isempty(airfoil)
    v_norm = [chordZ(1,1) - chordZ(end,1) chordX(end,1) - chordX(1,1)];
    v_norm = v_norm./sqrt(sum(v_norm.^2,2));
    v_dir = [chordX(2:end,1) - chordX(1:end-1,1) chordZ(2:end,1) - chordZ(1:end-1,1)];
    v_dir = cumsum(sqrt(sum(v_dir.^2,2)));
    chord = v_dir(end);
    v_dir = v_dir./max(v_dir);
    
    tmp_us = fcnLINSPLINE(airfoil(i,1).coord(1:airfoil(i,1).le_idx,1), airfoil(i,1).coord(1:airfoil(i,1).le_idx,2));
    tmp_ls = fcnLINSPLINE(airfoil(i,1).coord(airfoil(i,1).le_idx+1:end,1), airfoil(i,1).coord(airfoil(i,1).le_idx+1:end,2));
    
    tmp_cmbr = ((tmp_us(v_dir) + tmp_ls(v_dir))./2).*chord;
    chordX(2:end,1) = chordX(2:end,1) + tmp_cmbr.*v_norm(:,1);
    chordZ(2:end,1) = chordZ(2:end,1) + tmp_cmbr.*v_norm(:,2);
end


% split Tip Chord to chordwise elements
chordX(:,2) = (panelX(2,2) - panelX(1,2)).*xspacing + panelX(1,2);
chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
chordZ(:,2) = (panelZ(2,2) - panelZ(1,2)).*xspacing + panelZ(1,2);

if ~isempty(airfoil)
    v_norm = [chordZ(1,2) - chordZ(end,2) chordX(end,2) - chordX(1,2)];
    v_norm = v_norm./sqrt(sum(v_norm.^2,2));
    v_dir = [chordX(2:end,2) - chordX(1:end-1,2) chordZ(2:end,2) - chordZ(1:end-1,2)];
    v_dir = cumsum(sqrt(sum(v_dir.^2,2)));
    chord = v_dir(end);
    v_dir = v_dir./max(v_dir);
    
    tmp_us = fcnLINSPLINE(airfoil(i,2).coord(1:airfoil(i,2).le_idx,1), airfoil(i,2).coord(1:airfoil(i,2).le_idx,2));
    tmp_ls = fcnLINSPLINE(airfoil(i,2).coord(airfoil(i,2).le_idx+1:end,1), airfoil(i,2).coord(airfoil(i,2).le_idx+1:end,2));
    
    tmp_cmbr = ((tmp_us(v_dir) + tmp_ls(v_dir))./2).*chord;
    chordX(2:end,2) = chordX(2:end,2) + tmp_cmbr.*v_norm(:,1);
    chordZ(2:end,2) = chordZ(2:end,2) + tmp_cmbr.*v_norm(:,2);
end

% split max camber location to chordwise elements
% epsilon = linspace(repsilon, tepsilon, N+1)';

% split chord into chordwise elements
% chord = linspace(rchord, tchord, N+1)';

% Spanwise spacing
yspacing = linspace(0,1,N+1);

% Preallocate the memories for point matrices
PX2 = zeros(M+1,N+1);
PY2 = zeros(M+1,N+1);
PZ2 = zeros(M+1,N+1);

for j = 1:M+1
    if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
        spanwise = (chordY(1,2) - chordY(1,1)).*yspacing + chordY(1,1);
        chordbase = chordY(1,:);
    else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
        spanwise = (chordZ(1,2) - chordZ(1,1)).*yspacing + chordZ(1,1);
        chordbase = chordZ(1,:);
    end
    PX2(j,:) = interp1(chordbase,chordX(j,:),spanwise); % X coordinates of all points
    PY2(j,:) = interp1(chordbase,chordY(j,:),spanwise); % Y coordinates of all points
    PZ2(j,:) = interp1(chordbase,chordZ(j,:),spanwise); % Z coordinates of all points
end

% DVE Parameters Calculation
% Calculate Control Points, stored in 3D matrix
CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],vecM(i),vecN(i),3);
vecULS = ones(size(CP,1), size(CP,2),1);
%     CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],vecM(i),vecN(i),3);
%LE_Mid no longer needed (Jan6,2017 -Alton)
%LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],vecM(i),vecN(i),3);
TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],vecM(i),vecN(i),3);
TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],vecM(i),vecN(i),3);

% Filter Leading Edge Point
LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],vecM(i),vecN(i),3);
LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],vecM(i),vecN(i),3);


end

