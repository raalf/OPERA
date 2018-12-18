function [ CP, LE_Left, LE_Right, TE_Left, TE_Right ] = fcnPANEL2DVE(strATYPE, panel4corners, i, vecN, vecM, repsilon, tepsilon, rchord, tchord, airfoil, strSPACING)
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

if strcmpi(strSPACING,'COSINE')
    spacing = ((1 - (cos(linspace(0,pi,M+1))))./2);
elseif strcmpi(strSPACING,'HALFCOSINE')
    spacing = (1 - ((cos(linspace(0,pi/2,M+1)))));
else
    spacing = linspace(0,1,M+1);
end

% split Root Chord to chordwise elements
chordX(:,1) = linspace(panelX(1,1),panelX(2,1),M+1)';
chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
if strcmpi(strATYPE{2},'PANEL')
    chordX(:,1) = (panelX(2,1) - panelX(1,1)).*spacing + panelX(1,1);
    % Upper surface Z-coordinates
    tmp = fcnLINSPLINE(airfoil(i,1).coord(1:airfoil(i,1).le_idx,1), airfoil(i,1).coord(1:airfoil(i,1).le_idx,2));
    chordZ(:,1,1) = tmp(chordX(:,1));
    % Lower surface Z-coordinates
    tmp = fcnLINSPLINE(airfoil(i,1).coord(airfoil(i,1).le_idx:end,1), airfoil(i,1).coord(airfoil(i,1).le_idx:end,2));
    chordZ(:,1,2) = tmp(chordX(:,1));
    chordZ(end,1,:) = mean(chordZ(end,1,:)); % Closing TE
else
    if isempty(airfoil)
        chordZ(:,1) = linspace(panelZ(1,1),panelZ(2,1),M+1)';
    else
        % Cosine spacing in the chordwise direction to improve LE resolution
        chordX(:,1) = (panelX(2,1) - panelX(1,1)).*spacing + panelX(1,1);
        % Upper surface Z-coordinates
        tmp = fcnLINSPLINE(airfoil(i,1).coord(1:airfoil(i,1).le_idx,1), airfoil(i,1).coord(1:airfoil(i,1).le_idx,2));
        tmpZ(:,1,1) = tmp(chordX(:,1));
        % Lower surface Z-coordinates
        tmp = fcnLINSPLINE(airfoil(i,1).coord(airfoil(i,1).le_idx:end,1), airfoil(i,1).coord(airfoil(i,1).le_idx:end,2));
        tmpZ(:,1,2) = tmp(chordX(:,1));
        chordZ(:,1) = mean(tmpZ(:,1,:),3);
    end
end
tmpZ = [];

% split Tip Chord to chordwise elements
chordX(:,2) = linspace(panelX(1,2),panelX(2,2),M+1)';
chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
if strcmpi(strATYPE{2},'PANEL')
    % Cosine spacing in the chordwise direction to improve LE resolution
    chordX(:,2) = (panelX(2,2) - panelX(1,2)).*spacing + panelX(1,2);
    % Upper surface Z-coordinates
    tmp = fcnLINSPLINE(airfoil(i,2).coord(1:airfoil(i,2).le_idx,1), airfoil(i,2).coord(1:airfoil(i,2).le_idx,2));
    chordZ(:,2,1) = tmp(chordX(:,2));
    % Lower surface Z-coordinates
    tmp = fcnLINSPLINE(airfoil(i,2).coord(airfoil(i,2).le_idx:end,1), airfoil(i,2).coord(airfoil(i,2).le_idx:end,2));
    chordZ(:,2,2) = tmp(chordX(:,2));
    chordZ(end,2,:) = mean(chordZ(end,2,:)); % Closing TE
else
    if isempty(airfoil)
        chordZ(:,2) = linspace(panelZ(1,2),panelZ(2,2),M+1)';
    else
        % Cosine spacing in the chordwise direction to improve LE resolution
        chordX(:,2) = (panelX(2,2) - panelX(1,2)).*spacing + panelX(1,2);
        % Upper surface Z-coordinates
        tmp = fcnLINSPLINE(airfoil(i,2).coord(1:airfoil(i,2).le_idx,1), airfoil(i,2).coord(1:airfoil(i,2).le_idx,2));
        tmpZ(:,2,1) = tmp(chordX(:,2));
        % Lower surface Z-coordinates
        tmp = fcnLINSPLINE(airfoil(i,2).coord(airfoil(i,2).le_idx:end,1), airfoil(i,2).coord(airfoil(i,2).le_idx:end,2));
        tmpZ(:,2,2) = tmp(chordX(:,2));
        chordZ(:,2) = mean(tmpZ(:,2,:),3);
    end
end

% split max camber location to chordwise elements
epsilon = linspace(repsilon, tepsilon, N+1)';

% split chord into chordwise elements
chord = linspace(rchord, tchord, N+1)';

% linspace spanwise elemenets at each chordwise station

if strcmpi(strATYPE{2},'THIN')
    % Preallocate the memories for point matrices
    PX2 = zeros(M+1,N+1);
    PY2 = zeros(M+1,N+1);
    PZ2 = zeros(M+1,N+1);
    
    for j = 1:M+1
        if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
            spanwise = linspace(chordY(1,1),chordY(1,2),N+1);
            chordbase = chordY(1,:);
        else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
            spanwise = linspace(chordZ(1,1),chordZ(1,2),N+1);
            chordbase = chordZ(1,:);
        end
        PX2(j,:) = interp1(chordbase,chordX(j,:),spanwise); % X coordinates of all points
        PY2(j,:) = interp1(chordbase,chordY(j,:),spanwise); % Y coordinates of all points
        PZ2(j,:) = interp1(chordbase,chordZ(j,:),spanwise); % Z coordinates of all points
    end
    
elseif strcmpi(strATYPE{2},'PANEL')
    % Preallocate the memories for point matrices
    PX2 = zeros(M+1,N+1,2);
    PY2 = zeros(M+1,N+1,2);
    PZ2 = zeros(M+1,N+1,2);
    
    for jj = 1:2
        for j = M+1:-1:1
            if diff(chordY(1,:)) ~= 0    %if: panel is NOT vertical
                spanwise = linspace(chordY(1,1),chordY(1,2),N+1);
                chordbase = chordY(1,:);
            else   %else: panel is vertical (SPECIAL CASE) (difference of Y coordinates is zero)
                spanwise = linspace(chordZ(1,1),chordZ(1,2),N+1);
                chordbase = chordZ(1,:);
            end
            PX2(j,:,jj) = interp1(chordbase,chordX(j,:),spanwise); % X coordinates of all points
            PY2(j,:,jj) = interp1(chordbase,chordY(j,:),spanwise); % Y coordinates of all points
            PZ2(j,:,jj) = interp1(chordbase,chordZ(j,:,jj),spanwise); % Z coordinates of all points
                        
%             hFig4 = figure(4);
%             clf(4);
%             scatter3(airfoil(i,1).coord(:,1), repmat(PY2(j,1,jj),size(airfoil(i,1).coord(:,2),1),1), airfoil(i,1).coord(:,2),'ok')
%             hold on
%             scatter3(airfoil(i,2).coord(:,1), repmat(PY2(j,end,jj),size(airfoil(i,2).coord(:,2),1),1), airfoil(i,2).coord(:,2),'^b')
%             scatter3(PX2(j,:,jj), PY2(j,:,jj), PZ2(j,:,jj), 'xr')
%             box on
%             grid minor
%             box on
%             axis equal
            
        end
    end
    tmp = flipud(PX2(:,:,1));
    PX2 = [tmp(1:end-1,:); PX2(:,:,2)];
    
    tmp = flipud(PY2(:,:,1));
    PY2 = [tmp(1:end-1,:); PY2(:,:,2)];
    
    tmp = flipud(PZ2(:,:,1));
    PZ2 = [tmp(1:end-1,:); PZ2(:,:,2)];
    
end

if strcmpi(strATYPE{2},'PANEL')
    % DVE Parameters Calculation
    % Calculate Control Points, stored in 3D matrix
    CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],(vecM(i)*2),vecN(i),3);
    %     CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],vecM(i),vecN(i),3);
    %LE_Mid no longer needed (Jan6,2017 -Alton)
    %LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],vecM(i),vecN(i),3);
    TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],(vecM(i)*2),vecN(i),3);
    TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],(vecM(i)*2),vecN(i),3);
    
    % Filter Leading Edge Point
    LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],(vecM(i)*2),vecN(i),3);
    LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],(vecM(i)*2),vecN(i),3);
    
else
    % DVE Parameters Calculation
    % Calculate Control Points, stored in 3D matrix
    CP = reshape([PX2(2:2:end,2:2:end),PY2(2:2:end,2:2:end),PZ2(2:2:end,2:2:end)],vecM(i),vecN(i),3);
    %     CP_Right = reshape([PX2(2:2:end,3:2:end),PY2(2:2:end,3:2:end),PZ2(2:2:end,3:2:end)],vecM(i),vecN(i),3);
    %LE_Mid no longer needed (Jan6,2017 -Alton)
    %LE_Mid = reshape([PX2(1:2:end-1,2:2:end),PY2(1:2:end-1,2:2:end),PZ2(1:2:end-1,2:2:end)],vecM(i),vecN(i),3);
    TE_Right = reshape([PX2(3:2:end,3:2:end),PY2(3:2:end,3:2:end),PZ2(3:2:end,3:2:end)],vecM(i),vecN(i),3);
    TE_Left = reshape([PX2(3:2:end,1:2:end-2),PY2(3:2:end,1:2:end-2),PZ2(3:2:end,1:2:end-2)],vecM(i),vecN(i),3);
    
    % Filter Leading Edge Point
    LE_Left = reshape([PX2(1:2:end-2,1:2:end-2),PY2(1:2:end-2,1:2:end-2),PZ2(1:2:end-2,1:2:end-2)],vecM(i),vecN(i),3);
    LE_Right = reshape([PX2(1:2:end-2,3:2:end),PY2(1:2:end-2,3:2:end),PZ2(1:2:end-2,3:2:end)],vecM(i),vecN(i),3);
end

end

