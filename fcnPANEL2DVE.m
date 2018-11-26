function [ CP, LE_Left, LE_Right, TE_Left, TE_Right ] = fcnPANEL2DVE( panel4corners, i, vecN, vecM, rcamber, tcamber, repsilon, tepsilon, rchord, tchord)
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
    
    % split Root Chord to chordwise elements
    chordX(:,1) = linspace(panelX(1,1),panelX(2,1),M+1)';
    chordY(:,1) = linspace(panelY(1,1),panelY(2,1),M+1)';
    chordZ(:,1) = linspace(panelZ(1,1),panelZ(2,1),M+1)';
    % split Tip Chord to chordwise elements
    chordX(:,2) = linspace(panelX(1,2),panelX(2,2),M+1)';
    chordY(:,2) = linspace(panelY(1,2),panelY(2,2),M+1)';
    chordZ(:,2) = linspace(panelZ(1,2),panelZ(2,2),M+1)';
    
    % split max camber location in spanwise dir
    camberLOC(:,1) = linspace(rcamber(1,1), tcamber(1,1), N+1)';
    camberMAX(:,1) = linspace(rcamber(1,2), tcamber(1,2), N+1)';

    % split max camber location to chordwise elements
    epsilon = linspace(repsilon, tepsilon, N+1)';
    
    % split chord into chordwise elements
    chord = linspace(rchord, tchord, N+1)';
    
    % linspace spanwise elemenets at each chordwise station
    
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
        
        span_dir = panel4corners(2,:) - panel4corners(1,:);
        span_dir = span_dir./norm(span_dir);
        % Adding in camber
        if j > 1
           x = PX2(j,:) - PX2(1,:);
           idx_for = x <= (camberLOC.*chord)';
           idx_aft = x > (camberLOC.*chord)';
           
           % Local airfoil z-translation for camber
           % http://airfoiltools.com/airfoil/naca4digit
           x = x./chord';
           norm_trans(idx_for) = (camberMAX(idx_for)'./(camberLOC(idx_for)'.^2)).*(2.*camberLOC(idx_for)'.*x(idx_for) - x(idx_for).^2);
           norm_trans(idx_aft) = (camberMAX(idx_aft)'./((1 - camberLOC(idx_aft)').^2)).*(1 - 2.*camberLOC(idx_aft)' + 2.*camberLOC(idx_aft)'.*x(idx_aft) - x(idx_aft).^2);
           norm_trans = norm_trans.*chord';
           
           chordwise_dir = [PX2(j,:) - PX2(1,:); PY2(j,:) - PY2(1,:); PZ2(j,:) - PZ2(1,:)]';
           chordwise_dir = chordwise_dir./(sqrt(sum(chordwise_dir.^2,2)));
           
           norm_dir = cross(chordwise_dir, repmat(span_dir, size(chordwise_dir,1),1), 2);
           
           PX2(j,:) = PX2(j,:) + norm_dir(:,1)'.*norm_trans;
           PY2(j,:) = PY2(j,:) + norm_dir(:,2)'.*norm_trans;
           PZ2(j,:) = PZ2(j,:) + norm_dir(:,3)'.*norm_trans;
        end
    end
    
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

