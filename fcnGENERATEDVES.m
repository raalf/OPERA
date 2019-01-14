function [matPOINTS, matTEPOINTS, matLEPOINTS, vecULS] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM, vecPANELTE, vecPANELLE, strATYPE, strAIRFOIL, strSPACING)

% This function has been taken from VAP2 and repurposed for OPERA. It takes panel corner points and discretizes
% the panel into quadralateral DVEs, which I then subdivide into triangles to match STL format.

%   V0 - before fixing spanwise interp
%   V1 - fixed vertical panel (90deg dihedral)
%      - Reprogrammed the DVE interpolation method
%  1.1 - Preallocate the memory for point matrices, do you know how memory is accessed?
%   V2 - Rework the Leading Edge Vector
%      - Caculate the Yaw angle of the DVE by assuming yaw=0 to rotate the xsi vector
%      - Rotate DVE normal vector by local roll, pitch, and yaw using 'glob_star_3D'
%      - Comptue LE Sweep
%      - Project TE to DVE, Rotate adn Comptue TE Sweep
%   V3 - Function overhaul for VAP2.0
%  3.5 - Modify non-planer VLST (Jan 6, 2017)
%        old matNPVLST will be not called matNTVLST to specify it holds dve
%        infomation of non-twisted wing
%      - new matNPVLST will now hold non-planer dve coordinates of
%        non-modified wing geometry specified in input file
%
% Fixed how DVEs matrix is converted from 2D grid to 1D array. 16/01/2016 (Alton)

% INPUT:
%   valPANELS - number of wing panels
%   matGEOM - 2 x 5 x valPANELS matrix, with (x,y,z) coords of edge points, and chord and twist at each edge
%   vecSYM - valPANELS x 1 vector of 0, 1, or 2 which denotes the panels with symmetry
%   vecN - valPANELS x 1 vector of spanwise elements per DVE
%   vecM - valPANELS x 1 vector of chordwise elements per DVE

% OUTPUT: (ALL OUTPUT ANGLES ARE IN RADIAN)
%   matCENTER - valNELE x 3 matrix of (x,y,z) locations of DVE control points
%   vecDVEHVSPN - valNELE x 1 vector of DVE half spans
%   vecDVEHVCRD - valNELE x 1 vector of DVE half chords
%   vecDVELESWP - valNELE x 1 vector of DVE leading edge sweep (radians)
%   vecDVEMCSWP - valNELE x 1 vector of DVE mid-chord sweep (radians)
%   vecDVETESWP - valNELE x 1 vector of DVE trailing-edge sweep (radians)
%   vecDVEROLL - valNELE x 1 vector of DVE roll angles (about x-axis) (radians)
%   vecDVEPITCH - valNELE x 1 vector of DVE pitch angles (about y-axis) (radians)
%   vecDVEYAW - valNELE x 1 vector of DVE yaw angles (about z-axis) (radians)
%   vecDVEAREA - valNELE x 1 vector of DVE area
%   vecDVENORM -  valNELE x 3 matrix of DVE normal vectors
%   matVLST - ? x 3 list of unique vertices, columns are (x,y,z) values
%   valNELE - total number of DVEs
%   matDVE - matrix of which DVE uses which vertices from the above list
%   matADJE - matADJE - ? x 3 adjacency matrix, where columns are: DVE | local edge | adjacent DVE
%   vecDVESYM - valNELE x 1 vector of which DVEs have symmetry on which edge (0 for no symmetry, 2 for local edge 2, 4 for local edge 4)
%   vecDVETIP - valNELE x 1 vector of which DVEs are at the wingtip. Similar format to vecDVESYM

% FUNCTIONS USED:
%   fcnPANELCORNERS
%   fcnPANEL2DVE
%   fcnGLOBSTAR

%% Preallocation
valNELE = sum(vecM.*vecN);
vecDVEPANEL   = nan(valNELE,1);
P1          = nan(valNELE,3);
P2          = nan(valNELE,3);
P3          = nan(valNELE,3);
P4          = nan(valNELE,3);
vecDVEWING  = nan(valNELE,1);
vecEnd      = cumsum(vecN.*vecM);

%% Assign Wing to Panel
panelEdges = reshape(permute(matGEOM,[1 3 2]),[],5);
[~,tempB,tempC] = unique(panelEdges,'rows','stable');
panelEdgesIdx = reshape(tempC,2,[])';
edge2wing = [(1:length(tempB))',nan(length(tempB),1)];
% Assign first edge to wing 1
edge2wing(1,2) = 1;

for n = 1:valPANELS
    curPanelEdges = panelEdgesIdx(n,:)';
    if max(edge2wing(curPanelEdges,2)) > 0
        wingIdx = max(edge2wing(curPanelEdges,2));
    else
        wingIdx = max(edge2wing(:,2))+1;
    end
    edge2wing(curPanelEdges,2)=wingIdx;
end
temp1 = reshape(edge2wing(tempC,2),2,[]);
panel2wing = temp1(1,:)';
clear tempB tempC temp1

%% Convert Panels to Corner Points to DVEs

points = [];
matTEPOINTS = [];
matLEPOINTS = [];

for i = 1:valPANELS
    temp_te = [];
    temp_le = [];
    
    rchord = matGEOM(1,4,i); repsilon = deg2rad(matGEOM(1,5,i));
    tchord = matGEOM(2,4,i); tepsilon = deg2rad(matGEOM(2,5,i));
    rLE = matGEOM(1,1:3,i);
    tLE = matGEOM(2,1:3,i);
    
    % Read panel corners
    % For DVE generation. Twist angle is handled along with dihedral angle
    panel4corners = reshape(fcnPANELCORNERS(rLE,tLE,rchord,tchord,repsilon,tepsilon),3,4)';
    
    if strcmpi(strATYPE{2}, 'PANEL') || size(strAIRFOIL,2) > 1
        % root
        airfoil(i,1).coord = dlmread(['airfoils/', strAIRFOIL{i,1}, '.dat'],'',1,0);
        airfoil(i,1).coord = airfoil(i,1).coord.*matGEOM(1,4,i);
        R = rotz(rad2deg(-repsilon));
        airfoil(i,1).coord = [R(1:2,1:2)*airfoil(i,1).coord']';
        [~,idx] = min(airfoil(i,1).coord(:,1));
        airfoil(i,1).le_idx = idx;
        ac_shift = -airfoil(i,1).coord(idx,:) + matGEOM(1,1:2:3,i);
        airfoil(i,1).coord = airfoil(i,1).coord + ac_shift;
        
        airfoil(i,2).coord = dlmread(['airfoils/', strAIRFOIL{i,2}, '.dat'],'',1,0);
        airfoil(i,2).coord = airfoil(i,2).coord.*matGEOM(2,4,i);
        R = rotz(rad2deg(-tepsilon));
        airfoil(i,2).coord = [R(1:2,1:2)*airfoil(i,2).coord']';
        [~,idx] = min(airfoil(i,2).coord(:,1));
        airfoil(i,2).le_idx = idx;
        ac_shift = -airfoil(i,2).coord(idx,:) + matGEOM(2,1:2:3,i);
        airfoil(i,2).coord = airfoil(i,2).coord + ac_shift;
        
        panel4corners(4,1:2:3) = (airfoil(i,1).coord(1,:) + airfoil(i,1).coord(end,:))./2;
        panel4corners(3,1:2:3) = (airfoil(i,2).coord(1,:) + airfoil(i,2).coord(end,:))./2;
    else
        airfoil = [];
    end
    
    % fcnPANEL2DVE takes four corners of a panel and outputs vertices of non-planer DVEs
    [~, LE_Left, LE_Right, TE_Left, TE_Right , tmpULS] = fcnPANEL2DVE(strATYPE, panel4corners, i, vecN, vecM, repsilon, tepsilon, rchord, tchord, airfoil, strSPACING);
    
    if strcmpi(strATYPE{2}, 'PANEL')
        vecM = vecM.*2;
        vecEnd = vecEnd.*2;
        valNELE = valNELE*2;
    end
    
    % Saving the TE and LE points so we can define those edges later on
    if vecPANELTE(i) == 1
        temp_te(:,:,1) = permute(TE_Left(end,:,:), [2 3 1]);
        temp_te(:,:,2) = permute(TE_Right(end,:,:), [2 3 1]);
    end
    
    if vecPANELLE(i) == 1
        temp_le(:,:,1) = permute(LE_Left(1,:,:), [2 3 1]);
        temp_le(:,:,2) = permute(LE_Right(1,:,:), [2 3 1]);
    end
    
    % WRITE RESULTS
    count = vecN(i)*vecM(i);
    idxStart = vecEnd(i)-count+1;
    idxEnd = vecEnd(i);
    
    vecDVEPANEL(idxStart:idxEnd,:) = repmat(i,count,1);
    
    % Write DVE WING Index
    vecDVEWING(idxStart:idxEnd,:) = repmat(panel2wing(i),count,1);
    
    % Write non-planer DVE coordinates
    P1(idxStart:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count,3);
    P2(idxStart:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count,3);
    P3(idxStart:idxEnd,:) = reshape(permute(TE_Right, [2 1 3]),count,3);
    P4(idxStart:idxEnd,:) = reshape(permute(TE_Left, [2 1 3]),count,3);
    vecULS(idxStart:idxEnd,:) = reshape(permute(tmpULS, [2 1 3]),count,1);
    
    % Creating "triangles" from the quadrilaterals on this panel. This is done in a way that SHOULD keep the normals "upwards" and not mixed
    vecULS = reshape(repmat(vecULS,1,2,1)', [],1,1);
    if strcmpi(strATYPE{2}, 'PANEL')
        temp_points = cat(3, permute(reshape([P1(idxStart:idxEnd,:) P3(idxStart:idxEnd,:) P2(idxStart:idxEnd,:)],[],3,3), [3 2 1]), ...
        permute(reshape([P1(idxStart:idxEnd,:) P4(idxStart:idxEnd,:) P3(idxStart:idxEnd,:)],[],3,3), [3 2 1]));
    else
        temp_points = cat(3, permute(reshape([P1(idxStart:idxEnd,:) P2(idxStart:idxEnd,:) P3(idxStart:idxEnd,:)],[],3,3), [3 2 1]), ...
        permute(reshape([P4(idxStart:idxEnd,:) P1(idxStart:idxEnd,:) P3(idxStart:idxEnd,:)],[],3,3), [3 2 1]));    
    end
    
    %     % Dealing with symmetry, assuming the symmetry plane is the YZ plane
    %     % This is done in a way that SHOULD keep the normals consistant
    %     if vecSYM(valPANELS) ~= 0
    %         len = length(P1(idxStart:idxEnd,1));
    %         temp_points_sym = cat(3, permute(reshape([P2(idxStart:idxEnd,:).*repmat([1 -1 1],len,1) P1(idxStart:idxEnd,:).*repmat([1 -1 1],len,1) P3(idxStart:idxEnd,:).*repmat([1 -1 1],len,1)],[],3,3), [3 2 1]), ...
    %             permute(reshape([P3(idxStart:idxEnd,:).*repmat([1 -1 1],len,1) P1(idxStart:idxEnd,:).*repmat([1 -1 1],len,1) P4(idxStart:idxEnd,:).*repmat([1 -1 1],len,1)],[],3,3), [3 2 1]));
    %         temp_points = cat(3, temp_points, temp_points_sym);
    %
    %         temp_shape = size(temp_le);
    %         temp_le = [temp_le; temp_le.*repmat([1 -1 1],temp_shape(1), 1, temp_shape(3))];
    %         temp_shape = size(temp_te);
    %         temp_te = [temp_te; temp_te.*repmat([1 -1 1],temp_shape(1), 1, temp_shape(3))];
    %     end
    
    points = cat(3, points, temp_points);
    
    matTEPOINTS = [matTEPOINTS; temp_te];
    matLEPOINTS = [matLEPOINTS; temp_le];
    
    
    clear LE_Left LE_Mid LE_Right TE_Right TE_Left ...
        imLEL imLER imTER imTEL ...
        idxStart idxEnd count temp_le temp_te temp_points temp_points_sym
end

verticeList = [reshape(permute(points,[2 1 3]),3,[],1)'];
[matVLST,idxVLST,matDVE] = unique(verticeList,'rows');
matDVE = reshape(matDVE,valNELE*2,3);

matPOINTS = permute(points, [3 2 1]);

end








