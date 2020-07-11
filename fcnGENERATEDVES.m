function [matPOINTS, matTEPOINTS, matLEPOINTS, vecDVESURFACE, vecDVEFLIP, vecDVEWING, vecDVEROTOR] = fcnGENERATEDVES(valPANELS, matGEOM, vecN, vecM, strAIRFOIL, vecPANELWING, vecPANELROTOR)

%% Preallocation
valNELE = sum(vecM.*vecN);
vecDVEPANEL   = nan(valNELE,1);
P1          = nan(valNELE,3);
P2          = nan(valNELE,3);
P3          = nan(valNELE,3);
P4          = nan(valNELE,3);
vecDVESURFACE  = [];
vecDVEWING  = [];
vecDVEROTOR  = [];
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
vecDVEFLIP = logical([]);

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
    
    if size(strAIRFOIL,2) > 1
        % root
        airfoil(i,1).coord = dlmread(['airfoils/', strAIRFOIL{i,1}, '.dat'],'',1,0);
        airfoil(i,1).coord = airfoil(i,1).coord.*(1./max(airfoil(i,1).coord(:,1)));
        
        [~,idx] = min(airfoil(i,1).coord(:,1));
        airfoil(i,1).le_idx = idx;
        
        airfoil(i,2).coord = dlmread(['airfoils/', strAIRFOIL{i,2}, '.dat'],'',1,0);
        airfoil(i,2).coord = airfoil(i,2).coord.*(1./max(airfoil(i,2).coord(:,1)));
        
        [~,idx] = min(airfoil(i,2).coord(:,1));
        airfoil(i,2).le_idx = idx;
    else
        airfoil = [];
    end
    
    % fcnPANEL2DVE takes four corners of a panel and outputs vertices of non-planer DVEs
    [~, LE_Left, LE_Right, TE_Left, TE_Right , tmpULS] = fcnPANEL2DVE(panel4corners, i, vecN, vecM, airfoil);
    
    % Saving the TE and LE points so we can define those edges later on
    temp_te(:,:,1) = permute(TE_Left(end,:,:), [2 3 1]);
    temp_te(:,:,2) = permute(TE_Right(end,:,:), [2 3 1]);
    
    temp_le(:,:,1) = permute(LE_Left(1,:,:), [2 3 1]);
    temp_le(:,:,2) = permute(LE_Right(1,:,:), [2 3 1]);
    
    
    % WRITE RESULTS
    count = vecN(i)*vecM(i);
    idxStart = vecEnd(i)-count+1;
    idxEnd = vecEnd(i);
    
    vecDVEPANEL(idxStart:idxEnd,:) = repmat(i,count,1);
    
    % Write DVE WING Index
    vecDVESURFACE = [vecDVESURFACE; repmat(panel2wing(i),count*2,1)];
    vecDVEWING = [vecDVEWING; repmat(vecPANELWING(i),count*2,1)];
    vecDVEROTOR = [vecDVEROTOR; repmat(vecPANELROTOR(i),count*2,1)];
    
    % Write non-planer DVE coordinates
    P1(idxStart:idxEnd,:) = reshape(permute(LE_Left, [2 1 3]),count,3);
    P2(idxStart:idxEnd,:) = reshape(permute(LE_Right, [2 1 3]),count,3);
    P3(idxStart:idxEnd,:) = reshape(permute(TE_Right, [2 1 3]),count,3);
    P4(idxStart:idxEnd,:) = reshape(permute(TE_Left, [2 1 3]),count,3);
    %     vecULS(idxStart:idxEnd,:) = reshape(permute(tmpULS, [2 1 3]),count,1);
    
    % Creating "triangles" from the quadrilaterals on this panel. This is done in a way that SHOULD keep the normals "upwards" and not mixed
    temp_points = cat(3, permute(reshape([P3(idxStart:idxEnd,:) P2(idxStart:idxEnd,:) P1(idxStart:idxEnd,:)],[],3,3), [3 2 1]), ...
        permute(reshape([P4(idxStart:idxEnd,:) P1(idxStart:idxEnd,:) P3(idxStart:idxEnd,:)],[],3,3), [3 2 1]));
    
    vecDVEFLIP = [vecDVEFLIP; true(length(idxStart:idxEnd),1); false(length(idxStart:idxEnd),1)];
    
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








