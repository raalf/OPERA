function [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE, valAZPERREV, valTIMESTEP)

%% Spanwise circulation in the front vertices of newest wake row
pts(:,:,1) = matWVLST(matWELST(vecWLE,1),:);
pts(:,:,2) = matWVLST(matWELST(vecWLE,2),:);
pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;

pts_loc(:,:,1) = fcnGLOBSTAR(pts(:,:,1) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
pts_loc(:,:,2) = fcnGLOBSTAR(pts(:,:,2) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
pts_loc(:,:,3) = fcnGLOBSTAR(pts(:,:,3) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
res1_1 = [   sum([0.5.*pts_loc(:,2,1).^2 pts_loc(:,2,1) 0.5.*pts_loc(:,1,1).^2 pts_loc(:,1,1) pts_loc(:,2,1).*pts_loc(:,1,1) ones(size(pts_loc(:,1,1)))].*matCOEFF(vecTEDVE(:,1),:),2); ...
    sum([0.5.*pts_loc(:,2,2).^2 pts_loc(:,2,2) 0.5.*pts_loc(:,1,2).^2 pts_loc(:,1,2) pts_loc(:,2,2).*pts_loc(:,1,2) ones(size(pts_loc(:,1,2)))].*matCOEFF(vecTEDVE(:,1),:),2); ...
    sum([0.5.*pts_loc(:,2,3).^2 pts_loc(:,2,3) 0.5.*pts_loc(:,1,3).^2 pts_loc(:,1,3) pts_loc(:,2,3).*pts_loc(:,1,3) ones(size(pts_loc(:,1,3)))].*matCOEFF(vecTEDVE(:,1),:),2)];

res1 = reshape(res1_1,[],3);

%% Assigning circulation values to leading edge of wake row (vertices and edge midpoints)
vecWVMU(matWELST(vecWLE,1)) = res1(:,1);
vecWVMU(matWELST(vecWLE,2)) = res1(:,2);
vecWEMU(vecWLE) = res1(:,3);

% vecWVMU(matWELST(vecWOTE,1)) = 0;
% vecWVMU(matWELST(vecWOTE,2)) = 0;
% vecWEMU(vecWOTE) = 0;

if (size(matWVGRID,1)) <= 2 || strcmpi(strATYPE{3},'STEADY')
    vecWVMU(matWVGRID(2:end,:)) = repmat(vecWVMU(matWVGRID(1,:))', size(matWVGRID,1) - 1, 1);
    vecWEMU(matWEGRID(2:end,:)) = repmat(vecWEMU(matWEGRID(1,:))', size(matWEGRID,1) - 1, 1);
    vecWEMU(matWE2GRID) = repmat(vecWVMU(matWVGRID(1,:))', size(matWE2GRID,1), 1);
else
    if valTIMESTEP > valAZPERREV
        % Applying wake vertex strengths over the last revolution to all
        % revolutions
        if mod(valTIMESTEP, valAZPERREV) == 0
            tmp = vecWVMU(matWVGRID(1,:))'; % Treated different if only 1 row (cause it comes out as a row vector, not a column vector)
        else
            tmp = vecWVMU(matWVGRID(1:mod(valTIMESTEP, valAZPERREV) + 1,:));
        end
        % Updating wake strengths by taking the last revolution, applying
        % to all revolution, and adding remainders all the way to the
        % oldest row
        vecWVMU(matWVGRID) = [repmat(vecWVMU(matWVGRID(1:valAZPERREV,:)), floor(valTIMESTEP/valAZPERREV), 1); tmp];
        
        % Updating edges
        vecWEMU(matWE2GRID) = (vecWVMU(matWVGRID(1:end-1,:)) + vecWVMU(matWVGRID(2:end,:)))./2;
        % Adding in 2nd edge from LE
        vecWEMU(matWEGRID(2,:)) = (vecWEMU(matWEGRID(1,:)) + vecWEMU(matWEGRID(3,:)))./2;

        % Updating edges all the way back to the oldest row of wake
        % elements
        if mod(valTIMESTEP, valAZPERREV) == 0
            tmp2 = vecWEMU(matWEGRID(1,:))';
        else
            tmp2 = vecWEMU(matWEGRID(1:mod(valTIMESTEP, valAZPERREV)*2 + 1,:));
        end
            vecWEMU(matWEGRID) = [repmat(vecWEMU(matWEGRID(1:(valAZPERREV*2),:)), floor(valTIMESTEP/valAZPERREV), 1); tmp2];
        
    else
        vecWEMU(matWE2GRID(1,:)) = (vecWVMU(matWVGRID(1,:)) + vecWVMU(matWVGRID(2,:)))./2;
        vecWEMU(matWEGRID(2,:)) = (vecWEMU(matWEGRID(1,:)) + vecWEMU(matWEGRID(3,:)))./2;
    end
end

end

