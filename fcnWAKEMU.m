function [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE)

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
if strcmpi(strATYPE{2},'PANEL')
    pts_loc(:,:,1) = fcnGLOBSTAR(pts(:,:,1) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    pts_loc(:,:,2) = fcnGLOBSTAR(pts(:,:,2) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    pts_loc(:,:,3) = fcnGLOBSTAR(pts(:,:,3) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    res1_2 = [   sum([0.5.*pts_loc(:,2,1).^2 pts_loc(:,2,1) 0.5.*pts_loc(:,1,1).^2 pts_loc(:,1,1) pts_loc(:,2,1).*pts_loc(:,1,1) ones(size(pts_loc(:,1,1)))].*matCOEFF(vecTEDVE(:,2),:),2); ...
        sum([0.5.*pts_loc(:,2,2).^2 pts_loc(:,2,2) 0.5.*pts_loc(:,1,2).^2 pts_loc(:,1,2) pts_loc(:,2,2).*pts_loc(:,1,2) ones(size(pts_loc(:,1,2)))].*matCOEFF(vecTEDVE(:,2),:),2); ...
        sum([0.5.*pts_loc(:,2,3).^2 pts_loc(:,2,3) 0.5.*pts_loc(:,1,3).^2 pts_loc(:,1,3) pts_loc(:,2,3).*pts_loc(:,1,3) ones(size(pts_loc(:,1,3)))].*matCOEFF(vecTEDVE(:,2),:),2)];
    
    res1 = res1_2 + res1_1;
else
    res1 = res1_1;
end
res1 = reshape(res1,[],3);

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
   vecWEMU(matWE2GRID(1,:)) = (vecWVMU(matWVGRID(1,:)) + vecWVMU(matWVGRID(2,:)))./2;
   vecWEMU(matWEGRID(2,:)) = (vecWEMU(matWEGRID(1,:)) + vecWEMU(matWEGRID(3,:)))./2;
end

end

