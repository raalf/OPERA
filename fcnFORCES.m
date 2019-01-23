function [CL, CDi, e, vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA)

% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 20);

% Forces
% Jank is in session, The Right Honourable FORLOOP is preciding
% VECTORIZE IT TRAVIS
vecDVELIFT = nan(valNELE,1);
vecDVEDRAG = nan(valNELE,1);
matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
matDVESIDEF_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);

liftfree = nan(valNELE,1);
liftind = nan(valNELE,1);
sidefree = nan(valNELE,1);
sideind = nan(valNELE,1);
dragind = nan(valNELE,1);

tmp1 = [];
tmp2 = [];
tmp3 = [];
for jj = 1:size(vecTE,1)
    % Lift
    pts = matVLST(matELST(vecTE(jj),:),:); % End points of TE edge segement
    pts_loc = fcnGLOBSTAR(pts - matCENTER(vecTEDVE(jj),:), repmat(matROTANG(vecTEDVE(jj),:),2,1)); % In element local (The TE DVE where this edge is from)
    
    len = 1000; % Number of divisions of this line (Jank)
    points = [linspace(pts_loc(1,1), pts_loc(2,1), len)' linspace(pts_loc(1,2), pts_loc(2,2), len)' linspace(pts_loc(1,3), pts_loc(2,3), len)'];
    distance = sqrt(sum((pts_loc(2,:) - pts_loc(1,:)).^2,2)); % Length of entire TE edge
    vec = (pts_loc(1,:) - pts_loc(2,:))./distance; % Direction of this edge (in local)
    
    % Circulation at the points along the edge (oriented along the edge)
    circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2).*vec;
    circ = fcnSTARGLOB(circ, repmat(matROTANG(vecTEDVE(jj),:), len, 1)); % Translate to global
    F = cross( repmat(valDENSITY.*matUINF(vecTEDVE(jj),:), len, 1), circ, 2); % A special guest mix on the track, Kutty J
    spans = repmat(distance/len, len, 1);
    spans([1 len],:) = spans([1 len],:).*0.5;
    liftfree(vecTEDVE(jj),1) = sqrt(sum((sum(F.*spans,1)).^2,2)); % Multiply by the length of the discretization, and sum
    vecDVELIFT(vecTEDVE(jj),1) = liftfree(vecTEDVE(jj),1);
    
    % Drag
    %         fpg = fcnSTARGLOB(points, repmat(matWROTANG(vecWLEDVE(jj),:),len,1)) + matWCENTER(vecWLEDVE(jj),:);
    %         fpg = fpg - repmat(matWDVECT(vecWLEDVE(jj),:,2),len,1).*1e-1;
    %         w_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    %
    %         vec_glob = (pts(1,:) - pts(2,:));
    %         vec_glob = vec_glob./sqrt(sum(vec_glob.^2,2));
    %
    %         tmp5 = fcnSTARGLOB(points, repmat(matWROTANG(vecWLEDVE(jj),:),len,1)) + matWCENTER(vecWLEDVE(jj),:);
    %         s1 = tmp5 + vec_glob.*(spans/2);
    %         s2 = tmp5 - vec_glob.*(spans/2);
    %
    %         circ_tmp = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2);
    %
    %         r1 = s1-fpg;
    %         r2 = s2-fpg;
    %         r0 = r1 - r2;
    %
    %         d1 = cross(r1,r2,2);
    %         d = sqrt(sum(d1.^2,2))./sqrt(sum(r0.^2,2));
    %         cosb1 = (dot(r0,r1,2)./(sqrt(sum(r0.^2,2)).*sqrt(sum(r0.^2,2))));
    %         cosb2 = (dot(r0,r2,2)./(sqrt(sum(r0.^2,2)).*sqrt(sum(r0.^2,2))));
    %
    %         q_theta = (circ_tmp./(4.*pi.*d)).*(cosb1-cosb2);
    %         vdirection = cross(r1,r2,2);
    %         vdirection = vdirection./sqrt(sum(vdirection.^2,2));
    %
    %         q_12 = q_theta.*vdirection;
    %
    %         w_ind = w_ind + q_12;
    %
    %
    %         w_ind(1,:) = w_ind(1,:).*0;
    %         w_ind(end,:) = w_ind(end,:).*0;
    %
    %         tmp1 = [tmp1; fpg];
    %         tmp2 = [tmp2; w_ind];
    %
    %         F_ind = cross( valDENSITY.*w_ind, circ, 2);
    %         spans(2,:) = spans(2,:) + 0.5.*(distance/len);
    %         spans(end-1,:) = spans(end-1,:) + 0.5.*(distance/len);
    %         force = sum(F_ind.*spans,1);
    %
    %         tmp3 = [tmp3; F_ind];
    %
    %         dragind(vecTEDVE(jj),1) = dot(force, matDVEDRAG_DIR(vecTEDVE(jj),:), 2);
    %         liftind(vecTEDVE(jj),1) = dot(force, matDVELIFT_DIR(vecTEDVE(jj),:), 2);
    %
    %         vecDVEDRAG(vecTEDVE(jj),1) = dragind(vecTEDVE(jj),1);
    %         vecDVELIFT(vecTEDVE(jj),1) = liftfree(vecTEDVE(jj),1) + liftind(vecTEDVE(jj),1);
    
end

CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*sum(vecDVEAREA));
e = (CL.^2)./(pi.*((max(matVLST(:,2))-min(matVLST(:,2))).^2)./sum(vecDVEAREA).*CDi);
fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

%     hold on
%     quiver3(tmp1(:,1), tmp1(:,2), tmp1(:,3), tmp2(:,1), tmp2(:,2), tmp2(:,3),'b')
%     quiver3(tmp1(:,1), tmp1(:,2), tmp1(:,3), tmp3(:,1), tmp3(:,2), tmp3(:,3),'r')
%     setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;
%     hold off
%     view([36 15]);

end