function [CL, CDi, CY, e, vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG)

vecDVELIFT = nan(valNELE,1);
vecDVEDRAG = nan(valNELE,1);
vecDVESIDE = nan(valNELE,1);
matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
matDVESIDE_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);

liftfree = nan(valNELE,1);
sidefree = nan(valNELE,1);
liftind = nan(valNELE,1);
sideind = nan(valNELE,1);
dragind = nan(valNELE,1);

%% To element local
xi_1 = permute(matPLEX(1,1,vecTEDVE),[3 2 1]);
xi_2 = permute(matPLEX(2,1,vecTEDVE),[3 2 1]);
xi_3 = permute(matPLEX(3,1,vecTEDVE),[3 2 1]);

eta_1 = permute(matPLEX(1,2,vecTEDVE),[3 2 1]);
eta_2 = permute(matPLEX(2,2,vecTEDVE),[3 2 1]);
eta_3 = permute(matPLEX(3,2,vecTEDVE),[3 2 1]);

idx_flp = abs(xi_2 - xi_3) < 1e-3 & abs(xi_1 - xi_2) > 1e-3;
xi_tmp(idx_flp) = xi_3(idx_flp);
xi_3(idx_flp) = xi_1(idx_flp);
xi_1(idx_flp) = xi_tmp(idx_flp);
eta_tmp(idx_flp) = eta_3(idx_flp);
eta_3(idx_flp) = eta_1(idx_flp);
eta_1(idx_flp) = eta_tmp(idx_flp);

idx_rrg = eta_2 < eta_1;
eta_tmp(idx_rrg) = eta_2(idx_rrg);
eta_2(idx_rrg) = eta_1(idx_rrg);
eta_1(idx_rrg) = eta_tmp(idx_rrg);

% Checking which elements are on the element
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and the
% density times U_INF. Circulation equation is put solely in terms of xi
% by defining eta in terms of xi. It is integrated from xi_1 to xi_3.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA
uvw = fcnGLOBSTAR(matUINF(vecTEDVE,:), matROTANG(vecTEDVE,:));

% For readability:
A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
B_2 = matCOEFF(vecTEDVE,4); C_3 = matCOEFF(vecTEDVE,5);
u = uvw(:,1); v = uvw(:,2); w = uvw(:,3);

% Analytically integrating circulation along the TE vector
% Projected circulation onto trailing edge
% term = ((A_1.*E.^2 + B_1).*xi_1.^2 + ((A_1.*E.^2 + B_1).*xi_3 + (3.*A_1.*D_TE + 3.*A_2).*E + 3.*B_2).*xi_1 + (A_1.*E.^2 + B_1).*xi_3.^2 + ((3.*A_1.*D_TE + 3.*A_2).*E + 3.*B_2).*xi_3 + 3.*A_1.*D_TE.^2 + 6.*A_2.*D_TE + 6.*C_3);
% Spanwise circulation
term = (B_1.*xi_1.^2 + (B_1.*xi_3 + 3.*B_2).*xi_1 + B_1.*xi_3.^2 + 3.*B_2.*xi_3 + 6.*C_3);
tmp21 = -(1/6).*((xi_1 - xi_3).*valDENSITY.*term.*w.*(eta_1 - eta_3));
tmp22 =  (1/6).*((xi_1 - xi_3).^2.*valDENSITY.*term.*w.*(eta_1 - eta_3));
tmp33 = -(1/6).*(xi_1.*v - xi_3.*v - u.*(eta_3 - eta_1)).*term.*(xi_1 - xi_3).*valDENSITY;

% This is the special force
F = fcnSTARGLOB(-[tmp21 tmp22 tmp33]./(xi_3 - xi_1), matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

vecDVELIFT(vecTEDVE,1) = liftfree(vecTEDVE,1);
vecDVESIDE(vecTEDVE,1) = sidefree(vecTEDVE,1);

%% Induced Forces
te_vec = matVLST(matELST(vecTE,1),:) - matVLST(matELST(vecTE,2),:);

for jj = 1:length(vecTE)
    
    dve = vecTEDVE(jj);
    
    locations = [0.05:0.1:0.95]';
    
    fpg_og = locations.*matVLST(matELST(vecTE(jj),1),:) + (1 - locations).*matVLST(matELST(vecTE(jj),2),:);
    
%     fpg_og = [0.50.*matVLST(matELST(vecTE(jj),1),:) + 0.50.*matVLST(matELST(vecTE(jj),2),:); ...
%         0.75.*matVLST(matELST(vecTE(jj),1),:) + 0.25.*matVLST(matELST(vecTE(jj),2),:);...
%         0.25.*matVLST(matELST(vecTE(jj),1),:) + 0.75.*matVLST(matELST(vecTE(jj),2),:);...
%         0.15.*matVLST(matELST(vecTE(jj),1),:) + 0.85.*matVLST(matELST(vecTE(jj),2),:);...
%         0.85.*matVLST(matELST(vecTE(jj),1),:) + 0.15.*matVLST(matELST(vecTE(jj),2),:)
%         0.15.*matVLST(matELST(vecTE(jj),1),:) + 0.85.*matVLST(matELST(vecTE(jj),2),:);...
%         0.85.*matVLST(matELST(vecTE(jj),1),:) + 0.15.*matVLST(matELST(vecTE(jj),2),:)];
    
    fpg = fpg_og;
    len = size(fpg,1);
%     fpg = fpg_og + [ 0 0 2e-1];
%     fpg = fpg_og + [ 0 0 -2e-2];
%     fpg = fpg_og + [2e-2 0 0 ];
%     fpg = fpg_og + [-2e-1 0 0 ];
    
    wind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    
    % Adding in bound induction call
    wind_b = zeros(len,3);
    for jjj = 1:length(vecTE)
        dve2 = vecTEDVE(jjj);


        p1 = fcnGLOBSTAR(matVLST(matELST(vecTE(jjj),1),:) - matCENTER(dve2,:), matROTANG(dve2,:));
        p2 = fcnGLOBSTAR(matVLST(matELST(vecTE(jjj),2),:) - matCENTER(dve2,:), matROTANG(dve2,:));
        
        hspan = abs(p1(1) - p2(1))./2;
        midpoint = mean([p1; p2],1);
        vec = p2 - p1;
        vec = vec./sqrt(sum(vec.^2,2));
        phi = acos(dot([1 0 0], vec, 2));
        if phi > pi/2
            phi = phi - pi;
        end
        fp_0 = fcnGLOBSTAR(fpg - matCENTER(dve2,:), matROTANG(dve2,:)) - midpoint;
        
        [aloc, bloc, cloc] = fcnBOUNDIND(repmat(hspan,len,1), repmat(-phi,len,1), [fp_0(:,2) fp_0(:,1) fp_0(:,3)]);
        aloc = [aloc(:,2) aloc(:,1) aloc(:,3)];
        bloc = [bloc(:,2) bloc(:,1) bloc(:,3)];
        cloc = [cloc(:,2) cloc(:,1) cloc(:,3)];
        
        D = [aloc bloc cloc]; % C is higher order
        D = reshape(reshape(D', 1, 9, []), 3, 3, len);
        
        % We need to translate the polynomial coefficients from xi = 0 to
        % xi = xi_midpoint
        xi_mp = midpoint(1);
        
        xi = linspace(p1(1,1), p2(1,1), 100);
        circ = 0.5.*matCOEFF(dve2,3).*xi.^2 + matCOEFF(dve2,4).*xi + matCOEFF(dve2,5);

        coeff = polyfit(xi - xi_mp, circ, 2);
        coeff = fliplr(coeff);
                
        w_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
        w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);
        
        wind_b = wind_b + fcnSTARGLOB(w_ind, repmat(matROTANG(dve2,:),len,1));
    end
    

%     hold on
%     quiver3(fpg(:,1), fpg(:,2), fpg(:,3), wind(:,1), wind(:,2), wind(:,3),'b');
%     quiver3(fpg(:,1), fpg(:,2), fpg(:,3), wind_b(:,1), wind_b(:,2), wind_b(:,3),'r')
%     hold off
    
    wind = wind - wind_b;
    points = fcnGLOBSTAR(fpg_og - matCENTER(vecTEDVE(jj),:), repmat(matROTANG(vecTEDVE(jj),:), len, 1));
    vec = points(1,:) - points(end,:);
%     span = abs(vec(:,1));
    span = sqrt(sum(vec.^2,2));
    vec = vec./sqrt(sum(vec.^2,2));
%     te_circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2).*vec;
    te_circ = sum([0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),3:end),2).*vec;
    te_circ = fcnSTARGLOB(te_circ, repmat(matROTANG(vecTEDVE(jj),:), len,1));
    F_ind = cross(valDENSITY.*wind, te_circ, 2).*(span/len);
    hold on
    quiver3(fpg(:,1), fpg(:,2), fpg(:,3), F_ind(:,1), F_ind(:,2), F_ind(:,3),'m')
    quiver3(fpg(:,1), fpg(:,2), fpg(:,3), wind(:,1), wind(:,2), wind(:,3),'b');
    hold off
    F_ind = sum(F_ind,1);
    dragind(vecTEDVE(jj),1) = dot(F_ind, matDVEDRAG_DIR(vecTEDVE(jj),:), 2);
    
end

vecDVEDRAG = dragind;

%% Output
CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*sum(vecDVEAREA));
e = (CL.^2)./(pi.*((max(matVLST(:,2))-min(matVLST(:,2))).^2)./sum(vecDVEAREA).*CDi);
fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

end