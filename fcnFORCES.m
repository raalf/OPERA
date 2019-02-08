function [CL, CDi, CY, e, vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN)
lim = 1e10;

%% Initializing
vecDVELIFT = nan(valNELE,1);
vecDVEDRAG = nan(valNELE,1);
vecDVESIDE = nan(valNELE,1);

% matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
% matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
% matDVESIDE_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);

matDVEDRAG_DIR = repmat(matUINF(1,:)./sqrt(sum(matUINF(1,:).^2,2)), valNELE, 1);
matDVELIFT_DIR = repmat(cross(matDVEDRAG_DIR(1,:), matSPANDIR(1,:), 2), valNELE, 1);
matDVESIDE_DIR = repmat(cross(matDVELIFT_DIR(1,:), matDVEDRAG_DIR(1,:), 2), valNELE,1);

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

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and U_INF.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA

% Velocity along the TE
locations = [0.1:0.4:0.9]';
uvw_te = locations.*reshape(matVUINF(matELST(vecTE,1),:)', 1, 3, []) + (1 - locations).*reshape(matVUINF(matELST(vecTE,2),:)', 1, 3, []);
fpl_te = locations.*reshape([xi_1 eta_1 xi_1.*0]', 1, 3, []) + (1 - locations).*reshape([xi_3 eta_3 xi_3.*0]', 1, 3, []);

tau = fpl_te(:,1,:);

for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(uvw_te(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);
end

A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
B_2 = matCOEFF(vecTEDVE,4); C_2 = matCOEFF(vecTEDVE,5); C_3 = matCOEFF(vecTEDVE,6);

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
F = [u_out v_out w_out];

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

%% Induced Forces
tmp51 = [];
tmp52 = [];
for jj = 1:length(vecTE)
    fpg_og = locations.*matVLST(matELST(vecTE(jj),1),:) + (1 - locations).*matVLST(matELST(vecTE(jj),2),:);
    
    fpg = fpg_og;
    len = size(fpg,1);
    
    wind(:,:,jj) = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    
    % Adding in bound induction call
    wind_b = zeros(len,3);

    for jjj = 1:length(vecWLE)
        dve2 = vecWLEDVE(jjj);

        p1 = fcnGLOBSTAR(matWVLST(matWELST(vecWLE(jjj),1),:) - matWCENTER(dve2,:), matWROTANG(dve2,:));
        p2 = fcnGLOBSTAR(matWVLST(matWELST(vecWLE(jjj),2),:) - matWCENTER(dve2,:), matWROTANG(dve2,:));
        
        midpoint = mean([p1; p2],1);
        vec = p2 - p1;
        hspan = abs(p1(1,1) - p2(1,1))./2;
        vec = vec./sqrt(sum(vec.^2,2));
        phi = acos(dot([1 0 0], vec, 2));
        if phi > pi/2
            phi = phi - pi;
        end

        xi = [linspace(p1(1,1), p2(1,1), 100)]';
        eta = [linspace(p1(1,2), p2(1,2), 100)]';       
        circ = 0.5.*matWCOEFF(dve2,1).*eta.^2 + matWCOEFF(dve2,2).*eta + 0.5.*matWCOEFF(dve2,3).*xi.^2 + matWCOEFF(dve2,4).*xi + matWCOEFF(dve2,5).*xi.*eta + matWCOEFF(dve2,6);
%         circ = 0.5.*matWCOEFF(dve2,1).*eta.^2 + matWCOEFF(dve2,2).*eta + matWCOEFF(dve2,5).*xi.*eta + matWCOEFF(dve2,6);

        if p2(1,1) > p1(1,1)
            tmp6 = [linspace(-hspan, hspan, 100)]';
        else
            tmp6 = [linspace(hspan, -hspan, 100)]';
        end
        fp_0 = fcnGLOBSTAR(fpg - matWCENTER(dve2,:), matWROTANG(dve2,:)) - midpoint;
        [aloc, bloc, cloc] = fcnBOUNDIND(repmat(hspan,len,1), repmat(phi,len,1), [-fp_0(:,2) fp_0(:,1) fp_0(:,3)]);
%         aloc = [aloc(:,2), -aloc(:,1), aloc(:,3)];
%         bloc = [bloc(:,2), -bloc(:,1), bloc(:,3)];
%         cloc = [cloc(:,2), -cloc(:,1), cloc(:,3)];
        
        D = [aloc bloc cloc]; % C is higher order
        D = reshape(reshape(D', 1, 9, []), 3, 3, len);
               
        coeff = polyfit(tmp6, circ, 2);
        coeff = fliplr(coeff);
        
%         hFig23 = figure(23);
%         clf(23);
% %         tmp6 = sqrt(sum(tmp6.^2,2));
%         plot(tmp6, circ, '-k');
% %       scatter3(xi, eta, circ, 'ok')
%         hold on
%         plot(tmp6, tmp6.^2.*coeff(3) + tmp6.*coeff(2) + coeff(1), '-.or','MarkerSize',2)
%         hold off
%         grid minor
%         box on
%         axis tight
% %         view([0 0]);
        
        w_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
        w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);
                
        % To global
        wind_b = wind_b + fcnSTARGLOB(w_ind, repmat(matWROTANG(dve2,:),len,1));
    
    end
    
    wind(:,:,jj) = wind(:,:,jj) + wind_b;
    
    tmp51 = [tmp51; fpg_og];
    tmp52 = [tmp52; wind(:,:,jj)];   
end

%     hold on
%     quiver3(tmp51(:,1), tmp51(:,2), tmp51(:,3), tmp52(:,1), tmp52(:,2), tmp52(:,3))
%     hold off


% Velocity along the TE
for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(wind(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);
    
%     hFig12 = figure(12);
%     clf(12);
%     x = linspace(0, dist(i), 20);
%     subplot(3,1,1);
%     plot(x, u(i,1).*x.^2 + u(i,2).*x + u(i,3), '-ok');
%     hold on
%     scatter(tau(:,1,i), uvw(:,1),'d')
%     hold off
%     
%     subplot(3,1,2);
%     plot(x, v(i,1).*x.^2 + v(i,2).*x + v(i,3), '-ok');
%     hold on
%     scatter(tau(:,1,i), uvw(:,2),'d')
%     hold off
%     
%     subplot(3,1,3);
%     plot(x, w(i,1).*x.^2 + w(i,2).*x + w(i,3), '-ok');
%     hold on
%     scatter(tau(:,1,i), uvw(:,3),'d')
%     hold off
    
end

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
F = [u_out v_out w_out];

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% Splitting induced force into lift and side forces
liftind(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sideind(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);
dragind(vecTEDVE,1) = dot(F, matDVEDRAG_DIR(vecTEDVE,:), 2);

% [~,b] = sort(tmp51(:,2),'ascend');
% tmp51 = tmp51(b,:);
% tmp52 = tmp52(b,:);
% % load('Stuff/downwash_goland_wing_vap3.mat');
% hFig34 = figure(34);
% clf(34)
% plot(tmp51(:,2), tmp52(:,3),'-^m')
% hold on
% % scatter(XDATA, YDATA, 'ok');
% hold off
% box on
% grid minor

vecDVELIFT = liftfree + liftind;
vecDVESIDE = sidefree + sideind;
vecDVEDRAG = dragind;

% vecDVELIFT = liftfree;
% vecDVESIDE = sidefree;
% vecDVEDRAG = dragind;

%% Output
CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*valAREA);
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*valAREA);
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*valAREA);
e = (CL.^2)./(pi.*((valSPAN.^2)./valAREA).*CDi);
fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

end