function [CL, CDi, CY, e, vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF)

%% Initializing
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

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and U_INF.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA

% Velocity along the TE
locations = [0.01:0.4:0.99]';
uvw_te = locations.*reshape(matVUINF(matELST(vecTE,1),:)', 1, 3, []) + (1 - locations).*reshape(matVUINF(matELST(vecTE,2),:)', 1, 3, []);
fpl_te = locations.*reshape([xi_1 eta_1 xi_1.*0]', 1, 3, []) + (1 - locations).*reshape([xi_3 eta_3 xi_3.*0]', 1, 3, []);
tau = fpl_te - reshape([xi_1 eta_1 xi_1.*0]', 1, 3, []);
tau = sqrt(sum(tau.^2, 2));

for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(uvw_te(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);
end

A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
B_2 = matCOEFF(vecTEDVE,4); C_2 = matCOEFF(vecTEDVE,5); C_3 = matCOEFF(vecTEDVE,6);

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
F = [u_out v_out w_out];

v_te = [xi_3 - xi_1, eta_3 - eta_1, eta_1.*0];
d_te = sqrt(sum(v_te.^2,2));

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% uvw = fcnGLOBSTAR(matUINF(vecTEDVE,:), matROTANG(vecTEDVE,:));
% % For readability:
% A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
% B_2 = matCOEFF(vecTEDVE,4); C_2 = matCOEFF(vecTEDVE,5); C_3 = matCOEFF(vecTEDVE,6);
% u = uvw(:,1); v = uvw(:,2); w = uvw(:,3); rho = valDENSITY;
% 
% % Analytically integrating circulation along the TE vector
% % Projected circulation onto trailing edge
% % Spanwise circulation
% tmp21 = -(((A_1.*E.^2+2.*C_2.*E+B_1).*xi_1.^2+((A_1.*E.^2+2.*C_2.*E+B_1).*xi_3+(3.*A_1.*E+3.*C_2).*D_TE+3.*A_2.*E+3.*B_2).*xi_1+(A_1.*E.^2+2.*C_2.*E+B_1).*xi_3.^2+((3.*A_1.*E+3.*C_2).*D_TE+3.*A_2.*E+3.*B_2).*xi_3+3.*A_1.*D_TE.^2+6.*A_2.*D_TE+6.*C_3).*w.*(xi_1-xi_3).*(eta_1-eta_3).*rho)./0.6e1;
% tmp22 = (((A_1.*E.^2+2.*C_2.*E+B_1).*xi_1.^2+((A_1.*E.^2+2.*C_2.*E+B_1).*xi_3+(3.*A_1.*D_TE+3.*A_2).*E+3.*C_2.*D_TE+3.*B_2).*xi_1+(A_1.*E.^2+2.*C_2.*E+B_1).*xi_3.^2+((3.*A_1.*D_TE+3.*A_2).*E+3.*C_2.*D_TE+3.*B_2).*xi_3+3.*A_1.*D_TE.^2+6.*A_2.*D_TE+6.*C_3).*w.*(xi_1-xi_3).^2.*rho)./0.6e1;
% tmp33 = -(((A_1.*E.^2+2.*C_2.*E+B_1).*xi_1.^2+((A_1.*E.^2+2.*C_2.*E+B_1).*xi_3+(3.*A_1.*D_TE+3.*A_2).*E+3.*C_2.*D_TE+3.*B_2).*xi_1+(A_1.*E.^2+2.*C_2.*E+B_1).*xi_3.^2+((3.*A_1.*D_TE+3.*A_2).*E+3.*C_2.*D_TE+3.*B_2).*xi_3+3.*A_1.*D_TE.^2+6.*A_2.*D_TE+6.*C_3).*(v.*xi_1-v.*xi_3-u.*(eta_1-eta_3)).*rho.*(xi_1-xi_3))./0.6e1;
% 
% % tmp21 = ((B_1./0.2e1+C_2.*E).*xi_1.^2+((B_1./0.2e1+C_2.*E).*xi_3+0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_1+(B_1./0.2e1+C_2.*E).*xi_3.^2+(0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_3+(3.*C_3)).*rho.*(eta_1-eta_3).*w.*(xi_3-xi_1)./0.3e1;
% % tmp22 = ((B_1./0.2e1+C_2.*E).*xi_1.^2+((B_1./0.2e1+C_2.*E).*xi_3+0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_1+(B_1./0.2e1+C_2.*E).*xi_3.^2+(0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_3+(3.*C_3)).*rho.*w.*(xi_3-xi_1).^2./0.3e1; 
% % tmp33 = -(v.*xi_1-v.*xi_3-u.*(eta_1-eta_3)).*rho.*((B_1./0.2e1+C_2.*E).*xi_1.^2+((B_1./0.2e1+C_2.*E).*xi_3+0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_1+(B_1./0.2e1+C_2.*E).*xi_3.^2+(0.3e1./0.2e1.*C_2.*D_TE+0.3e1./0.2e1.*B_2).*xi_3+(3.*C_3)).*(xi_1-xi_3)./0.3e1;
% 
% v_te = [xi_3 - xi_1, eta_3 - eta_1, eta_1.*0];
% d_te = sqrt(sum(v_te.^2,2));
% 
% % This is the special force
% F = fcnSTARGLOB(-[tmp21 tmp22 tmp33]./d_te, matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

%% Induced Forces
tmp51 = [];
tmp52 = [];
for jj = 1:length(vecTE)
%     locations = [0.1:0.4:0.9]';
    fpg_og = locations.*matVLST(matELST(vecTE(jj),1),:) + (1 - locations).*matVLST(matELST(vecTE(jj),2),:);
    
    fpg = fpg_og;
    len = size(fpg,1);
    
    wind(:,:,jj) = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    
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

        xi_mp = midpoint(1);
        
        xi = linspace(p1(1,1), p2(1,1), 100);
        eta = linspace(p1(1,2), p2(1,2), 100);
        circ = 0.5.*matCOEFF(dve2,3).*xi.^2 + matCOEFF(dve2,4).*xi + matCOEFF(dve2,5).*xi.*eta + matCOEFF(dve2,6);
        
        coeff = polyfit(xi - xi_mp, circ, 2);
        coeff = fliplr(coeff);
        
        w_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
        w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);
        
        wind_b = wind_b + fcnSTARGLOB(w_ind, repmat(matROTANG(dve2,:),len,1));
    
    end
    wind(:,:,jj) = wind(:,:,jj) - wind_b;
    
    tmp51 = [tmp51; fpg_og];
    tmp52 = [tmp52; wind(:,:,jj)];
    
%         wcoeff = polyfit(fpg(:,2), wind(:,3), 2);
%         hFig30 = figure(30);
%         clf(30);
%         scatter(fpg_og(:,2), wind(:,3),'mo')
%         hold on
%         x = linspace(fpg_og(1,2), fpg_og(end,2),100);
%         plot(x, wcoeff(1).*x.^2 + wcoeff(2).*x + wcoeff(3), '--^k')
%         hold off
%         box on
%         grid minor

%     hspan = fcnGLOBSTAR(matVLST(matELST(vecTE(jj),1),:) - matVLST(matELST(vecTE(jj),2),:), matROTANG(vecTEDVE(jj),:));
%     eta8 = hspan.*0.8;
%     hspan = abs(hspan(1))/2;
%     s = (matVLST(matELST(vecTE(jj),1),:) - matVLST(matELST(vecTE(jj),2),:))./(hspan.*2);
%     points = fcnGLOBSTAR(fpg_og - matCENTER(vecTEDVE(jj),:), repmat(matROTANG(vecTEDVE(jj),:), len, 1));
%     te_circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) points(:,1).*points(:,2) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2);
% %     te_circ = sum([0.5.*points(:,1).^2 points(:,1) points(:,1).*points(:,2) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),3:end),2);
%     tempr = cross(valDENSITY.*wind, repmat(s,len,1), 2).*te_circ;
%     R = (tempr(1,:) + 4.*tempr(2,:) + tempr(3,:)).*(eta8)/3;
%     R = R + 2.*(7.*tempr(1,:) - 8.*tempr(2,:) + 7.*tempr(3,:)).*(hspan-eta8)./3;        
%     
%     dragind(vecTEDVE(jj),1) = dot(R, matDVEDRAG_DIR(vecTEDVE(jj),:), 2);
%     liftind(vecTEDVE(jj),1) = dot(R, matDVELIFT_DIR(vecTEDVE(jj),:), 2);
%     sideind(vecTEDVE(jj),1) = dot(R, matDVESIDE_DIR(vecTEDVE(jj),:), 2);
end

% Velocity along the TE
for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(wind(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);
end

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY);
F = [u_out v_out w_out];

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% Splitting induced force into lift and side forces
liftind(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sideind(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);
dragind(vecTEDVE,1) = dot(F, matDVEDRAG_DIR(vecTEDVE,:), 2);

%%
% [~,b] = sort(tmp51(:,2),'ascend');
% tmp51 = tmp51(b,:);
% tmp52 = tmp52(b,:);
% load('Stuff/downwash_goland_wing_vap3.mat');
% hFig34 = figure(34);
% clf(34)
% plot(tmp51(:,2), tmp52(:,3),'-^m')
% hold on
% scatter(XDATA, YDATA, 'ok');
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
CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*sum(vecDVEAREA));
e = (CL.^2)./(pi.*((max(matVLST(:,2))-min(matVLST(:,2))).^2)./sum(vecDVEAREA).*CDi);
fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

end