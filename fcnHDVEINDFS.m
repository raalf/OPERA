function [infl_glob] = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matCENTER)

fpl = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

%%
xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

% Checking which elements are on the element, moving them off by a small
% amount
margin = 1e-5;
le_eta = eta_2 + (fpl(:,1) - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
te_eta = eta_1 + (fpl(:,1) - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));
idx_on_element = fpl(:,2) >= te_eta & fpl(:,2) <= le_eta & fpl(:,1) >= xi_1 & fpl(:,1) <= xi_3 & abs(fpl(:,3)) <= margin;
idx_on_edge = idx_on_element & ((le_eta - margin <= fpl(:,2) & fpl(:,2) <= le_eta + margin) | (te_eta - margin <= fpl(:,2) & fpl(:,2) <= te_eta + margin));

fpl(idx_on_element,3) = fpl(idx_on_element,3) + 100*margin;
fpl(idx_on_edge,2) = fpl(idx_on_edge,2) + 100*margin;

xi_p = fpl(:,1);
eta_p = fpl(:,2);
zeta_p = fpl(:,3);

syms x y
syms xi eta

d_inf = [...
    0,...
    -eta.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -eta.*(eta-eta_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    ...
    0,...
    -zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    -(eta-eta_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    ...
    -xi.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    0, ...
    -xi.*(xi-xi_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    0.0,...
    -(xi-xi_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -eta.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),-xi.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -(xi.*(eta-eta_p)+eta.*(xi-xi_p)).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    ...
    0.0,0.0,0.0];

d_inf = subs(d_inf, 'eta', 'y');
d_inf = subs(d_inf, 'xi', 'x');

eta_le = eta_2 + (x - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
eta_te = eta_1 + (x - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));

count = 1;
h = waitbar(count/length(d_inf(:)),'Please wait...');
infl_loc = zeros(3,6,size(fpl,1));
for i = 1:size(d_inf,1)
    for j = 1:size(d_inf,2)
        if isempty(find(d_inf(i,j)==0,1))
            infl_loc(mod(j-1,3) + 1, floor(j/3) + 1, i) = integral2(matlabFunction(d_inf(i,j)), xi_1(i), xi_3(i), matlabFunction(eta_te(i)), matlabFunction(eta_le(i)));
        end
        waitbar(count/length(d_inf(:)), h)
        count = count + 1;
    end
end
close(h);

% Only using the tangential velocities for points on the element
infl_loc(1:2,:,idx_on_element) = infl_loc(1:2,:,idx_on_element).*0;

%%
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));
infl_tot(isnan(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);


end
