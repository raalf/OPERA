function [infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER)
warning('on')
fpl = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,:));
len = size(fpl,1);

%%
xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

% Checking which elements are on the element, moving them off by a small
% amount
le_eta = eta_2 + (fpl(:,1) - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
te_eta = eta_1 + (fpl(:,1) - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));

% tester = abs(fpl(:,3)) < 1e-3;
% fpl(tester,3) = 1e-5
margin_edge = 1e-10;
margin_on_element = 1e-10;
idx_on_element = fpl(:,2) >= te_eta - margin_edge & fpl(:,2) <= le_eta + margin_edge & fpl(:,1) >= xi_1 - margin_edge & fpl(:,1) <= xi_3 + margin_edge & abs(fpl(:,3)) <= margin_on_element;

% margin_above_below_element = 1e-3;
% idx_above_below_element = fpl(:,2) >= te_eta - margin_edge & fpl(:,2) <= le_eta + margin_edge & fpl(:,1) >= xi_1 - margin_edge & fpl(:,1) <= xi_3 + margin_edge & abs(fpl(:,3)) > margin_on_element & abs(fpl(:,3)) <= margin_above_below_element;
idx_on_edge = idx_on_element & ((xi_1 - margin_edge <= fpl(:,1) & fpl(:,1) <= xi_1 + margin_edge) | (le_eta - margin_edge <= fpl(:,2) & fpl(:,2) <= le_eta + margin_edge) | (te_eta - margin_edge <= fpl(:,2) & fpl(:,2) <= te_eta + margin_edge));

% fpl(idx_on_element,3) = fpl(idx_on_element,3).*0 + 1e-4;
% fpl(idx_on_edge,:) = fpl(idx_on_edge,:) - margin.*fpl(idx_on_edge,:);
% fpl(idx_on_edge,3) = fpl(idx_on_edge,3) + margin;

xi_p = fpl(:,1);
eta_p = fpl(:,2);
zeta_p = fpl(:,3);

%%
tic
x_m = xi_p;
y_m = eta_p;
z_m = zeta_p;
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

J_1 = fcnJ_1(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE);
J_2 = fcnJ_2(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE);
J_3 = fcnJ_3(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE);
J_4 = fcnJ_4(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE);
J_5 = [];
J_6 = fcnJ_6(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE);
infl_new = zeros(3,6,len);

infl_new(1,3,:) = reshape(J_2.*z_m,1,1,[]);
infl_new(1,4,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(1,5,:) = reshape(J_4.*z_m,1,1,[]);

infl_new(2,1,:) = reshape(J_4.*z_m,1,1,[]);
infl_new(2,2,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(2,5,:) = reshape(J_2.*z_m,1,1,[]);

infl_new(3,1,:) = (reshape(-J_4.*y_m ,1,1,[]));% + reshape(J_5,1,1,[]));
infl_new(3,2,:) = (reshape(-J_1.*y_m ,1,1,[]) + reshape(J_4,1,1,[]));
infl_new(3,3,:) = (reshape(-J_2.*x_m,1,1,[]) + reshape(J_3,1,1,[]));
infl_new(3,4,:) = (reshape(-J_1.*x_m,1,1,[]) + reshape(J_2,1,1,[]));
infl_new(3,5,:) = (reshape(-J_2.*y_m,1,1,[]) - reshape(J_4.*x_m,1,1,[]) +  reshape(2.*J_6,1,1,[]));

infl_new(isnan(infl_new)) = 0;
infl_loc = real(infl_new);

toc
%%
% count = 1;
% D = parallel.pool.DataQueue;
% h = waitbar(0, 'Please wait ...');
% afterEach(D, @nUpdateWaitbar);
% 
% AbsTol = 1e-7;
% RelTol = 1e-7;
% 
% infl_loc = zeros(3,6,len);
% 
% parfor r = 1:len
% % for i = 1:len
%     le_eta = @(x) eta_2(r) + (x - xi_2(r)).*((eta_3(r) - eta_2(r))./(xi_3(r) - xi_2(r)));
%     te_eta = @(x) eta_1(r) + (x - xi_1(r)).*((eta_3(r) - eta_1(r))./(xi_3(r) - xi_1(r)));
%     
%     tmp = zeros(3,6);
%     
%     if idx_on_element(r) == false
%         %{
%         denom = @(x,y) (((abs(zeta_p(r)).^2) + (abs(y - eta_p(r)).^2) + (abs(xi_p(r) - x).^2)).^(3/2));
% %         disp(['Off Element Totally: ', num2str(r)])
%                 
%         % A1
%         term = @(x,y) (y.*zeta_p(r))./denom(x,y);
%         tmp(2,1) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % A2
%         term = @(x,y) zeta_p(r)./denom(x,y);
%         tmp(2,2) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % B1
%         term = @(x,y) (x.*zeta_p(r))./denom(x,y);
%         tmp(1,3) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % B2
%         term = @(x,y) zeta_p(r)./denom(x,y);
%         tmp(1,4) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % C2
%         term = @(x,y) (y.*zeta_p(r))./denom(x,y);
%         tmp(1,5) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         term = @(x,y) (x.*zeta_p(r))./denom(x,y);
%         tmp(2,5) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         
%         % A1
%         term = @(x,y) -(y.*(eta_p(r)-y))./denom(x,y);
%         tmp(3,1) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % A2
%         term = @(x,y) (y - eta_p(r))./denom(x,y);
%         tmp(3,2) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % B1
%         term = @(x,y) -(x.*(xi_p(r)-x))./denom(x,y);
%         tmp(3,3) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % B2
%         term = @(x,y) (x - xi_p(r))./denom(x,y);
%         tmp(3,4) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         % C2
%         term = @(x,y) (-x.*(eta_p(r) - y) - y.*(xi_p(r) - x))./denom(x,y);
%         tmp(3,5) = integral2(term, xi_1(r), xi_3(r), le_eta, te_eta,'AbsTol',AbsTol','RelTol',RelTol);
%         %}
%     elseif idx_on_element(r) == true && idx_on_edge(r) == false
%         denom = @(x,y) ((abs(y-eta_p(r)).^2 + abs(xi_p(r)-x).^2).^(3/2));
% %         margin_edge = 2e-2;
%         margin_edge = 1e-8;
% %         disp(['On Element: ', num2str(r)])
%         % A1
% %         term = @(x,y) -(y.*(eta_p(r)-y))./denom(x,y);        
% %         tmp(3,1) =  integral2(term, xi_1(r), xi_p(r) - margin_edge, le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_1(r), xi_p(r) - margin_edge, eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated');
%         % A2
% %         term = @(x,y) (y - eta_p(r))./denom(x,y);
% %         tmp(3,2) =  integral2(term, xi_1(r), xi_p(r) - margin_edge, le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_1(r), xi_p(r) - margin_edge, eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated');
%         % B1
% %         term = @(x,y) -(x.*(xi_p(r)-x))./denom(x,y);
% %         tmp(3,3) =  integral2(term, xi_1(r), xi_p(r) - margin_edge, le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_1(r), xi_p(r) - margin_edge, eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated');
%         % B2
% %         term = @(x,y) (x - xi_p(r))./denom(x,y);
% %         tmp(3,4) =  integral2(term, xi_1(r), xi_p(r) - margin_edge, le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_1(r), xi_p(r) - margin_edge, eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated');
%         % C2
% %         term = @(x,y) (-x.*(eta_p(r) - y) - y.*(xi_p(r) - x))./denom(x,y);
% %         tmp(3,5) =  integral2(term, xi_1(r), xi_p(r) - margin_edge, le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), le_eta, eta_p(r) + margin_edge,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_1(r), xi_p(r) - margin_edge, eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated') + ...
% %             integral2(term, xi_p(r) + margin_edge, xi_3(r), eta_p(r) - margin_edge, te_eta,'AbsTol',AbsTol','RelTol',RelTol, 'method', 'iterated');
%     end
%     
%     infl_loc(:,:,r) = tmp;
%     
%     send(D, r);
% end
% close(h);
% 
%     function nUpdateWaitbar(~)
%         waitbar(count/len, h);
%         count = count + 1;
%     end
% 
% % Only using the normal velocities for points on the element
% % infl_loc(1:2,:,tester) = infl_loc(1:2,:,tester).*0;

%%
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);

infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));
infl_tot(isnan(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);
% infl_glob(isnan(infl_glob)) = 0;
end

function J_1 = fcnJ_1(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE)
t2 = (C .* x_m);
t4 = abs(z_m);
t5 = t4 .* (t2 + D_LE - y_m);
t7 = x_m + z_m;
t10 = x_m - z_m;
t13 = (C .* t10 + D_LE - y_m) .* (C .* t7 + D_LE - y_m);
t15 = sqrt(-2.*i .* t5 .* C + t13);
t19 = x_m .* E;
t21 = t4 .* (t19 + D_TE - y_m);
t27 = (E .* t10 + D_TE - y_m) .* (E .* t7 + D_TE - y_m);
t29 = sqrt(-2.*i .* t21 .* E + t27);
t35 = sqrt(2.*i .* t5 .* C + t13);
t40 = sqrt(2.*i .* t21 .* E + t27);
t43 = t40 .* t35;
t44 = t4 .* E;
t46 = t15 .* (i .* t44 - t19 - D_TE + y_m);
t47 = -i .* t29;
t48 = (xi_1 .^ 2);
t49 = E .^ 2;
t51 = D_TE - y_m;
t55 = D_TE .^ 2;
t57 = 2 .* D_TE .* y_m;
t58 = x_m .^ 2;
t60 = 2 .* xi_1 .* x_m;
t61 = y_m .^ 2;
t62 = z_m .^ 2;
t64 = sqrt((2 .* E .* xi_1 .* t51 + t49 .* t48 + t48 + t55 - t57 + t58 - t60 + t61 + t62));
t66 = xi_1 .* t49;
t68 = -E .* t51;
t71 = -i .* x_m;
t72 = t66 .* t71;
t74 = -i .* (x_m + xi_1);
t75 = E .* t51;
t76 = t75 .* t74;
t77 = -i .* t55;
t79 = 2.*i .* D_TE .* y_m;
t80 = -i .* t61;
t81 = -i .* t62;
t83 = -i .* xi_1;
t85 = 0.1e1 ./ (t4 + i .* x_m + t83);
t87 = log(t85 .* (t64 .* t47 + (t4 .* (-t66 + t68 + x_m - xi_1)) + t72 + t76 + t77 + t79 + t80 + t81));
t90 = (xi_3 .^ 2);
t96 = 2 .* xi_3 .* x_m;
t98 = sqrt((2 .* E .* t51 .* xi_3 + t90 .* t49 + t55 - t57 + t58 + t61 + t62 + t90 - t96));
t100 = xi_3 .* t49;
t104 = x_m + xi_3;
t105 = -i .* t104;
t108 = -i .* xi_3;
t110 = 0.1e1 ./ (t4 + i .* x_m + t108);
t112 = log(t110 .* (t98 .* t47 + (t4 .* (-t100 + t68 - xi_3 + x_m)) + t100 .* t71 + t75 .* t105 + t77 + t79 + t80 + t81));
t115 = t4 .* C;
t117 = -i .* t15;
t118 = C .^ 2;
t120 = D_LE - y_m;
t124 = D_LE .^ 2;
t126 = 2 .* D_LE .* y_m;
t128 = sqrt((2 .* C .* xi_1 .* t120 + t118 .* t48 + t124 - t126 + t48 + t58 - t60 + t61 + t62));
t130 = xi_1 .* t118;
t132 = -t120 .* C;
t136 = t118 .* x_m .* t83;
t137 = t120 .* C;
t138 = t137 .* t74;
t139 = -i .* t124;
t141 = 2.*i .* D_LE .* y_m;
t144 = log(t85 .* (t128 .* t117 + (t4 .* (-t130 + t132 + x_m - xi_1)) + t136 + t138 + t139 + t141 + t80 + t81));
t150 = sqrt((2 .* C .* t120 .* xi_3 + t90 .* t118 + t124 - t126 + t58 + t61 + t62 + t90 - t96));
t152 = xi_3 .* t118;
t155 = t152 .* t71;
t156 = t137 .* t105;
t159 = log(t110 .* (t150 .* t117 + (t4 .* (-t152 + t132 - xi_3 + x_m)) + t155 + t156 + t139 + t141 + t80 + t81));
t163 = -i .* t35;
t169 = 0.1e1 ./ (i .* x_m + t83 - t4);
t171 = log(t169 .* (t128 .* t163 + (t4 .* (t130 + t137 - x_m + xi_1)) + t136 + t138 + t139 + t141 + t80 + t81));
t177 = 0.1e1 ./ (i .* x_m + t108 - t4);
t179 = log(t177 .* (t150 .* t163 + (t4 .* (t152 + t137 - x_m + xi_3)) + t155 + t156 + t139 + t141 + t80 + t81));
t192 = t55 - t57 + t61 + t62;
t196 = log(t169 .* (-i .* t64 .* t40 + (t4 .* (t66 + t75 - x_m + xi_1)) + t72 + t76 + -i .* t192));
t213 = log(t177 .* (-2.*i .* t40 .* t98 + (t4 .* (2 .* E .* t51 + 2 .* t100 - 2 .* x_m + 2 .* xi_3)) + -2.*i .* t100 .* x_m + -2.*i .* t75 .* t104 + -2.*i .* t192));
t214 = log(0.2e1);
J_1 = -0.1e1 ./ 0.2e1.*i ./ t4 .* (t87 .* t46 .* t43 - t112 .* t46 .* t43 + t29 .* (t40 .* (-t35 .* (t144 - t159) .* (i .* t115 - t2 - D_LE + y_m) - t15 .* (i .* t115 + t2 + D_LE - y_m) .* (t171 - t179)) + t15 .* (t196 - t213 + t214) .* (i .* t44 + t19 + D_TE - y_m) .* t35)) ./ t40 ./ t35 ./ t29 ./ t15;
end

function J_2 = fcnJ_2(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE)
t2 = (C .* x_m);
t4 = abs(z_m);
t5 = t4 .* (t2 + D_LE - y_m);
t7 = x_m + z_m;
t10 = x_m - z_m;
t13 = (C .* t10 + D_LE - y_m) .* (C .* t7 + D_LE - y_m);
t15 = sqrt(2.*i .* t5 .* C + t13);
t21 = sqrt(-2.*i .* t5 .* C + t13);
t24 = x_m .* E;
t26 = t4 .* (t24 + D_TE - y_m);
t32 = (E .* t10 + D_TE - y_m) .* (E .* t7 + D_TE - y_m);
t34 = sqrt(-2.*i .* t26 .* E + t32);
t41 = sqrt(2.*i .* t26 .* E + t32);
t44 = -0.1e1 ./ 0.2e1.*i .* y_m;
t46 = t4 .* (i .* t24 + 0.1e1 ./ 0.2e1.*i .* D_TE + t44);
t47 = x_m .^ 2;
t49 = (t47 .* E) ./ 0.2e1;
t50 = -D_TE + y_m;
t51 = t50 ./ 0.2e1;
t53 = z_m .^ 2;
t55 = (t53 .* E) ./ 0.2e1;
t58 = C .^ 2;
t60 = sqrt((t58 + 1));
t61 = t60 .* ((x_m .* t51) + t46 - t49 + t55) .* t21;
t62 = E .^ 2;
t64 = sqrt((t62 + 1));
t65 = t41 .* t64;
t66 = -i .* t34;
t67 = (xi_1 .^ 2);
t69 = -t50;
t73 = D_TE .^ 2;
t75 = 2 .* D_TE .* y_m;
t77 = 2 .* xi_1 .* x_m;
t78 = y_m .^ 2;
t80 = sqrt((2 .* E .* xi_1 .* t69 + t62 .* t67 + t47 + t53 + t67 + t73 - t75 - t77 + t78));
t82 = xi_1 .* t62;
t83 = E .* t50;
t86 = -i .* x_m;
t87 = t82 .* t86;
t89 = -i .* (x_m + xi_1);
t90 = E .* t69;
t91 = t90 .* t89;
t92 = -i .* t73;
t94 = 2.*i .* D_TE .* y_m;
t95 = -i .* t78;
t96 = -i .* t53;
t98 = -i .* xi_1;
t100 = 0.1e1 ./ (t4 + i .* x_m + t98);
t102 = log(t100 .* (t80 .* t66 + (t4 .* (-t82 + t83 + x_m - xi_1)) + t87 + t91 + t92 + t94 + t95 + t96));
t106 = (xi_3 .^ 2);
t112 = 2 .* xi_3 .* x_m;
t114 = sqrt((2 .* E .* t69 .* xi_3 + t106 .* t62 + t106 - t112 + t47 + t53 + t73 - t75 + t78));
t116 = xi_3 .* t62;
t120 = x_m + xi_3;
t121 = -i .* t120;
t124 = -i .* xi_3;
t126 = 0.1e1 ./ (t4 + i .* x_m + t124);
t128 = log(t126 .* (t114 .* t66 + (t4 .* (-t116 + t83 - xi_3 + x_m)) + t116 .* t86 + t90 .* t121 + t92 + t94 + t95 + t96));
t133 = D_LE - y_m;
t137 = D_LE .^ 2;
t139 = 2 .* D_LE .* y_m;
t141 = sqrt((2 .* C .* xi_1 .* t133 + t58 .* t67 + t137 - t139 + t47 + t53 + t67 - t77 + t78));
t143 = xi_1 .* t58;
t144 = t133 .* C;
t146 = 0.1e1 ./ t60;
t148 = log(t146 .* (t60 .* t141 + t143 + t144 - x_m + xi_1));
t154 = sqrt((2 .* C .* t133 .* xi_3 + t106 .* t58 + t106 - t112 + t137 - t139 + t47 + t53 + t78));
t156 = xi_3 .* t58;
t159 = log(t146 .* (t60 .* t154 + t144 + t156 - x_m + xi_3));
t165 = 0.1e1 ./ t64;
t167 = log(t165 .* (t64 .* t80 + t82 + t90 - x_m + xi_1));
t171 = log(t165 .* (t64 .* t114 + t116 + t90 - x_m + xi_3));
t180 = t4 .* (i .* t2 + t44 + 0.1e1 ./ 0.2e1.*i .* D_LE);
t182 = (t47 .* C) ./ 0.2e1;
t183 = -t133;
t184 = t183 ./ 0.2e1;
t187 = (t53 .* C) ./ 0.2e1;
t190 = -i .* t21;
t192 = C .* t183;
t196 = t58 .* x_m .* t98;
t197 = t144 .* t89;
t198 = -i .* t137;
t200 = 2.*i .* D_LE .* y_m;
t203 = log(t100 .* (t141 .* t190 + (t4 .* (-t143 + t192 + x_m - xi_1)) + t196 + t197 + t198 + t200 + t95 + t96));
t207 = t156 .* t86;
t208 = t144 .* t121;
t211 = log(t126 .* (t154 .* t190 + (t4 .* (-t156 + t192 - xi_3 + x_m)) + t207 + t208 + t198 + t200 + t95 + t96));
t222 = -i .* t15;
t228 = 0.1e1 ./ (i .* x_m + t98 - t4);
t230 = log(t228 .* (t141 .* t222 + (t4 .* (t143 + t144 - x_m + xi_1)) + t196 + t197 + t198 + t200 + t95 + t96));
t236 = 0.1e1 ./ (i .* x_m + t124 - t4);
t238 = log(t236 .* (t154 .* t222 + (t4 .* (t156 + t144 - x_m + xi_3)) + t207 + t208 + t198 + t200 + t95 + t96));
t263 = t73 - t75 + t78 + t53;
t267 = log(t236 .* (-2.*i .* t41 .* t114 + (t4 .* (2 .* E .* t69 + 2 .* t116 - 2 .* x_m + 2 .* xi_3)) + -2.*i .* t116 .* x_m + -2.*i .* t90 .* t120 + -2.*i .* t263));
t275 = log(t228 .* (-i .* t80 .* t41 + (t4 .* (t82 + t90 - x_m + xi_1)) + t87 + t91 + -i .* t263));
t276 = log(0.2e1);
J_2 = -i ./ t4 .* t146 .* t165 .* (t102 .* t15 .* t65 .* t61 - t128 .* t15 .* t65 .* t61 + (t41 .* (t15 .* (i .* t21 .* (t64 .* (t148 - t159) .* C - E .* t60 .* (t167 - t171)) .* t4 - (t203 - t211) .* t64 .* t60 .* ((x_m .* t184) + t180 - t182 + t187)) - (t230 - t238) .* (-(x_m .* t184) + t180 + t182 - t187) .* t64 .* t60 .* t21) - (t267 - t275 - t276) .* t15 .* t64 .* t60 .* (-(x_m .* t51) + t46 + t49 - t55) .* t21) .* t34) ./ t41 ./ t34 ./ t21 ./ t15;
end

function J_3 = fcnJ_3(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE)
t4 = abs(z_m);
t5 = t4 .* (x_m .* E + D_TE - y_m);
t7 = x_m + z_m;
t10 = x_m - z_m;
t13 = (E .* t10 + D_TE - y_m) .* (E .* t7 + D_TE - y_m);
t15 = sqrt(-2.*i .* t5 .* E + t13);
t21 = sqrt(2.*i .* t5 .* E + t13);
t23 = (C .^ 2);
t24 = t23 + 1;
t25 = sqrt(t24);
t26 = t25 .* t24;
t33 = t4 .* (C .* x_m + D_LE - y_m);
t39 = (C .* t10 + D_LE - y_m) .* (C .* t7 + D_LE - y_m);
t41 = sqrt(-2.*i .* t33 .* C + t39);
t43 = E .^ 2;
t44 = t43 + 1;
t45 = sqrt(t44);
t46 = t45 .* t44;
t49 = t26 .* t21;
t50 = t46 .* t49;
t54 = sqrt(2.*i .* t33 .* C + t39);
t55 = t15 .* t54;
t56 = x_m .^ 2;
t59 = -0.2e1 ./ 0.3e1.*i .* y_m;
t63 = z_m .^ 2;
t66 = t4 .* (i .* t56 .* C + x_m .* (0.2e1 ./ 0.3e1.*i .* D_LE + t59) + -0.1e1 ./ 0.3e1.*i .* t63 .* C);
t67 = x_m .* t56;
t69 = (t67 .* C) ./ 0.3e1;
t70 = -D_LE + y_m;
t71 = t70 ./ 0.3e1;
t74 = x_m .* t63 .* C;
t75 = -t70;
t77 = (t75 .* t63) ./ 0.3e1;
t78 = (t56 .* t71) + t66 - t69 + t74 + t77;
t79 = -i .* t41;
t80 = (xi_1 .^ 2);
t85 = D_LE .^ 2;
t87 = 2 .* D_LE .* y_m;
t89 = 2 .* xi_1 .* x_m;
t90 = y_m .^ 2;
t92 = sqrt((2 .* C .* xi_1 .* t75 + t23 .* t80 + t56 + t63 + t80 + t85 - t87 - t89 + t90));
t94 = xi_1 .* t23;
t95 = C .* t70;
t98 = -i .* xi_1;
t100 = t23 .* x_m .* t98;
t102 = -i .* (x_m + xi_1);
t103 = t75 .* C;
t104 = t103 .* t102;
t105 = -i .* t85;
t107 = 2.*i .* D_LE .* y_m;
t108 = -i .* t90;
t109 = -i .* t63;
t112 = 0.1e1 ./ (t4 + i .* x_m + t98);
t114 = log(t112 .* (t92 .* t79 + (t4 .* (-t94 + t95 + x_m - xi_1)) + t100 + t104 + t105 + t107 + t108 + t109));
t118 = (xi_3 .^ 2);
t124 = 2 .* xi_3 .* x_m;
t126 = sqrt((2 .* C .* t75 .* xi_3 + t118 .* t23 + t118 - t124 + t56 + t63 + t85 - t87 + t90));
t128 = xi_3 .* t23;
t131 = -i .* x_m;
t132 = t128 .* t131;
t133 = x_m + xi_3;
t134 = -i .* t133;
t135 = t103 .* t134;
t137 = -i .* xi_3;
t139 = 0.1e1 ./ (t4 + i .* x_m + t137);
t141 = log(t139 .* (t126 .* t79 + (t4 .* (-t128 + t95 - xi_3 + x_m)) + t132 + t135 + t105 + t107 + t108 + t109));
t149 = t26 .* (-(t56 .* t71) + t66 + t69 - t74 - t77) .* t21;
t150 = t15 .* t46;
t151 = -i .* t54;
t157 = 0.1e1 ./ (i .* x_m + t98 - t4);
t159 = log(t157 .* (t92 .* t151 + (t4 .* (t94 + t103 - x_m + xi_1)) + t100 + t104 + t105 + t107 + t108 + t109));
t167 = 0.1e1 ./ (i .* x_m + t137 - t4);
t169 = log(t167 .* (t126 .* t151 + (t4 .* (t128 + t103 - x_m + xi_3)) + t132 + t135 + t105 + t107 + t108 + t109));
t179 = t4 .* (i .* t56 .* E + x_m .* (0.2e1 ./ 0.3e1.*i .* D_TE + t59) + -0.1e1 ./ 0.3e1.*i .* t63 .* E);
t181 = (t67 .* E) ./ 0.3e1;
t182 = D_TE - y_m;
t183 = t182 ./ 0.3e1;
t186 = x_m .* t63 .* E;
t188 = (t182 .* t63) ./ 0.3e1;
t189 = (t56 .* t183) + t179 + t181 - t186 - t188;
t190 = t189 .* t26;
t195 = D_TE .^ 2;
t197 = 2 .* D_TE .* y_m;
t199 = sqrt((2 .* E .* t182 .* xi_3 + t118 .* t43 + t118 - t124 + t195 - t197 + t56 + t63 + t90));
t202 = xi_3 .* t43;
t213 = E .* t182;
t215 = t195 - t197 + t90 + t63;
t219 = log(t167 .* (-2.*i .* t21 .* t199 + (t4 .* (2 .* E .* t182 + 2 .* t202 - 2 .* x_m + 2 .* xi_3)) + -2.*i .* t202 .* x_m + -2.*i .* t213 .* t133 + -2.*i .* t215));
t223 = -i .* t21;
t229 = sqrt((2 .* E .* xi_1 .* t182 + t43 .* t80 + t195 - t197 + t56 + t63 + t80 - t89 + t90));
t231 = xi_1 .* t43;
t234 = t231 .* t131;
t235 = t213 .* t102;
t239 = log(t157 .* (t229 .* t223 + (t4 .* (t231 + t213 - x_m + xi_1)) + t234 + t235 + -i .* t215));
t246 = t46 .* (-(t56 .* t183) + t179 - t181 + t186 + t188);
t247 = -i .* t15;
t250 = -E .* t182;
t253 = -i .* t195;
t255 = 2.*i .* D_TE .* y_m;
t258 = log(t112 .* (t229 .* t247 + (t4 .* (-t231 + t250 + x_m - xi_1)) + t234 + t235 + t253 + t255 + t108 + t109));
t269 = log(t139 .* (t199 .* t247 + (t4 .* (-t202 + t250 - xi_3 + x_m)) + t202 .* t131 + t213 .* t134 + t253 + t255 + t108 + t109));
t279 = y_m ./ 0.2e1;
t281 = (x_m .* ((C .* t23) + 0.3e1 ./ 0.2e1 .* C) + D_LE ./ 0.2e1 - t279) .* t46;
t284 = 0.1e1 ./ t25;
t286 = log(t284 .* (t25 .* t92 + t103 + t94 - x_m + xi_1));
t295 = (x_m .* ((E .* t43) + 0.3e1 ./ 0.2e1 .* E) + D_TE ./ 0.2e1 - t279) .* t4;
t298 = 0.1e1 ./ t45;
t300 = log(t298 .* (t45 .* t229 + t213 + t231 - x_m + xi_1));
t307 = log(t284 .* (t25 .* t126 + t103 + t128 - x_m + xi_3));
t315 = log(t298 .* (t45 .* t199 + t202 + t213 - x_m + xi_3));
t320 = t44 .* C;
t323 = -0.2e1 ./ 0.3e1.*i .* t21;
t336 = log(0.2e1);
J_3 = -0.3e1 ./ 0.2e1.*i ./ t4 ./ t54 .* (-t114 .* t78 .* t55 .* t50 + t141 .* t78 .* t55 .* t50 + (-t159 .* t150 .* t149 + t169 .* t150 .* t149 - (0.4e1 ./ 0.3e1) .* t54 .* ((0.3e1 ./ 0.4e1) .* t219 .* t150 .* t190 - (0.3e1 ./ 0.4e1) .* t239 .* t150 .* t190 - (0.3e1 ./ 0.4e1) .* t258 .* t246 .* t49 + (0.3e1 ./ 0.4e1) .* t269 .* t246 .* t49 + t15 .* (t286 .* t281 .* t4 .* t223 + i .* t300 .* t295 .* t49 + i .* t307 .* t281 .* t4 .* t21 - (0.3e1 ./ 0.4e1) .* t25 .* (0.4e1 ./ 0.3e1.*i .* t315 .* t295 .* t24 .* t21 + (0.2e1 ./ 0.3e1.*i .* t92 .* t320 .* t4 .* t21 + t229 .* E .* t4 .* t24 .* t323 + t126 .* t320 .* t4 .* t323 + (0.2e1 ./ 0.3e1.*i .* t4 .* t21 .* t199 .* E + t336 .* t44 .* t189) .* t24) .* t45)))) .* t41) ./ t46 ./ t41 ./ t26 ./ t21 ./ t15;
end

function J_4 = fcnJ_4(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE)
t2 = (C .* x_m);
t4 = abs(z_m);
t5 = t4 .* (t2 + D_LE - y_m);
t7 = x_m + z_m;
t10 = x_m - z_m;
t13 = (C .* t10 + D_LE - y_m) .* (C .* t7 + D_LE - y_m);
t15 = sqrt(2.*i .* t5 .* C + t13);
t21 = sqrt(-2.*i .* t5 .* C + t13);
t22 = C .^ 2;
t24 = sqrt((t22 + 1));
t26 = (E .^ 2);
t28 = sqrt((t26 + 1));
t30 = y_m .* t28 .* t24 .* t21;
t32 = x_m .* E;
t34 = t4 .* (t32 + D_TE - y_m);
t40 = (E .* t10 + D_TE - y_m) .* (E .* t7 + D_TE - y_m);
t42 = sqrt(2.*i .* t34 .* E + t40);
t43 = t15 .* t42;
t44 = t4 .* E;
t45 = i .* t44 - t32 - D_TE + y_m;
t49 = sqrt(-2.*i .* t34 .* E + t40);
t50 = -i .* t49;
t51 = (xi_1 .^ 2);
t53 = D_TE - y_m;
t57 = D_TE .^ 2;
t59 = 2 .* D_TE .* y_m;
t60 = x_m .^ 2;
t62 = 2 .* xi_1 .* x_m;
t63 = y_m .^ 2;
t64 = z_m .^ 2;
t66 = sqrt((2 .* E .* xi_1 .* t53 + t26 .* t51 + t51 + t57 - t59 + t60 - t62 + t63 + t64));
t68 = xi_1 .* t26;
t70 = -E .* t53;
t73 = -i .* x_m;
t74 = t68 .* t73;
t76 = -i .* (x_m + xi_1);
t77 = E .* t53;
t78 = t77 .* t76;
t79 = -i .* t57;
t81 = 2.*i .* D_TE .* y_m;
t82 = -i .* t63;
t83 = -i .* t64;
t85 = -i .* xi_1;
t87 = 0.1e1 ./ (t4 + i .* x_m + t85);
t89 = log(t87 .* (t66 .* t50 + (t4 .* (-t68 + t70 + x_m - xi_1)) + t74 + t78 + t79 + t81 + t82 + t83));
t93 = (xi_3 .^ 2);
t99 = 2 .* xi_3 .* x_m;
t101 = sqrt((2 .* E .* t53 .* xi_3 + t93 .* t26 + t57 - t59 + t60 + t63 + t64 + t93 - t99));
t103 = xi_3 .* t26;
t107 = x_m + xi_3;
t108 = -i .* t107;
t111 = -i .* xi_3;
t113 = 0.1e1 ./ (t4 + i .* x_m + t111);
t115 = log(t113 .* (t101 .* t50 + (t4 .* (-t103 + t70 - xi_3 + x_m)) + t103 .* t73 + t77 .* t108 + t79 + t81 + t82 + t83));
t121 = D_LE - y_m;
t125 = D_LE .^ 2;
t127 = 2 .* D_LE .* y_m;
t129 = sqrt((2 .* C .* xi_1 .* t121 + t22 .* t51 + t125 - t127 + t51 + t60 - t62 + t63 + t64));
t131 = xi_1 .* t22;
t132 = t121 .* C;
t134 = 0.1e1 ./ t24;
t136 = log(t134 .* (t24 .* t129 + t131 + t132 - x_m + xi_1));
t142 = sqrt((2 .* C .* t121 .* xi_3 + t93 .* t22 + t125 - t127 + t60 + t63 + t64 + t93 - t99));
t144 = xi_3 .* t22;
t147 = log(t134 .* (t24 .* t142 + t132 + t144 - x_m + xi_3));
t152 = 0.1e1 ./ t28;
t154 = log(t152 .* (t28 .* t66 + t68 + t77 - x_m + xi_1));
t158 = log(t152 .* (t28 .* t101 + t103 + t77 - x_m + xi_3));
t164 = t4 .* C;
t167 = -i .* t21;
t170 = -t121 .* C;
t174 = t22 .* x_m .* t85;
t175 = t132 .* t76;
t176 = -i .* t125;
t178 = 2.*i .* D_LE .* y_m;
t181 = log(t87 .* (t129 .* t167 + (t4 .* (-t131 + t170 + x_m - xi_1)) + t174 + t175 + t176 + t178 + t82 + t83));
t185 = t144 .* t73;
t186 = t132 .* t108;
t189 = log(t113 .* (t142 .* t167 + (t4 .* (-t144 + t170 - xi_3 + x_m)) + t185 + t186 + t176 + t178 + t82 + t83));
t196 = -i .* t15;
t202 = 0.1e1 ./ (i .* x_m + t85 - t4);
t204 = log(t202 .* (t129 .* t196 + (t4 .* (t131 + t132 - x_m + xi_1)) + t174 + t175 + t176 + t178 + t82 + t83));
t210 = 0.1e1 ./ (i .* x_m + t111 - t4);
t212 = log(t210 .* (t142 .* t196 + (t4 .* (t144 + t132 - x_m + xi_3)) + t185 + t186 + t176 + t178 + t82 + t83));
t235 = t57 - t59 + t63 + t64;
t239 = log(t210 .* (-2.*i .* t42 .* t101 + (t4 .* (2 .* E .* t53 + 2 .* t103 - 2 .* x_m + 2 .* xi_3)) + -2.*i .* t103 .* x_m + -2.*i .* t77 .* t107 + -2.*i .* t235));
t247 = log(t202 .* (-i .* t66 .* t42 + (t4 .* (t68 + t77 - x_m + xi_1)) + t74 + t78 + -i .* t235));
t248 = log(0.2e1);
J_4 = 0.1e1 ./ 0.2e1.*i ./ t4 .* t152 ./ t42 ./ t49 .* t134 ./ t21 .* (-t89 .* t45 .* t43 .* t30 + t115 .* t45 .* t43 .* t30 + t49 .* (t42 .* (t15 .* (2.*i .* t21 .* (t28 .* (t136 - t147) - (t154 - t158) .* t24) .* t4 + y_m .* (t181 - t189) .* t28 .* t24 .* (i .* t164 - t2 - D_LE + y_m)) + y_m .* (i .* t164 + t2 + D_LE - y_m) .* t28 .* t24 .* (t204 - t212) .* t21) + t15 .* y_m .* (i .* t44 + t32 + D_TE - y_m) .* t28 .* t24 .* t21 .* (t239 - t247 - t248))) ./ t15;
end

function J_6 = fcnJ_6(x_m, y_m, z_m, xi_1, xi_3, C, D_LE, E, D_TE)
t2 = (x_m .* E);
t4 = abs(z_m);
t5 = t4 .* (t2 + D_TE - y_m);
t7 = x_m + z_m;
t10 = x_m - z_m;
t13 = (E .* t10 + D_TE - y_m) .* (E .* t7 + D_TE - y_m);
t15 = sqrt(2.*i .* t5 .* E + t13);
t19 = C .* x_m;
t21 = t4 .* (t19 + D_LE - y_m);
t27 = (C .* t10 + D_LE - y_m) .* (C .* t7 + D_LE - y_m);
t29 = sqrt(-2.*i .* t21 .* C + t27);
t31 = E .^ 2;
t32 = t31 + 1;
t33 = sqrt(t32);
t34 = t33 .* t32;
t41 = sqrt(2.*i .* t21 .* C + t27);
t43 = C .^ 2;
t44 = t43 + 1;
t45 = sqrt(t44);
t46 = t45 .* t44;
t52 = sqrt(-2.*i .* t5 .* E + t13);
t54 = t46 .* t52;
t56 = t41 .* y_m .* t54;
t57 = t15 .* t34;
t58 = -0.1e1 ./ 0.2e1.*i .* y_m;
t61 = t4 .* (i .* t19 + t58 + 0.1e1 ./ 0.2e1.*i .* D_LE);
t63 = (x_m .* y_m) ./ 0.2e1;
t64 = x_m .^ 2;
t65 = z_m .^ 2;
t67 = -t64 ./ 0.2e1 + t65 ./ 0.2e1;
t70 = (x_m .* D_LE) ./ 0.2e1;
t71 = (C .* t67) + t61 + t63 - t70;
t72 = -i .* t29;
t73 = (xi_1 .^ 2);
t75 = D_LE - y_m;
t79 = D_LE .^ 2;
t81 = 2 .* D_LE .* y_m;
t83 = 2 .* xi_1 .* x_m;
t84 = y_m .^ 2;
t86 = sqrt((2 .* C .* xi_1 .* t75 + t43 .* t73 + t64 + t65 + t73 + t79 - t81 - t83 + t84));
t88 = xi_1 .* t43;
t90 = -t75 .* C;
t93 = -i .* xi_1;
t95 = t43 .* x_m .* t93;
t97 = -i .* (x_m + xi_1);
t98 = t75 .* C;
t99 = t98 .* t97;
t100 = -i .* t79;
t102 = 2.*i .* D_LE .* y_m;
t103 = -i .* t84;
t104 = -i .* t65;
t107 = 0.1e1 ./ (t4 + i .* x_m + t93);
t109 = log(t107 .* (t86 .* t72 + (t4 .* (-t88 + t90 + x_m - xi_1)) + t95 + t99 + t100 + t102 + t103 + t104));
t113 = (xi_3 .^ 2);
t119 = 2 .* xi_3 .* x_m;
t121 = sqrt((2 .* C .* t75 .* xi_3 + t113 .* t43 + t113 - t119 + t64 + t65 + t79 - t81 + t84));
t123 = xi_3 .* t43;
t126 = -i .* x_m;
t127 = t123 .* t126;
t128 = x_m + xi_3;
t129 = -i .* t128;
t130 = t98 .* t129;
t132 = -i .* xi_3;
t134 = 0.1e1 ./ (t4 + i .* x_m + t132);
t136 = log(t134 .* (t121 .* t72 + (t4 .* (-t123 + t90 - xi_3 + x_m)) + t127 + t130 + t100 + t102 + t103 + t104));
t140 = y_m .* t54;
t141 = -t67;
t144 = ((C .* t141) + t61 - t63 + t70) .* t34;
t145 = -i .* t41;
t151 = 0.1e1 ./ (i .* x_m + t93 - t4);
t153 = log(t151 .* (t86 .* t145 + (t4 .* (t88 + t98 - x_m + xi_1)) + t95 + t99 + t100 + t102 + t103 + t104));
t162 = 0.1e1 ./ (i .* x_m + t132 - t4);
t164 = log(t162 .* (t121 .* t145 + (t4 .* (t123 + t98 - x_m + xi_3)) + t127 + t130 + t100 + t102 + t103 + t104));
t170 = t4 .* (i .* t2 + 0.1e1 ./ 0.2e1.*i .* D_TE + t58);
t173 = (D_TE .* x_m) ./ 0.2e1;
t174 = (E .* t141) + t170 + t173 - t63;
t175 = t174 .* t34;
t177 = D_TE - y_m;
t181 = D_TE .^ 2;
t183 = 2 .* D_TE .* y_m;
t185 = sqrt((2 .* E .* t177 .* xi_3 + t113 .* t31 + t113 - t119 + t181 - t183 + t64 + t65 + t84));
t188 = xi_3 .* t31;
t199 = E .* t177;
t201 = t181 - t183 + t84 + t65;
t205 = log(t162 .* (-2.*i .* t15 .* t185 + (t4 .* (2 .* E .* t177 + 2 .* t188 - 2 .* x_m + 2 .* xi_3)) + -2.*i .* t188 .* x_m + -2.*i .* t199 .* t128 + -2.*i .* t201));
t214 = sqrt((2 .* E .* xi_1 .* t177 + t31 .* t73 + t181 - t183 + t64 + t65 + t73 - t83 + t84));
t216 = xi_1 .* t31;
t219 = t216 .* t126;
t220 = t199 .* t97;
t224 = log(t151 .* (-i .* t214 .* t15 + (t4 .* (t216 + t199 - x_m + xi_1)) + t219 + t220 + -i .* t201));
t228 = t34 .* y_m .* t46;
t231 = ((E .* t67) + t170 - t173 + t63) .* t15;
t232 = -i .* t52;
t235 = -E .* t177;
t238 = -i .* t181;
t240 = 2.*i .* D_TE .* y_m;
t243 = log(t107 .* (t214 .* t232 + (t4 .* (-t216 + t235 + x_m - xi_1)) + t219 + t220 + t238 + t240 + t103 + t104));
t253 = log(t134 .* (t185 .* t232 + (t4 .* (-t188 + t235 - xi_3 + x_m)) + t188 .* t126 + t199 .* t129 + t238 + t240 + t103 + t104));
t256 = -i .* t4;
t260 = y_m .* C .* t43 + D_LE .* C - x_m;
t264 = 0.1e1 ./ t45;
t266 = log(t264 .* (t45 .* t86 + t88 + t98 - x_m + xi_1));
t272 = y_m .* E .* t31 + D_TE .* E - x_m;
t274 = t15 .* t4;
t277 = 0.1e1 ./ t33;
t279 = log(t277 .* (t33 .* t214 + t199 + t216 - x_m + xi_1));
t286 = log(t264 .* (t45 .* t121 + t123 + t98 - x_m + xi_3));
t293 = log(t277 .* (t33 .* t185 + t188 + t199 - x_m + xi_3));
t308 = log(0.2e1);
J_6 = -i ./ t4 .* (-t109 .* t71 .* t57 .* t56 + t136 .* t71 .* t57 .* t56 + (-t153 .* t15 .* t144 .* t140 + t164 .* t15 .* t144 .* t140 - (t205 .* t175 .* t140 - t224 .* t175 .* t140 - t243 .* t231 .* t228 + t253 .* t231 .* t228 + (t266 .* t57 .* t260 .* t256 + i .* t279 .* t274 .* t272 .* t46 + i .* t286 .* t57 .* t260 .* t4 - t45 .* (i .* t293 .* t274 .* t272 .* t44 + t33 .* (t86 .* t15 .* t32 .* t256 + i .* t214 .* t15 .* t4 .* t44 + i .* t121 .* t15 .* t32 .* t4 + (-i .* t274 .* t185 + t308 .* t174 .* t32 .* y_m) .* t44))) .* t52) .* t41) .* t29) ./ t52 ./ t46 ./ t41 ./ t34 ./ t29 ./ t15;
end

