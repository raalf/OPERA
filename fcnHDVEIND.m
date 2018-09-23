function [infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER)
warning('on')
fpl = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,:));

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

len = size(fpl,1);

count = 1;
D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);

AbsTol = 1e-9;
RelTol = 1e-9;

infl_loc = zeros(3,6,len);
parfor i = 1:len
% for i = 1:len
    le_eta = @(x) eta_2(i) + (x - xi_2(i)).*((eta_3(i) - eta_2(i))./(xi_3(i) - xi_2(i)));
    te_eta = @(x) eta_1(i) + (x - xi_1(i)).*((eta_3(i) - eta_1(i))./(xi_3(i) - xi_1(i)));
    
    tmp = zeros(3,6);
    
    if idx_on_element(i) == false
        denom = @(x,y) ((abs(zeta_p(i)).^2 + abs(y-eta_p(i)).^2 + abs(xi_p(i)-x).^2).^(3/2));
%         disp(['Off Element Totally: ', num2str(i)])
        % A1
        term = @(x,y) (y.*zeta_p(i))./denom(x,y);
        tmp(2,1) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % A2
        term = @(x,y) zeta_p(i)./denom(x,y);
        tmp(2,2) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B1
        term = @(x,y) (x.*zeta_p(i))./denom(x,y);
        tmp(1,3) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B2
        term = @(x,y) zeta_p(i)./denom(x,y);
        tmp(1,4) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % C2
        term = @(x,y) (y.*zeta_p(i))./denom(x,y);
        tmp(1,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        term = @(x,y) (x.*zeta_p(i))./denom(x,y);
        tmp(2,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        
        % A1
        term = @(x,y) -(y.*(eta_p(i)-y))./denom(x,y);
        tmp(3,1) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % A2
        term = @(x,y) (y - eta_p(i))./denom(x,y);
        tmp(3,2) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B1
        term = @(x,y) -(x.*(xi_p(i)-x))./denom(x,y);
        tmp(3,3) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B2
        term = @(x,y) (x - xi_p(i))./denom(x,y);
        tmp(3,4) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % C2
        term = @(x,y) (-x.*(eta_p(i) - y) - y.*(xi_p(i) - x))./denom(x,y);
        tmp(3,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        
    elseif idx_on_element(i) == true && idx_on_edge(i) == false
        denom = @(x,y) ((abs(y-eta_p(i)).^2 + abs(xi_p(i)-x).^2).^(3/2));
        margin_edge = 2e-10;
%         disp(['On Element: ', num2str(i)])
        % A1
        term = @(x,y) -(y.*(eta_p(i)-y))./denom(x,y);
        tmp(3,1) =  integral2(term, xi_1(i), xi_p(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_1(i), xi_p(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % A2
        term = @(x,y) (y - eta_p(i))./denom(x,y);
        tmp(3,2) =  integral2(term, xi_1(i), xi_p(i), te_eta, eta_p(i) - margin_edge,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), te_eta, eta_p(i) - margin_edge,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_1(i), xi_p(i), eta_p(i) + margin_edge, le_eta,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), eta_p(i) + margin_edge, le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B1
        term = @(x,y) -(x.*(xi_p(i)-x))./denom(x,y);
        tmp(3,3) =  integral2(term, xi_1(i), xi_p(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_1(i), xi_p(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % B2
        term = @(x,y) (x - xi_p(i))./denom(x,y);
        tmp(3,4) =  integral2(term, xi_1(i), xi_p(i) - margin_edge, te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i) + margin_edge, xi_3(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_1(i), xi_p(i) - margin_edge, eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i) + margin_edge, xi_3(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol);
        % C2
        term = @(x,y) (-x.*(eta_p(i) - y) - y.*(xi_p(i) - x))./denom(x,y);
        tmp(3,5) =  integral2(term, xi_1(i), xi_p(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), te_eta, eta_p(i),'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_1(i), xi_p(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol) + ...
            integral2(term, xi_p(i), xi_3(i), eta_p(i), le_eta,'AbsTol',AbsTol','RelTol',RelTol);
    end
    
    infl_loc(:,:,i) = tmp;
    
    send(D, i);
end
close(h);

    function nUpdateWaitbar(~)
        waitbar(count/len, h);
        count = count + 1;
    end

% Only using the normal velocities for points on the element
% infl_loc(1:2,:,tester) = infl_loc(1:2,:,tester).*0;

%%
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);

infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));
infl_tot(isnan(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);


end

