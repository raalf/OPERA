function [infl_glob] = fcnHDVEINDFS(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER)

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
margin = 1e-3;
le_eta = eta_2 + (fpl(:,1) - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
te_eta = eta_1 + (fpl(:,1) - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));
idx_on_element = fpl(:,2) >= te_eta - margin & fpl(:,2) <= le_eta + margin & fpl(:,1) >= xi_1 - margin & fpl(:,1) <= xi_3 + margin & abs(fpl(:,3)) <= margin;
idx_on_edge = idx_on_element & ((xi_1 - margin <= fpl(:,1) & fpl(:,1) <= xi_1 + margin) | (le_eta - margin <= fpl(:,2) & fpl(:,2) <= le_eta + margin) | (te_eta - margin <= fpl(:,2) & fpl(:,2) <= te_eta + margin));

fpl(idx_on_element,3) = fpl(idx_on_element,3) + margin;
% fpl(idx_on_edge,:) = fpl(idx_on_edge,:) - margin.*fpl(idx_on_edge,:);
fpl(idx_on_edge,3) = fpl(idx_on_edge,3) + margin;

xi_p = fpl(:,1);
eta_p = fpl(:,2);
zeta_p = fpl(:,3);

len = size(fpl,1);

count = 1;
D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);

infl_loc = zeros(3,6,len);
parfor i = 1:len
% for i = 1:len

    le_eta = @(x) eta_2(i) + (x - xi_2(i)).*((eta_3(i) - eta_2(i))./(xi_3(i) - xi_2(i)));
    te_eta = @(x) eta_1(i) + (x - xi_1(i)).*((eta_3(i) - eta_1(i))./(xi_3(i) - xi_1(i)));
    
    try 
        tmp = zeros(3,6);
        for j = 1:15
            % a1
            if j == 2
                term = @(x,y) -y.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(2,1) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);            
            elseif j == 3
                term = @(x,y) -y.*(y-eta_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(3,1) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
                % a2
            elseif j == 5
                term = @(x,y) -zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(2,2) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            elseif j == 6
                term = @(x,y) -(y-eta_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(3,2) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
                % b1
            elseif j == 7
                term = @(x,y) -x.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(1,3) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            elseif j == 9
                term = @(x,y) -x.*(x-xi_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(3,3) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
                % b2
            elseif j == 10
                term = @(x,y) -zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(1,4) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            elseif j == 12
                term = @(x,y) -(x-xi_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(3,4) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
                % c2
            elseif j == 13
                term = @(x,y) -y.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(1,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            elseif j == 14
                term = @(x,y) -x.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(2,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            elseif j == 15
                term = @(x,y) -(x.*(y-eta_p(i))+y.*(x-xi_p(i))).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
                tmp(3,5) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
            end
            
            infl_loc(:,:,i) = tmp;
        end
        
        send(D, i);
    catch
        i
        j
        le_eta
        te_eta
        xi_1(i)
        xi_3(i)
    end
    
%     waitbar(count/len, h)
%     count = count + 1;
    
end
close(h);

    function nUpdateWaitbar(~)
        waitbar(count/len, h);
        count = count + 1;
    end

% Only using the tangential velocities for points on the element
% infl_loc(1:2,:,idx_on_element) = infl_loc(1:2,:,idx_on_element).*0;
% Zero it all for singular edges
% infl_loc(:,:,idx_on_edge) = infl_loc(:,:,idx_on_edge).*0;

%%
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));
infl_tot(isnan(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);


end
