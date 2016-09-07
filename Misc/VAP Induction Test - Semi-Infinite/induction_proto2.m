clc
clear

endpoints(:,:,1) = [-0.5 -0.5 0; 0.5 -0.5 0];
endpoints(:,:,2) = [-0.5 0.5 0; 0.5 0.5 0];

phi = [0; 0];

eta_VS = [1 0 0; 1 0 0];

% fpl = repmat([1 5 5],2,1);

k = 0;

% endpoints(:,:,1) = [-0.5 -0.5 0];
% endpoints(:,:,2) = [-0.5 0.5 0];
%
% phi = [0];
%
% eta_VS = [1 0 0];
%
% fpl = repmat([1 5 5],1,1);
%
% k = 0;

COEFF = [0 1 0]';

gran = 0.5;
x_fp = 0;
y_fp = 0;
z_fp = 1;
count = 1;

for z_fp = 0.1:gran:3
    for y_fp = -3:gran:3
        fpl = [x_fp y_fp z_fp];
        [al, bl, cl] = fcnVSIND(endpoints, phi, eta_VS, fpl, k);
        a = al(1,:) - al(2,:);
        b = bl(1,:) - bl(2,:);
        c = cl(1,:) - cl(2,:);
        D = [a' b' c'];
        q = D*COEFF;
        velo(count,:) = [x_fp y_fp z_fp q(1) q(2) q(3)];
        count = count + 1;
    end
end

hFig10 = figure(10);
clf(10);
patch([-0.5 0.5 0.5 -0.5],[-0.5 -0.5 0.5 0.5],[0 0 0 0],'m')
alpha(0.5);
xlabel('X-dir','FontSize',15);
ylabel('Y-dir','FontSize',15);
zlabel('Z-dir','FontSize',15);
grid on
box on
axis equal
hold on
quiver3(velo(:,1), velo(:,2), velo(:,3), velo(:,4), velo(:,5), velo(:,6), 0)
hold off


%%
load('quad.mat');

% % input two end points in local, phi, eta_VS direction and field point in local, singularity factor?

dvenum = 1;

endpoints(:,:,1) = PLEX(2,:,dvenum);
endpoints(:,:,2) = PLEX(3,:,dvenum);

fpg = [1 1 0.5];

fpl = fcnTOLOC(dvenum, fpg, DVE, DVECT, VLST, DNORM);

% eta_VS direction
eta_VS = [0 1 0];

% midpoint = mean(endpoints,3)
% r1 = endpoints(:,:,2) - midpoint
% r2 = midpoint + eta_VS.*endpoints(:,:,2)
% % leading edge vector of vortex sheet
R1 = (endpoints(:,:,2) - endpoints(:,:,1));
% unit vector along leading edge of vortex sheet
r1 = R1./sqrt(sum(abs(R1).^2,2));

% if we want eta direction sheet, then
dot(r1, eta_VS,2)
phi = acos(dot(r1,eta_VS,2)); % FIX THISSSSSSSSSSSSSSSSSSSSSSS

k = 0;



num = length(al(:,1));
temp = [al; bl; cl];
dvenum = zeros(length(temp(:,1)),1) + dvenum;

% this probably isn't right, cause a, b and c are in terms of the VS local ref frame, not the DVE
% and even if that did work, this function may only work for points
% temp2 = fcnTOGLOB(dvenum, temp, DVE, DVECT, VLST);
temp2 = temp;

a = temp2(1:num,:);
b = temp2(num+1:2*num,:);
c = temp2(2*num+1:3*num,:);

D = [a' b' c'];
COEFF = [1 1 1]';

q_ind = D*COEFF

%%
% % Plotting global and local to visualize
test_num = 1;
hFig2 = figure(2);
clf(2);
subplot(2,1,1)
patch(PLEX(:,1,test_num),PLEX(:,2,test_num),'b')
hold on
scatter3(fpl(1), fpl(2), fpl(3),100,'ok','filled')
hold off
alpha(0.5);
xlabel('eta-direction','FontSize',15);
ylabel('xi-direction','FontSize',15);
axis equal
grid on
box on
subplot(2,1,2)
patch(VLST(DVE(test_num,:,1),1),VLST(DVE(test_num,:,1),2),VLST(DVE(test_num,:,1),3),'r')
hold on
scatter3(fpg(1), fpg(2), fpg(3),100,'ok','filled')
verbose = 1;
if verbose == 1
    for ii = 1:NELE
        str = sprintf('%d',ii);
        text(DVE(ii,1,3),DVE(ii,2,3),DVE(ii,3,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(VLST(:,1))
        str = sprintf('%d',ii);
        text(VLST(ii,1),VLST(ii,2),VLST(ii,3),str,'Color','g','FontSize',20);
    end
    
    edge1 = VLST(ELST(:,1),:);
    edge2 = VLST(ELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end
    
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,1), DVECT(:,2,1), DVECT(:,3,1), 0.25, 'b') % eta
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,2), DVECT(:,2,2), DVECT(:,3,2), 0.25, 'k') % xi
    quiver3(DVE(:,1,3), DVE(:,2,3), DVE(:,3,3), DVECT(:,1,3), DVECT(:,2,3), DVECT(:,3,3), 0.25,'m') % zeta (normal)
    
end

hold off
alpha(0.5);
xlabel('X-direction','FontSize',15);
ylabel('Y-direction','FontSize',15);
zlabel('Z-direction','FontSize',15);
axis equal
grid on
box on




































