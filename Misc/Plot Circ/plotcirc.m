clc
clear

load('quad.mat')

test_num = 1;

lambda = [1 0 0];

x1 = PLEX(1,1,test_num);
x2 = PLEX(2,1,test_num);
x3 = PLEX(3,1,test_num);

y1 = PLEX(1,2,test_num);
y2 = PLEX(2,2,test_num);
y3 = PLEX(3,2,test_num);

z1 = PLEX(1,3,test_num);
z2 = PLEX(2,3,test_num);
z3 = PLEX(3,3,test_num);

count = 1;
for i = 0:0.05:1
   for j = 1-i:-0.05:0
       lambda(count,1) = i;
       lambda(count,2) = j;
       lambda(count,3) = 1-i-j;
       
       count = count + 1;
   end
end

eta = (x1.*lambda(:,1) + x2.*lambda(:,2) + x3.*lambda(:,3));
xsi = (y1.*lambda(:,1) + y2.*lambda(:,2) + y3.*lambda(:,3));

circ = matCOEFF(test_num,1).*(eta.^2) + matCOEFF(test_num,2).*eta + matCOEFF(test_num,3).*(xsi.^2) ...
    + matCOEFF(test_num,4).*xsi + matCOEFF(test_num,5);

vort = matCOEFF(test_num,1).*eta + matCOEFF(test_num,2) + matCOEFF(test_num,3).*xsi + matCOEFF(test_num,4);

hFig10 = figure(10);
clf(10);
patch(PLEX(:,1,test_num),PLEX(:,2,test_num),'b')
alpha(0.5);
xlabel('eta-direction','FontSize',15);
ylabel('xsi-direction','FontSize',15);
axis equal
grid on
box on
hold on
% surf(eta, xsi, circ)
% scatter3(eta, xsi, vort)
DT = delaunay(eta, xsi)
trisurf(DT, eta, xsi, circ,'EdgeColor','k','FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0.5)







