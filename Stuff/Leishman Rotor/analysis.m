clc
clear

load('n10_165.mat')
cd ./../../

r_R = [];
circ = [];

omega = 2100.*(pi/30);
R = 0.4064;
c = 0.0425;

div = 10;
for i = 1:length(vecTEDVE)
    pts = matVLST(matELST(vecTE(i),:),:);
    
    pts_glob(:,:,i) = [linspace(pts(1,1), pts(2,1), div)', linspace(pts(1,2), pts(2,2), div)', linspace(pts(1,3), pts(2,3), div)'];
    pts_loc = fcnGLOBSTAR(pts_glob(:,:,i)-repmat(matCENTER(vecTEDVE(i),:),div,1), repmat(matROTANG(vecTEDVE(i),:),div,1));
    
    r_R = [r_R; sqrt(sum(pts_glob(:,1:2,i).^2,2))];
    
    circ = [circ; sum([0.5.*pts_loc(:,2).^2 pts_loc(:,2) 0.5.*pts_loc(:,1).^2 pts_loc(:,1) pts_loc(:,1).*pts_loc(:,2) ones(size(pts_loc(:,1)))].*matCOEFF(vecTEDVE(i),:),2)];
end

r_R = r_R./R;

hFig5 = figure(5);
clf(5);

plot(r_R, circ/(omega*R*c), '-k');
xlabel('r/R');
ylabel('\Gamma/\Omega Rc')
grid minor
axis tight

load('Stuff/Leishman Rotor/n20_123.mat')
r_R = [];
circ = [];

omega = 2100.*(pi/30);
R = 0.4064;
c = 0.0425;

div = 10;
for i = 1:length(vecTEDVE)
    pts = matVLST(matELST(vecTE(i),:),:);
    
    pts_glob(:,:,i) = [linspace(pts(1,1), pts(2,1), div)', linspace(pts(1,2), pts(2,2), div)', linspace(pts(1,3), pts(2,3), div)'];
    pts_loc = fcnGLOBSTAR(pts_glob(:,:,i)-repmat(matCENTER(vecTEDVE(i),:),div,1), repmat(matROTANG(vecTEDVE(i),:),div,1));
    
    r_R = [r_R; sqrt(sum(pts_glob(:,1:2,i).^2,2))];
    
    circ = [circ; sum([0.5.*pts_loc(:,2).^2 pts_loc(:,2) 0.5.*pts_loc(:,1).^2 pts_loc(:,1) pts_loc(:,1).*pts_loc(:,2) ones(size(pts_loc(:,1)))].*matCOEFF(vecTEDVE(i),:),2)];
end
r_R = r_R./R;
hold on
plot(r_R, circ/(omega*R*c), '--m');
hold off


load('Stuff/Leishman Rotor/single_rotor.mat');
hold on
scatter(single_rotor(:,1), single_rotor(:,2), 'sb');
hold off

legend('OPERA (n = 10)','OPERA (n = 20)','Bhagwat & Leishman','Location','NorthWest') 