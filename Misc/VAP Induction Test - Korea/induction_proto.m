% clc
clear

% Vertices of the HDVE num
% testDVE = VLST(DVE(num,:,1),:);
testDVE = [0.5 -0.5 0; 0.5 0.5 0; -0.5 0.5 0; -0.5 -0.5 0];
DNORM = [0 0 1];
DVECT(:,:,1) = [1 0 0];
DVECT(:,:,2) = [0 1 0];
DVECT(:,:,3) = [0 0 1];

% FP = [0 -2 0];
% FP = [0 3 0];
% FP = [0.55 -3 0.56]
% FP = [0.5 -1 0.5];
% FP = [0.5 3 0.5];
% FP = [0.5 1 0.5];
% FP = [0 5 5];
% FP = [0 -5 5];
FP = [0 0 .5];

COEFF = [1 1 1 1 0]';

D = INFLUENCECOEFF2(testDVE, DVECT, FP);

q = D*COEFF

%%
% y_p = -3;
% z_p = -3;
count = 1;
granularity = 0.5;
x_p = 0;
y_p = 0;
z_p = 0;
for y_p = -5:granularity:5
%     for x_p = -5:granularity:5
        for z_p = -5:granularity:5
            
            FP = [x_p y_p z_p];
            
            D = INFLUENCECOEFF2(testDVE, DVECT, FP);
            
            D*COEFF
            
            velo(count,1:3) = [x_p, y_p, z_p];
            velo(count,4:6) = transpose(D*COEFF);
            
            count = count + 1;
            
        end
%     end
end

% x_p = 0;
% y_p = 0;
% z_p = 0;
% for x_p = -5:granularity:5
% %     for x_p = -5:granularity:5
%         for z_p = -5:granularity:5
%             
%             FP = [x_p y_p z_p];
%             
%             D = INFLUENCECOEFF2(testDVE, DVECT, FP);
%             
%             D*COEFF
%             
%             velo(count,1:3) = [x_p, y_p, z_p];
%             velo(count,4:6) = transpose(D*COEFF);
%             
%             count = count + 1;
%             
%         end
% %     end
% end

hFig1 = figure(1);
clf(1)
patch(testDVE(:,1),testDVE(:,2),testDVE(:,3),'r','LineWidth',2)
alpha(0.5);
hold on
% scatter3(FP(1),FP(2),FP(3),100,'ok','filled');
quiver3(velo(:,1), velo(:,2), velo(:,3), velo(:,4), velo(:,5), velo(:,6),0);
hold off
grid on
axis equal
axis tight
box on
xlabel('X-dir','FontSize',15);
ylabel('Y-dir','FontSize',15);
zlabel('Z-dir','FontSize',15);

