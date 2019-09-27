clc
clear

%% Getting tunnel data
A = dlmread('Alpha_15_550.txt', '', 9, 0);
A(A(:,4) > 360,:) = [];

valDIAM = 0.4572;
valRPM = 3000;
valJ = mean(A(:,3))/((valRPM.*(pi/30)).*(valDIAM/2));

vecPOS_TUNNEL_OG = A(:,4);
vecPOS_TUNNEL = A(:,4);

idx_go = vecPOS_TUNNEL(2:end) < vecPOS_TUNNEL(1:end-1);
offset = [0 360:360:(length(find(idx_go))*360)]';
offset_idx = 1;
count = 1;
for i = 1:length(idx_go)
    vecPOS_TUNNEL(i) = vecPOS_TUNNEL(i) + offset(offset_idx);
    if idx_go(i) == true
        offset_idx = offset_idx + 1;
    end
end

vecDENSITY = (3386.39*29.23)./(287.058.*((A(:,2)-32).*(5/9) + 273.15));
CT_tunnel = A(:,7)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
error_tunnel = (1/4)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));


%%

hFig99 = figure(99);
clf(99);

[vecPOS_TUNNEL_OG, idx] = sort(vecPOS_TUNNEL_OG, 'ascend');
CT_tunnel = CT_tunnel(idx);
plot(vecPOS_TUNNEL_OG, smooth(CT_tunnel), 'sk')
grid minor
box on
axis tight

%% Getting tunnel data
A = dlmread('Alpha_5_550.txt', '', 9, 0);
A(A(:,4) > 360,:) = [];

valDIAM = 0.4572;
valRPM = 3000;
valJ = mean(A(:,3))/((valRPM.*(pi/30)).*(valDIAM/2));

vecPOS_TUNNEL_OG = A(:,4);
vecPOS_TUNNEL = A(:,4);

idx_go = vecPOS_TUNNEL(2:end) < vecPOS_TUNNEL(1:end-1);
offset = [0 360:360:(length(find(idx_go))*360)]';
offset_idx = 1;
count = 1;
for i = 1:length(idx_go)
    vecPOS_TUNNEL(i) = vecPOS_TUNNEL(i) + offset(offset_idx);
    if idx_go(i) == true
        offset_idx = offset_idx + 1;
    end
end

vecDENSITY = (3386.39*29.23)./(287.058.*((A(:,2)-32).*(5/9) + 273.15));
CT_tunnel = A(:,7)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
error_tunnel = (1/4)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));


%%

[vecPOS_TUNNEL_OG, idx] = sort(vecPOS_TUNNEL_OG, 'ascend');
CT_tunnel = CT_tunnel(idx);
hold on
plot(vecPOS_TUNNEL_OG, smooth(CT_tunnel), 'dr')
hold off

%% Getting tunnel data
A = dlmread('Alpha_0_550.txt', '', 9, 0);
A(A(:,4) > 360,:) = [];

valDIAM = 0.4572;
valRPM = 3000;
valJ = mean(A(:,3))/((valRPM.*(pi/30)).*(valDIAM/2));

vecPOS_TUNNEL_OG = A(:,4);
vecPOS_TUNNEL = A(:,4);

idx_go = vecPOS_TUNNEL(2:end) < vecPOS_TUNNEL(1:end-1);
offset = [0 360:360:(length(find(idx_go))*360)]';
offset_idx = 1;
count = 1;
for i = 1:length(idx_go)
    vecPOS_TUNNEL(i) = vecPOS_TUNNEL(i) + offset(offset_idx);
    if idx_go(i) == true
        offset_idx = offset_idx + 1;
    end
end

vecDENSITY = (3386.39*29.23)./(287.058.*((A(:,2)-32).*(5/9) + 273.15));
CT_tunnel = A(:,7)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
error_tunnel = (1/4)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));


%%

[vecPOS_TUNNEL_OG, idx] = sort(vecPOS_TUNNEL_OG, 'ascend');
CT_tunnel = CT_tunnel(idx);
hold on
plot(vecPOS_TUNNEL_OG, smooth(CT_tunnel), 'ob')
hold off

legend('Alpha 15','Alpha 5','Alpha 0','Location','SouthWest')