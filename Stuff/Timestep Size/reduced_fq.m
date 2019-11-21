clc
clear

valDIAM = 0.4572;
valRPM = 3000;
valJ = 0.3;
valALPHA = 0;
valDELTIME = [0.00025];

r_R = linspace(0.2, 1, 20)';

vecHUB = [0 0 0];
vecROTORRADPS = valRPM.*2.*pi./60;

vecUINF = [cosd(valALPHA)*cosd(0) sind(0) sind(valALPHA)*cosd(0)];
translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*vecUINF;

load('chord.mat')
c = interp1(chord(:,1), chord(:,2), r_R).*(valDIAM/2);
omega = valRPM.*(pi/30);

locs = [r_R.*0 r_R.*(valDIAM/2) r_R.*0];
fordir = [-1 0 0];
locs_og = locs;
azs = linspace(0, 360, 36)';

for i = 1:length(azs)
    dcmROTORSTEP = angle2dcm(deg2rad(azs(i)),0,0,'ZXY');
    locs = locs_og*dcmROTORSTEP;
    forward = fordir*dcmROTORSTEP;

    uinf = cross(repmat([0, 0, -vecROTORRADPS], length(r_R),1), locs) - translation;
    tmp = dot(uinf, repmat(forward, size(uinf,1),1), 2);
    reversal = ones(size(uinf,1),1);
    reversal(tmp > 0) = -1;
%     flow_reversed = 
    z(:,i,1) = (0.5.*c.*omega)./(sqrt(sum(uinf.^2,2)).*reversal);
    z(:,i,2) = (0.5.*c.*omega)./(sqrt(sum(uinf.^2,2)).*reversal);
    z(:,i,3) = (0.5.*c.*omega)./(sqrt(sum(uinf.^2,2)).*reversal);
    z(:,i,4) = (0.5.*c.*omega)./(sqrt(sum(uinf.^2,2)).*reversal);
    x(:,i) = locs(:,1);
    y(:,i) = locs(:,2);
end

hFig1 = figure(1);
clf(1);
% subplot(2,2,1)
contourf(x,y,z(:,:,1));
% title(['\DeltaT = ', num2str(valDELTIME(1))])
axis equal
colormap(jet)
colorbar;
axis off

% subplot(2,2,2)
% contourf(x,y,z(:,:,2));
% title(['\DeltaT = ', num2str(valDELTIME(2))])
% axis equal
% colormap(jet)
% colorbar;
% axis off
% 
% subplot(2,2,3)
% contourf(x,y,z(:,:,3));
% title(['\DeltaT = ', num2str(valDELTIME(3))])
% axis equal
% colormap(jet)
% colorbar;
% axis off
% 
% subplot(2,2,4)
% contourf(x,y,z(:,:,4));
% title(['\DeltaT = ', num2str(valDELTIME(4))])
% axis equal
% colormap(jet)
% colorbar;
% axis off

sgtitle(['Reduced Frequency vs. Azimuth Location (\alpha = 0, \mu = ', num2str(valJ), ', RPM = ', num2str(valRPM), ', VINF RIGHT TO LEFT, CCW)'])

