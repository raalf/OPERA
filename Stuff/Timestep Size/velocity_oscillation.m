clc
clear

valDIAM = 0.4572;
valRPM = 3000;
valJ = [0.15 0.3];
valALPHA = 0;
valDELTIME = 0.0005;

r_R = linspace(0.2, 1, 20)';

vecHUB = [0 0 0];
vecROTORRADPS = valRPM.*2.*pi./60;

vecUINF = [cosd(valALPHA)*cosd(0) sind(0) sind(valALPHA)*cosd(0)];
translation(:,:,1) = valJ(1).*(valRPM.*(pi/30)).*(valDIAM/2).*vecUINF;
translation(:,:,2) = valJ(2).*(valRPM.*(pi/30)).*(valDIAM/2).*vecUINF;

load('chord.mat')
c = interp1(chord(:,1), chord(:,2), r_R).*(valDIAM/2);

locs = [r_R.*0 r_R.*(valDIAM/2) r_R.*0];
fordir = [-1 0 0];
locs_og = locs;
azs = linspace(0, 360, 36)';

for i = 1:length(azs)
    dcmROTORSTEP = angle2dcm(deg2rad(azs(i)),0,0,'ZXY');
    locs = locs_og*dcmROTORSTEP;
    forward = fordir*dcmROTORSTEP;

    uinf(:,:,1) = cross(repmat([0, 0, -vecROTORRADPS], length(r_R),1), locs) - translation(:,:,1);
    uinf(:,:,2) = cross(repmat([0, 0, -vecROTORRADPS], length(r_R),1), locs) - translation(:,:,2);
    tmp = dot(uinf, repmat(forward, size(uinf,1), 1, 2), 2);
    reversal = ones(size(uinf,1),1,2);
    reversal(tmp > 0) = -1;

    z(:,i,1) = sqrt(sum(uinf(:,:,1).^2,2)).*reversal(:,:,1);
    z(:,i,2) = sqrt(sum(uinf(:,:,2).^2,2)).*reversal(:,:,2);
    x(:,i) = locs(:,1);
    y(:,i) = locs(:,2);
end

amplitude(:,:,1) = valJ(1)./(r_R.*(valDIAM/2));
amplitude(:,:,2) = valJ(2)./(r_R.*(valDIAM/2));
amplitude = repmat(amplitude, 1, length(azs), 1);
% amplitude = 100.*abs((z - repmat(mean(z,1), size(z,1), 1, 1))./repmat(mean(z,1), size(z,1), 1, 1));

hFig1 = figure(1);
clf(1);
subplot(1,2,1)
contourf(x,y,amplitude(:,:,1));
title(['\mu = ', num2str(valJ(1))])
axis equal
colormap(jet)
colorbar;
axis off

subplot(1,2,2)
contourf(x,y,amplitude(:,:,2));
title(['\mu = ', num2str(valJ(2))])
axis equal
colormap(jet)
colorbar;
axis off

sgtitle(['Percent V/V_{mean} vs. Azimuth Location (\alpha = 0, \DeltaT = ', num2str(valDELTIME), ', RPM = ', num2str(valRPM), ', VINF RIGHT TO LEFT, CCW)'])

