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

valJ = reshape(valJ, 1, 1, []);
vecUINF = [cosd(valALPHA)*cosd(0) sind(0) sind(valALPHA)*cosd(0)];
translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*vecUINF;

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

    uinf = cross(repmat([0, 0, -vecROTORRADPS], length(r_R),1), locs) - translation;

    tmp = dot(uinf, repmat(forward, size(uinf,1), 1, 2), 2);
    reversal = ones(size(uinf,1),1,2);
    reversal(tmp > 0) = -1;

    z(:,i,:) = sqrt(sum(uinf.^2,2)).*reversal;
    x(:,i) = locs(:,1);
    y(:,i) = locs(:,2);
end

% hFig2 = figure(2);
% clf(2);
% 
% z2 = mean(z,1);
% z3 = z2(:,azs < 180,:) + z2(:,azs >= 180,:);
% azs2 = azs(azs < 180);
% scatter(azs2, z3(:,:,1), 'ok');
% hold on
% scatter(azs2, z3(:,:,2), '^b');
% hold off
% grid minor
% box on
% axis tight


hFig1 = figure(1);
clf(1);
h1 = subplot(1,2,1);
contourf(x,y,z(:,:,1));
title(['\mu = ', num2str(valJ(1))])
axis equal
colormap(jet)
colorbar;
c1 = caxis;
axis off

h2 = subplot(1,2,2);
contourf(x,y,z(:,:,2));
title(['\mu = ', num2str(valJ(2))])
axis equal
colormap(jet)
colorbar;
c2 = caxis;
axis off
lims = [min([c1 c2]) max([c1 c2])];
caxis(lims);
set(hFig1, 'currentaxes', h1);
caxis(lims);

sgtitle(['Percent V/V_{mean} vs. Azimuth Location (\alpha = 0, \DeltaT = ', num2str(valDELTIME), ', RPM = ', num2str(valRPM), ', VINF RIGHT TO LEFT, CCW)'])

