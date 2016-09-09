clc
clear

% Input DVE number and field point in global, get back influence in global on that field point

load('quad.mat')

dvenum = 1;

fpg = [0 0 0.5];

fpl = fcnTOLOC(dvenum, fpg, DVE, DVECT, VLST, DNORM);
k = [1 1]';

% Left to Right

endpoints(1,:,1) = PLEX(1,:,1);
endpoints(1,:,2) = PLEX(2,:,1);

endpoints(2,:,1) = PLEX(2,:,1);
endpoints(2,:,2) = PLEX(3,:,1);

phi = [atan(PLEX(2,1,1)./PLEX(2,2,1)); atan((PLEX(3,1,1)-PLEX(2,1,1))./PLEX(2,2,1))];

psi = [90 90]';

[al, bl, cl] = fcnVSIND(endpoints, phi, deg2rad(psi), [fpl; fpl], k);

a3 = fcnROTVECT(1, al(1,:)-al(2,:),DVECT);
a2 = fcnROTVECT(1, bl(1,:)-bl(2,:),DVECT);
a1 = fcnROTVECT(1, cl(1,:)-cl(2,:),DVECT);


% Up to Down

endpoints(1,:,1) = PLEX(1,:,1);
endpoints(1,:,2) = PLEX(2,:,1);

endpoints(2,:,1) = PLEX(2,:,1);
endpoints(2,:,2) = PLEX(3,:,1);

endpoints(3,:,1) = PLEX(1,:,1);
endpoints(3,:,2) = PLEX(3,:,1);

phi = [atan(PLEX(2,2,1)./PLEX(2,1,1)); atan((PLEX(2,2,1)./PLEX(3,1,1)-PLEX(2,1,1))); 0];

psi = [0 0 0]';
k = [1 1 1]';

[al, bl, cl] = fcnVSIND(endpoints, phi, deg2rad(psi), [fpl; fpl; fpl], k);

b3 = fcnROTVECT(1, al(1,:)+al(2,:)-al(3,:),DVECT);
b2 = fcnROTVECT(1, bl(1,:)+bl(2,:)-bl(3,:),DVECT);
b1 = fcnROTVECT(1, cl(1,:)+cl(2,:)-cl(3,:),DVECT);

D = [a1' a2' b1' b2' (a3+b3)'];

COEFF = [0 1 0 1 0]';

q = D*COEFF;

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT)
hold on
% scatter3(fpg(1), fpg(2), fpg(3),100,'ok','filled');
quiver3(fpg(1), fpg(2), fpg(3), q(1), q(2), q(3));
hold off



