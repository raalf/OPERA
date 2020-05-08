function mirrored_points = fcnMIRRORPTS(points, vecHUB, valWOFF, vecWOFF, vecWNORM)
plane_orig = vecHUB + valWOFF.*vecWOFF;
d = -vecWNORM(1).*plane_orig(1) - vecWNORM(2).*plane_orig(2) - vecWNORM(3).*plane_orig(3);
D = (vecWNORM(1).*points(:,1) + vecWNORM(2).*points(:,2) + vecWNORM(3).*points(:,3) + d);
mirrored_points = points - 2.*D.*vecWNORM;
end