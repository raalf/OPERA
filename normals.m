clc
clear

p0 = [0 0 0];
p1 = [1 1 0];
p2 = [1 0.8 0];

vectarrow(p0,p1)
hold on
vectarrow(p0,p2)
hold off
box on
grid minor
axis tight
axis equal
view(2)

a = (p1-p0)./norm(p1-p0);
b = (p2-p0)./norm(p2-p0);

angle = (atan2d(b(2),b(1)) - atan2d(a(2),a(1)))
angle(angle < 0) = angle + rad2deg(2*pi)