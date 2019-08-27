function [R] = fcnROTMAX(u,theta)
%ROTATION Summary of this function goes here
%   Detailed explanation goes here
ux = u(1);
uy = u(2);
uz = u(3);

cost = cosd(theta);
sint = sind(theta);

R = [cost + ux.^2*(1 - cost) ux*uy*(1 - cost) - uz*sint ux*uz*(1 - cost) + uy*sint; ...
    uy*ux*(1 - cost) + uz*sint cost + uy.^2*(1 - cost) uy*uz*(1 - cost) - ux*sint; ...
    uz*ux*(1 - cost) - uy*sint uz*uy*(1 - cost) + ux*sint cost + uz.^2*(1 - cost); ...
    ];

end

