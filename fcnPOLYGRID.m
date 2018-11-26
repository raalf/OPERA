function [inPoints] = fcnPOLYGRID( xv, yv, N)

%Find the bounding rectangle
lower_x = min(xv);
higher_x = max(xv);

lower_y = min(yv);
higher_y = max(yv);
%Create a grid of points within the bounding rectangle

inc_x = (higher_x - lower_x)/N;
inc_y = (higher_y - lower_y)/N;


interval_x = lower_x:inc_x:higher_x;
interval_y = lower_y:inc_y:higher_y;
[bigGridX, bigGridY] = meshgrid(interval_x, interval_y);

%Filter grid to get only points in polygon
[in,on] = inpolygon(bigGridX(:), bigGridY(:), xv, yv);
in = in | on;

%Return the co-ordinates of the points that are in the polygon
inPoints = [bigGridX(in), bigGridY(in)];

end


