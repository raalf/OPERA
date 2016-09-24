function [matVUINF, matCUINF] = fcnUINFWING(valALPHA, valBETA, matVLST, matCENTER)

% INPUT:
%   valALPHA - Angle of attack for this move (radians)
%   valBETA - Sideslip angle for this move (radians)
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
% OUTPUT:
%   matVUINF - Freestream values at each of the vertices
%   matCUINF - Freestream values at each of the in-centers

uinf = 1;

uinf = [uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)];

% Saving velocities at all points (uniform for every vertex with a wing)
matVUINF = repmat(uinf, length(matVLST(:,1)), 1);
matCUINF = repmat(uinf, length(matCENTER(:,1)), 1);

end

