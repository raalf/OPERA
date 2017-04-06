function [matVUINF, matCUINF] = fcnUINFROTOR(vecROTORIG, valROTRPM, matVLST, matCENTER)
% INPUT:
%   vecROTORIG - Origin of the rotation. Rotation happens in the z-direction.
%   valROTRPM - Rotation speed of the rotor (rev/s)
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
% OUTPUT:
%   matVUINF - Freestream values at each of the vertices
%   matCUINF - Freestream values at each of the in-centers

% UINF for vertex list
radius = matVLST - repmat(vecROTORIG, length(matVLST(:,1)),1);
radius = sum(abs(sqrt(radius(:,1:2).^2)),2);

matVUINF = 2.*pi.*valROTRPM.*radius;

% UINF for control point locations
radius = matCENTER - repmat(vecROTORIG, length(matCENTER(:,1)),1);
radius = sum(abs(sqrt(radius(:,1:2).^2)),2);

matCUINF = 2.*pi.*valROTRPM.*radius;

end

