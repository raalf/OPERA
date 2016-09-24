function [matVLST, matCENTER, matVUINF, matCUINF, matNEWWAKE] = fcnMOVEROTOR(vecROTORIG, valDELROT, valROTRPM, matVLST, matCENTER, matELST, vecTE)

% This function rotates a rotor about the z-axis by valDELROT angle at valROTRPM timestep size.
% INPUT:
%   vecROTORIG - Origin of the rotation. Rotation happens in the z-direction.
%   valDELROT - Angle by which to rotate the rotor (radians)
%   valROTRPM - Rotation speed of the rotor (rev/s)
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
%   matELST - List of edges from fcnTRIANG
%   vecTE - List of trailing edge edges
% OUTPUT:
%   matVLST - Updated vertex list
%   matCENTER - Updated in-center list
%   matVUINF - Freestream values at each of the vertices
%   matCUINF - Freestream values at each of the in-centers
%   matNEWWAKE - Outputs a n x 3 x 3 matrix of points for the wake triangulation

% UINF for vertex list
radius = matVLST - repmat(vecROTORIG, length(matVLST(:,1)),1);
radius = sum(abs(sqrt(radius(:,1:2).^2)),2);

matVUINF = 2.*pi.*valROTRPM.*radius;

% UINF for control point locations
radius = matCENTER - repmat(vecROTORIG, length(matCENTER(:,1)),1);
radius = sum(abs(sqrt(radius(:,1:2).^2)),2);

matCUINF = 2.*pi.*valROTRPM.*radius;

% Rotation
% 1. Translate space so that the rotation axis passes through the origin
% 2. Perform the desired rotation by theta about the z axis
% 3. Apply the inverse of step (1)

tempVLST = matVLST - repmat(vecROTORIG, length(matVLST(:,1)),1);
tempCENTER = matCENTER - repmat(vecROTORIG, length(matCENTER(:,1)),1);

ROT = [cos(valDELROT) -sin(valDELROT) 0; sin(valDELROT) cos(valDELROT) 0; 0 0 1];

vlst2 = (ROT*tempVLST')' + repmat(vecROTORIG, length(matVLST(:,1)),1);
center2 = (ROT*tempCENTER')' + repmat(vecROTORIG, length(matCENTER(:,1)),1);

% Old trailing edge of wing
old_te = matVLST(matELST(vecTE,:),:);

% New trailing edge of wing
new_te = vlst2(matELST(vecTE,:),:);

% These vertices will be used to calculate the wake HDVE geometry
matNEWWAKE(:,:,1) = [new_te(1:end/2,:); old_te((end/2)+1:end,:)];
matNEWWAKE(:,:,2) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,3) = [old_te(1:end/2,:); new_te((end/2)+1:end,:)];

% Replacing old values with new locations
matVLST = vlst2;
matCENTER = center2;


end

