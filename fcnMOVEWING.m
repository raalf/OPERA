function [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matELST, vecTE)
% This function moves a wing (NOT rotor) by translating all of the vertices
% in the VLST and the in-centers of each triangle in CENTER.

% INPUT:
%   valALPHA - Angle of attack for this move (radians)
%   valBETA - Sideslip angle for this move (radians)
%   valDELTIME - Timestep size for this move
%   matVLST - Vertex list from fcnTRIANG
%   matCENTER - In-center list from fcnTRIANG
%   matELST - List of edges from fcnTRIANG
%   vecTE - List of trailing edge edges
% OUTPUT:
%   matVLST - New vertex list with moved points
%   matCENTER - New in-center list with moved points
%   matVUINF - Freestream values at each of the vertices
%   matCUINF - Freestream values at each of the in-centers
%   matNEWWAKE - Outputs a n x 3 x 3 matrix of points for the wake triangulation

% In order to ensure the local coordinate system of the wake DVEs is aligned
% with the global streamwise and spanwise directions, we need to know which trailing
% edges to "flip," if they naturally go from outwards to inwards. If we don't do
% this flip, then the normals of those elements will be facing downards

te_flip = permute(reshape(matVLST(matELST(vecTE,:),:)',[],3,2),[2 1 3]); % Just creating a 3d matrix of wake vertices
idx_flip = te_flip(:,2,1) > te_flip(:,2,2); % Finding out which trailing edges go from out to in, so we know what to flip
idx_flip2 = logical([zeros(length(idx_flip),1); idx_flip]); % idx_flip1 is the first point, idx_flip2 is the second point

% Won't work for rotors
uinf = 1;
uinf = [uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)];

translation = valDELTIME.*uinf;

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

% Flipping the necessary vertices
temp = old_te(idx_flip,:);
old_te(idx_flip,:) = old_te(idx_flip2,:);
old_te(idx_flip2,:) = temp;

matVLST = matVLST - repmat(translation, length(matVLST(:,1)), 1);
matCENTER = matCENTER - repmat(translation, length(matCENTER(:,1)), 1);

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

% Flipping the necessary vertices
temp = new_te(idx_flip,:);
new_te(idx_flip,:) = new_te(idx_flip2,:);
new_te(idx_flip2,:) = temp;

% Now that everything is flipped, the wake SHOULD generate with eta in the streamwise direction and
% the normals pointed upwards. If it doesn't, then this needs to be rewritten to ensure that is the case

% These vertices will be used to calculate the wake HDVE geometry
matNEWWAKE(:,:,1) = [new_te(1:end/2,:); old_te((end/2)+1:end,:)];
matNEWWAKE(:,:,2) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,3) = [old_te(1:end/2,:); new_te((end/2)+1:end,:)];