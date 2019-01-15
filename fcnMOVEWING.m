function [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P)
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

% 3D matrix of TE points, 3x?x2, rows are (x,y,z), columns are points, and depth is first and second TE point
te_flip(:,:,1) = matVLST(matELST(vecTE,1),:)';
te_flip(:,:,2) = matVLST(matELST(vecTE,2),:)';

idx_flip = [te_flip(2,:,1) > te_flip(2,:,2)]'; % Finding out which trailing edges go from out to in, so we know what to flip
idx_flip2 = logical([zeros(length(idx_flip),1); idx_flip]); % idx_flip1 is the first point, idx_flip2 is the second point

% NO-NO for RO-RO-ROTORSSSSSSSSS
translation = valDELTIME.*matUINF(1,:);

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

% Flipping the necessary vertices
temp = old_te(idx_flip,:);
old_te(idx_flip,:) = old_te(idx_flip2,:);
old_te(idx_flip2,:) = temp;

matVLST = matVLST - translation;
matCENTER = matCENTER - translation;
matCONTROL = matCONTROL - translation;
matKINCON_P = matKINCON_P - translation;

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

% Flipping the necessary vertices
temp = new_te(idx_flip,:);
new_te(idx_flip,:) = new_te(idx_flip2,:);
new_te(idx_flip2,:) = temp;

% Now that everything is flipped, the wake SHOULD generate with eta in the streamwise direction and
% the normals pointed upwards. If it doesn't, then this needs to be rewritten to ensure that is the case

% These vertices will be used to calculate the wake HDVE geometry
matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te((end/2)+1:end,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te(1:end/2,:)];




