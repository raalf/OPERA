function [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P, vecWDVEFLIP] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P, vecWDVEFLIP)
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

% NO-NO for RO-RO-ROTORSSSSSSSSS
translation = valDELTIME.*matUINF(1,:);

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

matVLST = matVLST - translation;
matCENTER = matCENTER - translation;
matCONTROL = matCONTROL - translation;
matKINCON_P = matKINCON_P - translation;

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

% These vertices will be used to calculate the wake HDVE geometry
matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te((end/2)+1:end,:)];

vecWDVEFLIP = [vecWDVEFLIP; false(size(matNEWWAKE,1)./2, 1); true(size(matNEWWAKE,1)./2, 1)]; 

end



