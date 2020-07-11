function [vecR] = fcnRWING(flgSTEADY, valZTOL, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, matKINCON_P, vecKINCON_DVE, matDVECT, w_ind_in)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = sparse(valDLEN,1);

normals = matDVECT(vecKINCON_DVE,:,3);
len = length(normals(:,1));

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK, normals,2);
elseif flgSTEADY == false && ~isempty(w_ind_in)
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK + w_ind_in, normals, 2);
else
    w_ind = (fcnSDVEVEL(matKINCON_P + matDVECT(vecKINCON_DVE,:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, [], 0) + ...
        fcnSDVEVEL(matKINCON_P - matDVECT(vecKINCON_DVE,:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, [], 0))./2;
       
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK + w_ind, normals, 2);
end

end

