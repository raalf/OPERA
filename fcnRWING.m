function [vecR] = fcnRWING(strATYPE, valZTOL, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, ...
    matVLST, matDVE, matCOEFF, matPLEX, vecDVEFLIP, matWVLST, matWDVE, vecWDVEFLIP, vecHUB, vecWOFF, matWOFF, matWNORM, w_ind_in)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = sparse(valDLEN,1);

normals = matDVECT(matKINCON_DVE,:,3);
len = length(normals(:,1));

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK, normals,2);
elseif strcmpi(strATYPE{3}, 'UNSTEADY') && ~isempty(w_ind_in)
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK + w_ind_in, normals, 2);
else
    
    w_ind = (fcnSDVEVEL(matKINCON_P + matDVECT(matKINCON_DVE,:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, [], [], 0) + ...
        fcnSDVEVEL(matKINCON_P - matDVECT(matKINCON_DVE,:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, [], [], 0))./2;
    
    if any(~isnan(vecWOFF))
        w_ind_mirror = zeros(size(w_ind));
        for k = 1:length(vecWOFF)
            [tmpROTANG, tmpCENTER, tmpWROTANG, tmpWCENTER] = ...
                fcnMIRROR(matVLST, matDVE, vecDVEFLIP, matWVLST, matWDVE, vecWDVEFLIP, ...
                vecHUB, vecWOFF(k), matWOFF(k,:), matWNORM(k,:));
            
            w_ind_mirror = w_ind_mirror + ...
                fcnSDVEVEL(matKINCON_P, valNELE, matCOEFF, matPLEX, tmpROTANG, tmpCENTER, [], [], 0) + ...
                fcnSDVEVEL(matKINCON_P, valWNELE, matWCOEFF, matWPLEX, tmpWROTANG, tmpWCENTER, [], [], 0);
        end 
        w_ind = w_ind + w_ind_mirror;
    end
    
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF_KK + w_ind, normals, 2);
end

end

