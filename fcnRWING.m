function [vecR] = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, vecWDVESYM)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valDLEN,1);

normals = matDVECT(matKINCON_DVE,:,3);
len = length(normals(:,1));

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = -(4*pi).*dot(matUINF(matKINCON_DVE,:), normals,2);
else
    % WAKE INDUCED SHIT HERE
    w_ind = fcnSDVEVEL(matKINCON_P, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM,[]);
    vecR(end-(len-1):end) = -(4*pi).*dot(matUINF(matKINCON_DVE,:) + w_ind, normals, 2);
end

end

