function [vecR] = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST, matWVSCOMB, matWCENTER)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valDLEN,1);

points = matCENTER;
normals = matDVECT(:,:,3);
% points = matVLST;
% normals = matVNORM;

len = length(normals(:,1));

if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF, normals,2);
%     vecR(end-(len-1):end) = (-4*pi).*dot(repmat(matUINF(1,:),size(normals,1),1), normals,2);
else
    [w_ind] = fcnWINDVEL(points, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matWVSCOMB, matWCENTER);
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF + w_ind, normals, 2);
end

end

