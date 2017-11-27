function [vecR] = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matELST, matCENTER, matDVECT, vecUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST, matWVSCOMB, matWCENTER)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valDLEN,1);

points = matCENTER;
normals = matDVECT(:,:,3);

% points = matVLST;
% normals = matVNORM;
% 
% points = [matCENTER; (matVLST(matELST(vecLE,1),:) + matVLST(matELST(vecLE,2),:))./2];
% normals = [matDVECT(:,:,3); matDVECT(vecLEDVE,:,3)];

len = length(normals(:,1));

if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1), normals,2);
else
    [w_ind] = fcnWINDVEL(points, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matWVSCOMB, matWCENTER);
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1) + w_ind, normals, 2);
end

end

