function [vecR] = fcnRWING(valDLEN, valTIMESTEP, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valDLEN,1);

points = matCENTER;
normals = matDVECT(:,:,3);

len = length(normals(:,1));

if valTIMESTEP < 1
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF, normals,2);
else
    % WAKE INDUCED SHIT HERE
    w_ind = fcnWINDVEL(points, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    vecR(end-(len-1):end) = (4*pi).*dot(matUINF + w_ind, normals, 2);
end

end

