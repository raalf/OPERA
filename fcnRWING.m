function [vecR] = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valDLEN,1);

len = length(matCENTER(:,1));

if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1), matDVECT(:,:,3),2);
    
%     % Trailing edge flow tangency goes just above the previous stuff in the resultant
% 
%     len2 = length(vecTE);
    
%     % If it is a panel code, we use the average of the norm of the 2 panels at the trailing
%     % edge. 
%     if strcmp(strATYPE,'PC') == 1
%         normals(:,:,1) = matDVECT(matEATT(vecTE,1),:,3);
%         normals(:,:,2) = matDVECT(matEATT(vecTE,2),:,3);
%         normals = mean(normals,3);
%     else
%         normals = matDVECT(nonzeros(matEATT(vecTE,:)),:,3);
%     end
%     
%     vecR(end-len-(len2-1):end-len) = (4*pi).*dot(repmat(vecUINF,len2,1), normals,2);
    
else
    [w_ind] = fcnWINDVEL(matCENTER, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX);

    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1) + w_ind, matDVECT(:,:,3), 2);
end

end

