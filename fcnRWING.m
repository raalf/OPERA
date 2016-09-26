function [R] = fcnRWING(ATYPE, valDLEN, valTIMESTEP, VLST, ELST, EATT, CENTER, DVECT, matVUINF, matCUINF, vecTE)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

R = zeros(valDLEN,1);


if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    len = length(CENTER(:,1));
    R(end-(len-1):end) = (4*pi).*dot(matCUINF, DVECT(:,:,3),2);
    
    % Trailing edge flow tangency goes just above the previous stuff in the resultant
    
    % Flow tangency was enforced at the midpoint of the trailing edge edges, so the UINF
    % at that point is the average of the two vertices
    vuinf(:,:,1) = matVUINF(ELST(vecTE,1),:);
    vuinf(:,:,2) = matVUINF(ELST(vecTE,2),:);
    vuinf = mean(vuinf,3);
    
    len2 = length(vuinf(:,1));
    
    % If it is a panel code, we use the average of the norm of the 2 panels at the trailing
    % edge. 
    if strcmp(ATYPE,'PC') == 1
        normals(:,:,1) = DVECT(EATT(vecTE,1),:,3);
        normals(:,:,2) = DVECT(EATT(vecTE,2),:,3);
        normals = mean(normals,3);
    else
        normals = DVECT(nonzeros(EATT(vecTE,:)),:,3);
    end
    
    R(end-len-(len2-1):end-len) = (4*pi).*dot(vuinf, normals,2);
    
else
    % some shit
end

end
