function [matUINF, matVUINF, matUINF_KK] = fcnUINF(valUINF, valALPHA, valBETA, valNELE, valROTORS, matVLST, matKINCON_P, ...
    vecROTORRPM, vecDVEROTOR, matCENTER, matDVE, matROTORHUB, matROTORAXIS, vecKINCON_DVE, valROTSTART, vecROTRATE, valTIMESTEP)

% Translational velocities
matUINF = repmat(fcnUINFWING(valALPHA, valBETA), valNELE, 1).*valUINF;
matVUINF = repmat(fcnUINFWING(valALPHA, valBETA), size(matVLST,1), 1).*valUINF;
matUINF_KK = repmat(fcnUINFWING(valALPHA, valBETA), size(matKINCON_P,1), 1).*valUINF;

% Rotational velocities
for n = 1:valROTORS
    tmpROTORRADPS = vecROTORRPM(n).*2.*pi./60;
    
    idxDVEROTOR = vecDVEROTOR == n;
    
    tmpC = matCENTER(idxDVEROTOR,:);
    tmpV = matVLST(matDVE(idxDVEROTOR,:),:);
    idx_k = ismember(vecKINCON_DVE,find(idxDVEROTOR));
    tmpK = matKINCON_P(idx_k,:);
    
    tmpC = tmpC - matROTORHUB(n,:);
    tmpV = tmpV - matROTORHUB(n,:);
    tmpK = tmpK - matROTORHUB(n,:);
    
    dcmXY2HUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], matROTORAXIS(n,:))));
    
    % transform rotor from hub plane to xy plane
    tmpC = tmpC/dcmXY2HUB;
    tmpV = tmpV/dcmXY2HUB;
    tmpK = tmpK/dcmXY2HUB;
    
    if valTIMESTEP > valROTSTART
        tempROTORUINF_C = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), -tmpROTORRADPS], length(tmpC(:,1)),1), tmpC);
        tempROTORUINF_V = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), -tmpROTORRADPS], length(tmpV(:,1)),1), tmpV);
        tempROTORUINF_K = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), -tmpROTORRADPS], length(tmpK(:,1)),1), tmpK);
        
%         temp1(idxDVEROTOR,:) = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), vecROTRATE(3)], length(tmpC(:,1)),1), tmpC)*dcmXY2HUB;
%         temp2(matDVE(idxDVEROTOR,:),:) = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), vecROTRATE(3)], length(tmpV(:,1)),1), tmpV)*dcmXY2HUB;
%         temp3(idx_k,:) = cross(repmat([-vecROTRATE(1), -vecROTRATE(2), vecROTRATE(3)], length(tmpK(:,1)),1), tmpK)*dcmXY2HUB;
    else
        % timestep rotor in local XY hub plane
        tempROTORUINF_C = cross(repmat([0,0, -tmpROTORRADPS], length(tmpC(:,1)),1), tmpC);
        tempROTORUINF_V = cross(repmat([0,0, -tmpROTORRADPS], length(tmpV(:,1)),1), tmpV);
        tempROTORUINF_K = cross(repmat([0,0, -tmpROTORRADPS], length(tmpK(:,1)),1), tmpK);
    end
    
    % transform rotor from xy plane to hub plane
    tempROTORUINF_C = tempROTORUINF_C*dcmXY2HUB;
    tempROTORUINF_V = tempROTORUINF_V*dcmXY2HUB;
    tempROTORUINF_K = tempROTORUINF_K*dcmXY2HUB;
    
    % write rotated rotor to matVLST
    matUINF(idxDVEROTOR,:) = matUINF(idxDVEROTOR,:) + tempROTORUINF_C;
    matVUINF(matDVE(idxDVEROTOR,:),:) = matVUINF(matDVE(idxDVEROTOR,:),:) + tempROTORUINF_V;
    matUINF_KK(idx_k,:) = matUINF_KK(idx_k,:) + tempROTORUINF_K;
end

% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), matUINF(:,1), matUINF(:,2), matUINF(:,3));
% quiver3(matVLST(:,1), matVLST(:,2), matVLST(:,3), matVUINF(:,1), matVUINF(:,2), matVUINF(:,3));
% quiver3(matKINCON_P(:,1), matKINCON_P(:,2), matKINCON_P(:,3), matUINF_KK(:,1), matUINF_KK(:,2), matUINF_KK(:,3));
% hold off

% if valTIMESTEP > valROTSTART
%     hold on
%     quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), temp1(:,1), temp1(:,2), temp1(:,3));
%     quiver3(matVLST(:,1), matVLST(:,2), matVLST(:,3), temp2(:,1), temp2(:,2), temp2(:,3));
%     quiver3(matKINCON_P(:,1), matKINCON_P(:,2), matKINCON_P(:,3), temp3(:,1), temp3(:,2), temp3(:,3));
%     hold off
% end

end

