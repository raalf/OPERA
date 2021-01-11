function [matUINF, matVUINF, matUINF_KK] = fcnUINF(valUINF, valALPHA, valBETA, valNELE, valROTORS, matVLST, matKINCON_P, ...
    vecROTORRPM, vecDVEROTOR, matCENTER, matDVE, matROTORHUB, matROTORAXIS, vecKINCON_DVE, valROTSTART, vecROTRATE, valTIMESTEP, vecVEHORIG)

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
    
    orig_vec = vecVEHORIG - matROTORHUB(n,:);
    orig_vec = orig_vec/dcmXY2HUB;
    
    % timestep rotor in local XY hub plane
    tempROTORUINF_C = cross(repmat([0,0, -tmpROTORRADPS], length(tmpC(:,1)),1), tmpC);
    tempROTORUINF_V = cross(repmat([0,0, -tmpROTORRADPS], length(tmpV(:,1)),1), tmpV);
    tempROTORUINF_K = cross(repmat([0,0, -tmpROTORRADPS], length(tmpK(:,1)),1), tmpK);
    
    if valTIMESTEP > valROTSTART
        % This code is the right block:
        tempROTORUINF_C = tempROTORUINF_C + cross(repmat([-vecROTRATE(1), -vecROTRATE(2), 0], length(tmpC(:,1)),1), tmpC - orig_vec);
        tempROTORUINF_V = tempROTORUINF_V + cross(repmat([-vecROTRATE(1), -vecROTRATE(2), 0], length(tmpV(:,1)),1), tmpV - orig_vec);
        tempROTORUINF_K = tempROTORUINF_K + cross(repmat([-vecROTRATE(1), -vecROTRATE(2), 0], length(tmpK(:,1)),1), tmpK - orig_vec);
        
        %         % This is the wrong block:
        %         const_trans = cross([-vecROTRATE(1), -vecROTRATE(2), 0], -orig_vec);
        %         tempROTORUINF_C = tempROTORUINF_C + const_trans;
        %         tempROTORUINF_V = tempROTORUINF_V + const_trans;
        %         tempROTORUINF_K = tempROTORUINF_K + const_trans;
    end
    
    % transform rotor from xy plane to hub plane
    tempROTORUINF_C = tempROTORUINF_C*dcmXY2HUB;
    tempROTORUINF_V = tempROTORUINF_V*dcmXY2HUB;
    tempROTORUINF_K = tempROTORUINF_K*dcmXY2HUB;
    
    %     if valTIMESTEP > valROTSTART
    %         %         temp1(idxDVEROTOR,:) = temp1(idxDVEROTOR,:)*dcmXY2HUB;
    %         %         temp2(matDVE(idxDVEROTOR,:),:) = temp2(matDVE(idxDVEROTOR,:),:)*dcmXY2HUB;
    %         %         temp3(idx_k,:) = temp3(idx_k,:)*dcmXY2HUB;
    %         const_trans = const_trans*dcmXY2HUB;
    %     end
    
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
%     pts = [matCENTER; matVLST; matKINCON_P];
%     %     vels = [temp1; temp2; temp3];
%     vels = repmat(const_trans, size(pts,1), 1);
%     hold on
%     scatter3(vecVEHORIG(:,1), vecVEHORIG(:,2), vecVEHORIG(:,3), 200, 'sr')
%     quiver3(pts(:,1), pts(:,2), pts(:,3), vels(:,1), vels(:,2), vels(:,3))
%
%     hold off
%
%     %     hold on
%     %     quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), temp1(:,1), temp1(:,2), temp1(:,3), 10);
%     %     quiver3(matVLST(:,1), matVLST(:,2), matVLST(:,3), temp2(:,1), temp2(:,2), temp2(:,3), 10);
%     %     quiver3(matKINCON_P(:,1), matKINCON_P(:,2), matKINCON_P(:,3), temp3(:,1), temp3(:,2), temp3(:,3), 10);
%     %     hold off
% end

end

