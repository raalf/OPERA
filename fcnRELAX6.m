function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
    fcnRELAX6(valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, matWE2GRID, matWEATT)

valTIMESTEP = valWNELE./(valWSIZE*2);
dvegrid = flipud(repmat([1 (1:valWSIZE)+valWSIZE], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');

tmp1 = fcnSDVEVEL(matWCENTER(dvegrid(:),:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 2e-1) + fcnSDVEVEL(matWCENTER(dvegrid(:),:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 2e-1);

tmp3 = zeros(size(matWELST,1),3);
tmp3(matWE2GRID(:),:) = tmp1;

hold on
quiver3(matWCENTER(dvegrid(:),1), matWCENTER(dvegrid(:),2), matWCENTER(dvegrid(:),3), tmp1(:,1), tmp1(:,2), tmp1(:,3));
hold off

% idx =  all(matWEATT(matWE2GRID,:),2);
% tmp3 = zeros(size(matWELST,1),3);
% 
% tmp3(matWE2GRID(idx),:) = (tmp1(matWEATT(matWE2GRID(idx),1),:) + tmp1(matWEATT(matWE2GRID(idx),2),:))./2;
% tmp3(matWE2GRID(~idx),:) = tmp1(nonzeros(matWEATT(matWE2GRID(~idx),:)),:);

tmp2 = zeros(size(matWVLST));

for j = 1:size(matWVGRID,1) - 1
    if j == 1
        tmp2(matWVGRID(j,:),:) = tmp3(matWE2GRID(j,:),:);
    else
        tmp2(matWVGRID(j,:),:) = (tmp3(matWE2GRID(j-1,:),:) + tmp3(matWE2GRID(j,:),:))./2;
    end
end


% Not moving some vertices
% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% % dont_move = unique([matWELST(vecWLE,:); matWELST(vecWTE,:); matWELST(oldest_edge,:)]);
% dont_move = unique([matWELST(vecWLE,:); matWELST(oldest_edge,:)]);

dont_move = unique([matWELST(vecWLE,:)]);

move = true(size(matWVLST,1),1);
move(dont_move) = false;
tmp2(dont_move,:) = tmp2(dont_move,:).*0;

%% Getting velocities at wake vertices
% tmp2 = zeros(size(matWVLST));
% tmp2(move,:) = fcnSDVEVEL(matWVLST(move,:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 5e-3) + fcnSDVEVEL(matWVLST(move,:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 5e-3);
% tmp2(matWELST(vecWSYM,:),2) = 0; % No y-component on symmetry line
%
% fore = [matWVGRID(1,:); matWVGRID(1:end-1,:)];
% aft = [matWVGRID(2:end,:); matWVGRID(end,:)];
%
% tmp2(matWVGRID,:) = (tmp2(fore,:) + 2.*tmp2(matWVGRID,:) + tmp2(aft,:))./4;

%% Moving
% Moving vertices
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
% matWCENTER
matWCENTER = (matWVLST(matWDVE(:,1),:) + matWVLST(matWDVE(:,2),:) + matWVLST(matWDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matWVLST(matWDVE(:,2),:) - matWVLST(matWDVE(:,3),:), matWVLST(matWDVE(:,1),:) - matWVLST(matWDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
[matWPLEX, matWDVECT, matWROTANG] = fcnTRITOLEX(P, DNORM, matWCENTER);

end


