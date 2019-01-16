function [matDVENFORCE, matDVEFDIST, matDVEFDIST_P] = fcnHDVENFORCE(strATYPE, matUINF, valDENSITY, matCONTROL, matDVECT, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, vecDVEAREA, matVLST, matDVE)
% Output is force per unit density

matDVEFDIST = nan;
matDVEFDIST_P = nan;

if strcmpi(strATYPE{2}, 'THIN') == 1
 
    fpg = matCONTROL;
    w_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    q_u = dot(matUINF + w_ind, matDVECT(:,:,3), 2).*matDVECT(:,:,3);

    for i = 1:valNELE
        corners = matPLEX(:,:,i);
        points = fcnPOLYGRID(corners(:,1), corners(:,2), 30);
        
        len = size(points,1);
        q_lm = -0.5.*[matCOEFF(i,3).*points(:,1) + matCOEFF(i,4), matCOEFF(i,1).*points(:,2) + matCOEFF(i,2), points(:,2).*0];
        q_lm = fcnSTARGLOB(q_lm, repmat(matROTANG(i,:),len,1));
        
        CPU = 1 - (sqrt(sum((q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
        CPL = 1 - (sqrt(sum((-q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
        
%         delta_p = valDENSITY.*sqrt(sum(matUINF(i,:).^2,2)).*sqrt(sum(q_lm.^2,2));
        delta_p = -(CPU - CPL).*0.5.*valDENSITY.*sqrt(sum(matUINF(i,:).^2,2)).^2;
%         delta_f = (delta_p.*(vecDVEAREA(i)./len).*[0 0 1]);
%         matDVEFDIST(:,:,i) = fcnSTARGLOB(delta_f, repmat(matROTANG(i,:),len,1));
%         matDVEFDIST_P(:,:,i) = fcnSTARGLOB([points points(:,1).*0], repmat(matROTANG(i,:),len,1)) + matCENTER(i,:);
        matDVENFORCE(i,:) = sum(delta_p.*(vecDVEAREA(i)./len),1).*matDVECT(i,:,3);
        
%         CPU = 1 - (sqrt(sum((q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
%         CPL = 1 - (sqrt(sum((-q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
%         matDVENFORCE(i,:) =  sum( -(CPU - CPL).*0.5.*valDENSITY.*(vecDVEAREA(i)./len).*sqrt(sum(matUINF(i,:).^2,2)).^2.*matDVECT(i,:,3) ,1);
    end

    %     tmp = -0.5.*[matCOEFF(:,4), matCOEFF(:,2), matCOEFF(:,5).*0];
    %     q_lm = fcnSTARGLOB(tmp, matROTANG);
    %
    %     fpg = matCONTROL;
    %
    %     q_u = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    %     q_u = dot(matDVECT(:,:,3), q_u, 2).*matDVECT(:,:,3);
    %     q_indw = q_u.*0;
    %
    %     q_ind = q_lm + q_u + q_indw + matUINF;
    % %     q_ind = q_lm + q_u + matUINF;
    %     CPU = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    % %     gca;
    % %     hold on
    % %     pt = matCENTER + 0.05.*matDVECT(:,:,3);
    % %     quiver3(pt(:,1), pt(:,2), pt(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
    % %     hold off
    %
    %     q_ind = -q_lm + q_u + q_indw + matUINF;
    % %     q_ind = -q_lm + q_u + matUINF;
    %     CPL = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    % %     gca;
    % %     hold on
    % %     pt = matCENTER - 0.05.*matDVECT(:,:,3);
    % %     quiver3(pt(:,1), pt(:,2), pt(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'r')
    % %     hold off
    %
    %     matDVENFORCE = -(CPU - CPL).*0.5.*vecDVEAREA.*matDVECT(:,:,3);
    
    %         delta = 1e-2;
    %         % Above DVE
    %         fpg = matCONTROL + delta.*matDVECT(:,:,3);
    %         q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    %         q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    %         q_ind = q_inds + q_indw + matUINF;
    %         CPU = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    %         % Below DVE
    %         fpg = matCONTROL - delta.*matDVECT(:,:,3);
    %         q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    %         q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    %         q_ind = q_inds + q_indw + matUINF;
    %         CPL = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    %         matDVENFORCE = -0.5.*(CPU - CPL).*(sqrt(sum(matUINF.^2,2))).*vecDVEAREA.*matDVECT(:,:,3);
    
elseif strcmpi(strATYPE{2}, 'PANEL') == 1
    fpg = matCONTROL;
    w_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
    q_u = dot(matUINF + w_ind, matDVECT(:,:,3), 2).*matDVECT(:,:,3);
    
    for i = 1:valNELE
        corners = matPLEX(:,:,2);
        points = fcnPOLYGRID(corners(:,1), corners(:,2), 10);
        
        len = size(points,1);
        vort = -0.5.*[matCOEFF(i,3).*points(:,1) + matCOEFF(i,4), matCOEFF(i,1).*points(:,2) + matCOEFF(i,2), points(:,2).*0];
        q_lm = fcnSTARGLOB(vort, repmat(matROTANG(i,:),len,1));
        CP = 1 - (sqrt(sum((q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
%         CPL = 1 - (sqrt(sum((-q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
        matDVENFORCE(i,:) =  sum( -(CP).*0.5.*valDENSITY.*(vecDVEAREA(i)./len).*sqrt(sum(matUINF(i,:).^2,2)).*matDVECT(i,:,3) ,1);
    end

%             delta = 1e-10;
%             % Above DVE
%             fpg = matCONTROL + delta.*matDVECT(:,:,3);
%             q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
%             q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
%             q_ind = q_inds + q_indw + matUINF;
%             CP = 1 - (sqrt(sum(q_ind.^2,2)).^2);
%     
%     
%             matDVENFORCE = -0.5.*(CP).*(sqrt(sum(matUINF.^2,2))).*vecDVEAREA.*matDVECT(:,:,3);

end

end