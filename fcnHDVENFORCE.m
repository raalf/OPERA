function matDVENFORCE = fcnHDVENFORCE(strATYPE, matUINF, matCONTROL, matDVECT, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, vecDVEAREA, matVLST, matDVE)
% Output is force per unit density

if strcmpi(strATYPE{2}, 'THIN') == 1
    
    
    fpg = matCONTROL;
%     q_u = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
%     q_u = dot(matDVECT(:,:,3), q_u, 2).*matDVECT(:,:,3);
    q_u = fpg.*0;
    
    for i = 1:valNELE
        corners = matPLEX(:,:,2);
        points = fcnPOLYGRID(corners(:,1), corners(:,2), 10);
        
        len = size(points,1);
        vort = -0.5.*[matCOEFF(i,3).*points(:,1) + matCOEFF(i,4), matCOEFF(i,1).*points(:,2) + matCOEFF(i,2), points(:,2).*0];
        q_lm = fcnSTARGLOB(vort, repmat(matROTANG(i,:),len,1));
        CPU = 1 - (sqrt(sum((q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
        CPL = 1 - (sqrt(sum((-q_lm + matUINF(i,:) + q_u(i,:)).^2,2))).^2; % WON'T WORK FOR ROTORS
        matDVENFORCE(i,:) =  sum(-(CPU - CPL).*0.5.*(vecDVEAREA(i)./len).*matDVECT(i,:,3),1);
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
    matDVENFORCE = nan(valNELE,1);
end

end