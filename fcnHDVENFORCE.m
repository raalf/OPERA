function matDVENFORCE = fcnHDVENFORCE(strATYPE, matUINF, matCONTROL, matDVECT, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, vecDVEAREA, matVLST, matDVE)
% Output is force per unit density

if strcmpi(strATYPE{2}, 'THIN') == 1
        
    tmp = -0.5.*[matCOEFF(:,4), matCOEFF(:,2), matCOEFF(:,5).*0];
    q_lm = fcnSTARGLOB(tmp, matROTANG);

    fpg = matCONTROL;
    q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    q_u = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    q_u = dot(matDVECT(:,:,3), q_u, 2).*matDVECT(:,:,3);
%     q_indw = q_lmu.*0;
    
    q_ind = q_lm + q_u + q_indw + matUINF;
%     q_ind = q_lm + q_u + matUINF;
    CPU = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    
    gca;
    hold on
    pt = matCENTER + 0.05.*matDVECT(:,:,3);
    quiver3(pt(:,1), pt(:,2), pt(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
    hold off
    
    q_ind = -q_lm + q_u + q_indw + matUINF;
%     q_ind = -q_lm + q_u + matUINF;
    CPL = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    
    gca;
    hold on
    pt = matCENTER - 0.05.*matDVECT(:,:,3);
    quiver3(pt(:,1), pt(:,2), pt(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'r')
    hold off
    
    matDVENFORCE = -(CPU - CPL).*0.5.*vecDVEAREA.*matDVECT(:,:,3);
    
    %     delta = 1e-2;
    %     % Above DVE
    %     fpg = matCONTROL + delta.*matDVECT(:,:,3);
    %     q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    %     q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    %     q_ind = q_inds + q_indw + matUINF;
    %     CPU = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    %     % Below DVE
    %     fpg = matCONTROL - delta.*matDVECT(:,:,3);
    %     q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
    %     q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    %     q_ind = q_inds + q_indw + matUINF;
    %     CPL = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    %
    %     matDVENFORCE = -0.5.*(CPU - CPL).*(sqrt(sum(matUINF.^2,2))).*vecDVEAREA.*matDVECT(:,:,3);
    
elseif strcmpi(strATYPE{2}, 'PANEL') == 1
    matDVENFORCE = nan(valNELE,1);
end

end