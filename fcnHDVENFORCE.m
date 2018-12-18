function matDVENFORCE = fcnHDVENFORCE(strATYPE, matUINF, matCONTROL, matDVECT, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, vecDVEAREA)
% Output is force per unit density

delta = 1e-8;

if strcmpi(strATYPE{2}, 'THIN') == 1

    % Above DVE
    fpg = matCONTROL + delta.*matDVECT(:,:,3);
    q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
        q_indw = q_inds.*0;
%     q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    q_ind = q_inds + q_indw + matUINF;
    CPU = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    
    % Below DVE
    fpg = matCONTROL - delta.*matDVECT(:,:,3);
    q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
%     q_indw = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
    q_ind = q_inds + q_indw + matUINF;
    CPL = 1 - (sqrt(sum(q_ind.^2,2)).^2);
    matDVENFORCE = -0.5.*(CPU - CPL).*(sqrt(sum(matUINF.^2,2))).*vecDVEAREA.*matDVECT(:,:,3);

elseif strcmpi(strATYPE{2}, 'PANEL') == 1
    matDVENFORCE = nan(valNELE,1);
end

end