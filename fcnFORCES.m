function [CL, CDi, e, CY, CLf, CLi, CYf, CYi, CT, CMx, CMy, CMz] = fcnFORCES(flgPRINT, valTIMESTEP, valWINGS, valROTORS, ...
    matF_FS, matF_IF, matF_ID, matF_NC, matLIFT_DIR, matSIDE_DIR, matDRAG_DIR, valDENSITY, valAREA, valSPAN, valUINF, ...
    vecROTORRPM, vecROTORDIAM, matROTORAXIS, matROTORHUB, matCENTER, vecTSITER, vecRSQUARED, vecDVEROTOR, vecDVEWING, matROTORTRANS)

if isempty(matF_NC)
    matF_NC = zeros(size(matF_FS));
end

out_str = sprintf('Timestep: %d\t\t', valTIMESTEP);

if valWINGS > 0
    idxwing = vecDVEWING > 0;
    [CL, CDi, e, CY, CLf, CLi, CYf, CYi] = fcnWFORCES(matF_FS(idxwing,:,valTIMESTEP), matF_IF(idxwing,:,valTIMESTEP), ...
        matF_ID(idxwing,:,valTIMESTEP), matF_NC(idxwing,:,valTIMESTEP), matLIFT_DIR(idxwing,:), matSIDE_DIR(idxwing,:), matDRAG_DIR(idxwing,:), valDENSITY, valAREA, valSPAN, valUINF);
    
    out_str = [out_str, sprintf('CL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\t\t', CL, CDi, e)];
    
else
    CL = nan; CDi = nan; e = nan; CY = nan; CLf = nan; CLi = nan; CYf = nan; CYi = nan;
end


if valROTORS == 0
    CT = nan;
    CMx = nan;
    CMy = nan;
    CMz = nan;
else
    for n = 1:valROTORS
        idxrotor = vecDVEROTOR == n;
        [CT(1,n), ~, CMx(1,n), CMy(1,n), CMz(1,n)] = fcnRFORCES(matF_FS(idxrotor,:,valTIMESTEP), matF_IF(idxrotor,:,valTIMESTEP), matF_ID(idxrotor,:,valTIMESTEP), ...
            matF_NC(idxrotor,:,valTIMESTEP), matROTORAXIS(n,:), matROTORHUB(n,:), matCENTER(idxrotor,:), valDENSITY, vecROTORDIAM(n), vecROTORRPM(n), matROTORTRANS(n,:,:));
        
        out_str = [out_str, sprintf('CT_%d = %0.5f\t\tCMx_%d = %0.5f\t\tCMy_%d = %0.5f\t\tCMz_%d = %0.5f\t\t', n, CT(1,n), n, CMx(1,n), n, CMy(1,n), n, CMz(1,n))];
    end
end

out_str = [out_str, sprintf('Iterations = (%d, %d)\t\tR^2 = (%g, %g)\n', vecTSITER(valTIMESTEP,1), vecTSITER(valTIMESTEP,2), vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2))];

if flgPRINT
    disp(out_str);
end

end