function [CL, CDi, e, CY, CLf, CLi, CYf, CYi] = fcnWFORCES(flagPRINT, valTIMESTEP, matF_FS, matF_IF, matF_ID, vecLIFT_DIR, vecSIDE_DIR, vecDRAG_DIR, valDENSITY, valAREA, valSPAN, vecTSITER, vecRSQUARED)

len = size(matF_FS, 1);

CLf = sum(dot(matF_FS(:,:,valTIMESTEP), repmat(vecLIFT_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);
CLi = sum(dot(matF_IF(:,:,valTIMESTEP), repmat(vecLIFT_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);

CL = CLf + CLi;

CYf = sum(dot(matF_FS(:,:,valTIMESTEP), repmat(vecSIDE_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);
CYi = sum(dot(matF_IF(:,:,valTIMESTEP), repmat(vecSIDE_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);

CY = CYf + CYi;

CDi = sum(dot(matF_ID(:,:,valTIMESTEP), repmat(vecDRAG_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);

e = (CL.^2)./(pi.*((valSPAN.^2)./valAREA).*CDi);

% special_drag = sum(dot(matF_IF(:,:,valTIMESTEP), repmat(vecDRAG_DIR, len, 1), 2))./(0.5.*valDENSITY.*valAREA);
% special_lift = special_drag.*tand(10);

if flagPRINT
    fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\t\tIterations = %d\t\tR^2(ie,kk) = (%g, %g)\n', valTIMESTEP, CL, CDi, e, vecTSITER(valTIMESTEP), vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2));
end
end