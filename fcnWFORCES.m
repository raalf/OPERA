function [CL, CDi, CY, e] = fcnWFORCES(valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, valDENSITY, valAREA, valSPAN, vecTSITER, vecRSQUARED)

CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*valAREA);
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*valAREA);
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*valAREA);
% e = (CL.^2)./(pi.*((valSPAN.^2)./valAREA).*CDi);
q_inf = 0.5.*valDENSITY;
e = ((nansum(vecDVELIFT)./q_inf).^2)./(pi.*(valSPAN.^2).*(nansum(vecDVEDRAG)./q_inf));

fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\t\tIterations = %d\t\tR^2(ie,kk) = (%g, %g)\n', valTIMESTEP, CL, CDi, e, vecTSITER(valTIMESTEP), vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2));

end