function [CL, CDi, e, CY, CLf, CLi, CYf, CYi] = fcnWFORCES(matF_FS, matF_IF, matF_ID, matLIFT_DIR, matSIDE_DIR, matDRAG_DIR, valDENSITY, valAREA, valSPAN)

CLf = sum(dot(matF_FS, matLIFT_DIR, 2))./(0.5.*valDENSITY.*valAREA);
CLi = sum(dot(matF_IF, matLIFT_DIR, 2))./(0.5.*valDENSITY.*valAREA);

CL = CLf + CLi;

CYf = sum(dot(matF_FS, matSIDE_DIR, 2))./(0.5.*valDENSITY.*valAREA);
CYi = sum(dot(matF_IF, matSIDE_DIR, 2))./(0.5.*valDENSITY.*valAREA);

CY = CYf + CYi;

CDi = sum(dot(matF_ID, matDRAG_DIR, 2))./(0.5.*valDENSITY.*valAREA);

e = (CL.^2)./(pi.*((valSPAN.^2)./valAREA).*CDi);

end