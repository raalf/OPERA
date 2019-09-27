clc
clear


load('Alpha 15 Results/fixed_coarse.mat')

matDGAMMADT = zeros(1,size(matLIFTFREE,2));
for i = 2:size(matLIFTFREE,1)-1
    matDGAMMADT(i,:) = (matINTCIRC(i+1,:) - matINTCIRC(i,:))./valDELTIME;
end
matDGAMMADT(i+1,:) = matDGAMMADT(i,:);

tmpLIFTFREE = matLIFTFREE + matDGAMMADT;
tmpDVELIFT = tmpLIFTFREE + matLIFTIND;
tmpDVEDRAG = matDRAGIND;
tmpDVESIDE = matSIDEFREE + matSIDEIND;
for i = 1:size(matLIFTFREE,1)
    tmpDVETHRUST(i,:) = dot(matDVELIFT_DIR(:,:,i).*tmpDVELIFT(i,:)', repmat([0 0 1], length(tmpDVELIFT(i,:)), 1), 2) ...
        + dot(matDVEDRAG_DIR(:,:,i).*tmpDVEDRAG(i,:)', repmat([0 0 1], length(tmpDVEDRAG(i,:)), 1), 2) ...
        + dot(matDVESIDE_DIR(:,:,i).*tmpDVESIDE(i,:)', repmat([0 0 1], length(tmpDVESIDE(i,:)), 1), 2);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    
end

save('coarse_fixed_2.mat');

figure(20);
plot(CT, '-k')
hold on
plot(CT_U, '--b');
hold off
grid minor
box on
axis tight