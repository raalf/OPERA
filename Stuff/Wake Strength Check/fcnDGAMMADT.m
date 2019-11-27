function [CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR)

% matDGAMMADT = zeros(size(matLIFTFREE));
% len = size(matLIFTFREE,1);
% for i = 1:skip:len
%     if i <= len-(2*skip)
%         matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
%     elseif i > len-(2*skip) && i <= len - skip
%         matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
%     else
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
%     end
% end
% matDGAMMADT(1,:) = matDGAMMADT(1,:).*0;

if size(matINTCIRC,1) > 10
for i = 1:size(matINTCIRC,2)
    
    matINTCIRC(:,i) = smooth(matINTCIRC(:,i),'rloess');
    
end
end

matDGAMMADT = nan(size(matLIFTFREE));
len = size(matLIFTFREE,1);
for i = 1:skip:len
    if i <= len-(2*skip)
        matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
    elseif i > len-(2*skip) && i <= len - skip
        matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
    else
        matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
    end
end
matDGAMMADT(1,:) = zeros(size(matDGAMMADT(1,:)));

% matDGAMMADT = zeros(size(matLIFTFREE));
% len = size(matLIFTFREE,1);
% for i = 1:skip:len
%     if i <= skip
%         matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
%     elseif i > skip && i <= len-(2*skip)
%         matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 6.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:) - 2.*matINTCIRC(i-skip,:))./(6.*valDELTIME.*skip);
%     elseif i > len-(2*skip) && i <= len - skip
%         matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
%     else
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
%     end
% end

% X = [1:size(matINTCIRC, 1)]'.*valDELTIME;
% Y = matINTCIRC(:,8);

% matDGAMMADT = nan(size(matLIFTFREE));
% for i = 1:size(matDGAMMADT, 2)
%     X = [1:size(matDGAMMADT, 1)]'.*valDELTIME;
%     [xData, yData] = prepareCurveData( X, matINTCIRC(:,i) );
%     
%     % Set up fittype and options.
%     ft = fittype( 'smoothingspline' );
%     opts = fitoptions( 'Method', 'SmoothingSpline' );
%     opts.SmoothingParam = 0.9999;
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     
%     matDGAMMADT(:,i) = differentiate(fitresult, X);
%     
% end
% matDGAMMADT(1,:) = zeros(size(matDGAMMADT(1,:)));

tmpLIFTFREE = matLIFTFREE + matDGAMMADT;
tmpDVELIFT = tmpLIFTFREE + matLIFTIND;
tmpDVEDRAG = matDRAGIND;
tmpDVESIDE = matSIDEFREE + matSIDEIND;

for i = 1:skip:size(matLIFTFREE,1)
    tmpDVETHRUST(i,:) = dot(matDVELIFT_DIR(:,:,i).*tmpDVELIFT(i,:)', repmat([0 0 1], length(tmpDVELIFT(i,:)), 1), 2) ...
        + dot(matDVEDRAG_DIR(:,:,i).*tmpDVEDRAG(i,:)', repmat([0 0 1], length(tmpDVEDRAG(i,:)), 1), 2) ...
        + dot(matDVESIDE_DIR(:,:,i).*tmpDVESIDE(i,:)', repmat([0 0 1], length(tmpDVESIDE(i,:)), 1), 2);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    
end

CL_U = sum(tmpLIFTFREE,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2));
CL_U = CL_U(~isnan(CL_U));

end



% len = size(matLIFTFREE,1);
% lambda = 0.5;
% matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
% for i = 3:len
%     if i <= len-2
%         matDGAMMADT(i,:) = ((-matINTCIRC(i+2,:) + 8.*matINTCIRC(i+1,:) - 8.*matINTCIRC(i-1,:) + matINTCIRC(i-2,:))./(12.*valDELTIME));
%     elseif i == len - 1
%         matDGAMMADT(i,:) = ((matINTCIRC(i+1,:) - matINTCIRC(i-1,:))./(2.*valDELTIME));
%     elseif i == len
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME);
%     end
% end

% len = size(matLIFTFREE,1);
% lambda = 0.5;
% matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
% for i = 2:len
%     matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME).*lambda + (1 - lambda).*matDGAMMADT(i-1,:);
% end