function [vort_edge] = fcnDVORT2(pts, dvenum, valNELE, matCENTER, matROTANG, direction)

nedg = length(dvenum);
zer = zeros(nedg,1);

pts = fcnGLOBSTAR(pts - matCENTER(dvenum,:), matROTANG(dvenum,:));

if strcmpi(direction,'B')
    dgamma1 = [zer, zer, pts(:,1), ones(nedg,1), zer];
elseif strcmpi(direction,'A')
    dgamma1 = [pts(:,2), ones(nedg,1), zer, zer, zer];
else
    disp('Issue in fcnDVORT1')
end

vort_edge = real(fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []));

end