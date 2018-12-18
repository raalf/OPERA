function [vort_edge] = fcnDVORT1(dvenum, valNELE, direction)

nedg = length(dvenum);
zer = zeros(nedg,1);

if strcmpi(direction,'B')
    dgamma1 = [zer, zer, ones(nedg,1), zer, zer];
    dgamma2 = [zer, zer zer, ones(nedg,1), zer];
elseif strcmpi(direction,'A')
    dgamma1 = [ones(nedg,1), zer, zer, zer, zer];
    dgamma2 = [zer, ones(nedg,1), zer, zer, zer];
elseif strcmpi(direction,'C')
    dgamma1 = [zer, zer, zer, zer, ones(nedg,1)];
    dgamma2 = [zer, zer, zer, zer, ones(nedg,1)];
else
    disp('Issue in fcnDVORT1')
end

vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma2, []);
vort_edge = real([vort_edge1; vort_edge2]);

end