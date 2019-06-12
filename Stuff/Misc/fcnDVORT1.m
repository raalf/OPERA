function [vort_edge] = fcnDVORT1(dvenum, valNELE, direction)

nedg = length(dvenum);
zer = zeros(nedg,1);

if strcmpi(direction,'B')
    dgamma1 = [zer, zer, ones(nedg,1), zer, zer, zer];
    dgamma2 = [zer, zer, zer, ones(nedg,1), zer, zer];
elseif strcmpi(direction,'A')
    dgamma1 = [ones(nedg,1), zer, zer, zer, zer, zer];
    dgamma2 = [zer, ones(nedg,1), zer, zer, zer, zer];
elseif strcmpi(direction,'C')
    dgamma1 = [zer, zer, zer, zer, zer, ones(nedg,1)];
    dgamma2 = [zer, zer, zer, zer, zer, ones(nedg,1)];
else
    disp('Issue in fcnDVORT1')
end

dgamma3 = [zer, zer, zer, zer, ones(nedg,1), zer];

vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum, [], dgamma1, []);
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum, [], dgamma2, []);
vort_edge3 = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum, [], dgamma3, []);
vort_edge = real([vort_edge1; vort_edge2; vort_edge3]);

end