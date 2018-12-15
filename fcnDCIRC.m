function [circ] = fcnDCIRC(pts, dvenum, valNELE, matROTANG, matCENTER)

nedg = size(pts,1);

pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));
pts_two = fcnGLOBSTAR(pts(:,:,2) - matCENTER(dvenum(:,2),:), matROTANG(dvenum(:,2),:));

a2 = [pts_one(:,2) pts_two(:,2)];
a1 = 0.5.*(a2.^2);

b2 = [pts_one(:,1) pts_two(:,1)];
b1 = 0.5.*(b2.^2);

c3 = ones(nedg,2);

gamma1 = [a1(:,1), a2(:,1), b1(:,1), b2(:,1), c3(:,1)];
gamma2 = [a1(:,2), a2(:,2), b1(:,2), b2(:,2), c3(:,2)].*-1;

if any(dvenum(:,2) - dvenum(:,1))
    circ = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum(:,1), dvenum(:,2), gamma1, gamma2);
elseif ~any(dvenum(:,2) - dvenum(:,1))   
    circ = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum(:,1), [], gamma1 + gamma2, []);  
else
    disp('Issue in fcnDCIRC');
end

end

