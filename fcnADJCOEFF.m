function [matCOEFF] = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE)

enum = permute(matEIDX, [2 3 1]);
tmpR = [vecVMU(permute(matDVE, [2 3 1])); vecEMU(enum)];

vloc = permute(reshape(matVLST(permute(matDVE, [2 3 1]),:)', 3, 3, []), [2 1 3]);

v1 = permute(reshape(matVLST(reshape(matELST(enum, 1), 3, 1, []),:)', 3, 3, []), [2 1 3]);
v2 = permute(reshape(matVLST(reshape(matELST(enum, 2), 3, 1, []),:)', 3, 3, []), [2 1 3]);
eloc = (v1 + v2)./2;

pts = [vloc; eloc];

% Similar to fcnDCIRC2
dvenum = repmat(reshape(1:valNELE, 1, 1, []), 6, 1, 1);
dvenum = reshape(permute(dvenum, [2 1 3]), size(dvenum,2), [])';
pts = reshape(permute(pts, [2 1 3]), size(pts,2), [])';

% circ = fcnDCIRC2(pts, dvenum, valNELE, matROTANG, matCENTER);
% tmpR = reshape(permute(tmpR, [2 1 3]), 1, [])';
% matCOEFF = fcnSOLVED(circ, tmpR, valNELE);

pts = fcnGLOBSTAR(pts - matCENTER(dvenum,:), matROTANG(dvenum,:));
pts = permute(reshape(pts', 3, 6, []), [2 1 3]);

a2 = pts(:,2,:);
a1 = 0.5.*(a2.^2);
b2 = pts(:,1,:);
b1 = 0.5.*(b2.^2);
c2 = pts(:,1,:).*pts(:,2,:);
c3 = ones(6,1,valNELE);
gamma = [a1, a2, b1, b2, c2, c3];

for i = 1:valNELE
    matCOEFF(i,:) = [gamma(:,:,i)\tmpR(:,:,i)]';
end

end