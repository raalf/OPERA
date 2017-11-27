function [infl_glob] = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER)

dvenum = reshape(dvenum, [], 1, 1); % Ensuring dvenum is a column vector

len = length(dvenum);

pa = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

p1 = permute(matPLEX(1,:,dvenum),[3 2 1]);
p2 = permute(matPLEX(2,:,dvenum),[3 2 1]);
p3 = permute(matPLEX(3,:,dvenum),[3 2 1]);

fn = cross(p3 - p1, p2 - p1);

N = zeros(len,3,3);
N(:,:,1) = cross(fn, p2-p1, 2);
N(:,:,2) = cross(fn, p3-p2, 2);
N(:,:,3) = cross(fn, p1-p3, 2);
N = N./sum(sqrt(N.^2),2);

%% (b1,b2,b3)
% b1 = (p3 - p1);
% b3 = cross(p3 - p1, p2 - p1);
% b2 = cross(b3, b1);
% b1 = b1./sqrt(sum(b1.^2,2));
% b2 = b2./sqrt(sum(b2.^2,2));
% b3 = b3./sqrt(sum(b3.^2,2));

b1 = matDVECT(dvenum,:,1);
b2 = matDVECT(dvenum,:,2);
b3 = matDVECT(dvenum,:,3);

%% Projection of P_A on plane
pb = pa - dot((pa - p1), b3, 2).*b3;

z = dot((pa - pb), b3, 2);

%% Influence from S1 S2 and S3
q1 = p1 - pb;
q2 = p2 - pb;
q3 = p3 - pb;

c3 = b3;

delta = 0.000001;
h = sqrt(z.^2 + delta.^2);

b(:,:,1) = b1;
b(:,:,2) = b2;
b(:,:,3) = b3;

% S1 (pb,p1,p2)
[infl_1] = fcnSNINF(q1, q2, c3, b, z, h);

% S2 (pb,p2,p3)
[infl_2] = fcnSNINF(q2, q3, c3, b, z, h);

% S3 (pb,p3,p1)
[infl_3] = fcnSNINF(q3, q1, c3, b, z, h);

%% Combining S1, S2, S3 to make S

comb = [dot(N(:,:,1), q1, 2) dot(N(:,:,2),q2,2) dot(N(:,:,3),q3,2)];
comb = comb./abs(comb);
comb = reshape(comb',1,3,[]);
comb(isnan(comb)) = 0;

infl_1 = infl_1.*repmat(permute(comb(:,1,:),[2 1 3]),3,6,1);
infl_2 = infl_2.*repmat(permute(comb(:,2,:),[2 1 3]),3,6,1);
infl_3 = infl_3.*repmat(permute(comb(:,3,:),[2 1 3]),3,6,1);

infl = infl_1 + infl_2 + infl_3;

dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);

infl_glob = fcnSTARGLOB(reshape(permute(infl,[2 3 1]),[],3,1), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));
infl_glob = reshape(infl_glob',3,6,[]);

end
