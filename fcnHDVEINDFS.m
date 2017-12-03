function [infl_glob] = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER)

flagPLOT = 0;

dvenum = reshape(dvenum, [], 1, 1); % Ensuring dvenum is a column vector

len = length(dvenum);

pa = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

p1 = permute(matPLEX(1,:,dvenum),[3 2 1]);
p2 = permute(matPLEX(2,:,dvenum),[3 2 1]);
p3 = permute(matPLEX(3,:,dvenum),[3 2 1]);

%% (b1,b2,b3)
b1 = matDVECT(dvenum,:,1);
b2 = matDVECT(dvenum,:,2);
b3 = matDVECT(dvenum,:,3);

%% Edge normals
N = zeros(len,3,3);
N(:,:,1) = cross(b3, p2-p1, 2);
N(:,:,2) = cross(b3, p3-p2, 2);
N(:,:,3) = cross(b3, p1-p3, 2);
N = N./sum(sqrt(N.^2),2);

%% Projection of P_A on plane
pb = pa - dot((pa - p1), b3, 2).*b3;

z = dot((pa - pb), b3, 2);

%% Influence from S1 S2 and S3
q(:,:,1) = p1 - pb;
q(:,:,2) = p2 - pb;
q(:,:,3) = p3 - pb;

c3 = b3;

delta = 0.2;
h = sqrt(z.^2 + delta.^2);

b(:,:,1) = b1;
b(:,:,2) = b2;
b(:,:,3) = b3;

% Whether or not to add or subtract the triangle
comb = [dot(N(:,:,1),q(:,:,1),2), dot(N(:,:,2),q(:,:,2),2), dot(N(:,:,3),q(:,:,3),2)];
comb = comb./abs(comb);
comb = reshape(comb',1,3,[]);
comb(isnan(comb)) = 0;

if flagPLOT == 1 && len == 1; temp_plt_s; end

infl = nan(3,6,len,size(q,3));
for i = 1:size(q,3)
    m = i;
    if i == size(q,3)
        n = 1;
    else
        n = i+1;
    end
    
    infl(:,:,:,i) = fcnSNINF(q(:,:,m), q(:,:,n), c3, b, z, h).*repmat(permute(comb(:,i,:),[2 1 3]),3,6,1);
    
end

% Combining S1, S2, S3 to make S

infl_tot = sum(infl,4);

dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);

infl_glob = fcnSTARGLOB(reshape(permute(infl_tot,[2 3 1]),[],3,1), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));
infl_glob = reshape(infl_glob',3,6,[]);

infl_glob(isnan(infl_glob)) = 0;

end
