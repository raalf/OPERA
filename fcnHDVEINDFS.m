function [infl_glob] = fcnHDVEINDFS(dvenum_all, fpg_all, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matCENTER)

flagPLOT = 0;

dvenum_all = reshape(dvenum_all, [], 1, 1); % Ensuring dvenum is a column vector

chunk_sz = 4e6;
num_pts = length(dvenum_all);

infl_glob = nan(3,6, num_pts);

for i = 1:chunk_sz:num_pts
    
    if num_pts <= chunk_sz
        idx_chunk = [1:num_pts];
    elseif i + chunk_sz - 1 <= num_pts
        idx_chunk = [i:i + chunk_sz - 1];
    else
        idx_chunk = [i:num_pts];
    end
    
    fpg = fpg_all(idx_chunk,:);
    dvenum = dvenum_all(idx_chunk,:);
    len = length(dvenum);

    pa = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

    p1 = permute(matPLEX(1,:,dvenum),[3 2 1]);
    p2 = permute(matPLEX(2,:,dvenum),[3 2 1]);
    p3 = permute(matPLEX(3,:,dvenum),[3 2 1]);

    %% (b1,b2,b3)
    b = [];
    b(:,:,3) = repmat([0 0 1], length(dvenum),1);
    b(pa(:,3) < 0,:,3) = b(pa(:,3) < 0,:,3).*-1;
    b(:,:,1) = repmat([1 0 0], length(dvenum),1); 
    b(:,:,2) = cross(b(:,:,3), b(:,:,1), 2);

    %% Projection of P_A on plane
    pb = pa - dot((pa - p1), b(:,:,3), 2).*b(:,:,3);
    z = dot((pa - pb), b(:,:,3), 2);
    
    %% Transform into Local-Local
    loc_rot = [-atan2(b(:,2,3), b(:,3,3)), asin(b(:,1,3)), atan2(b(:,2,1), b(:,1,1))];

%     ROLL = -atan2(b(:,2,3), b(:,3,3));    
%     PITCH = asin(b(:,1,3));
%     YAW = atan2(b(:,2,1), b(:,1,1));
    
    p1 = fcnGLOBSTAR(p1 - pb, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
    p2 = fcnGLOBSTAR(p2 - pb, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
    p3 = fcnGLOBSTAR(p3 - pb, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
        
    %% Edge normals
    N = nan(len,3,3);
    N(:,:,1) = cross(b(:,:,3), p2-p1, 2);
    N(:,:,2) = cross(b(:,:,3), p3-p2, 2);
    N(:,:,3) = cross(b(:,:,3), p1-p3, 2);
    N = N./sum(sqrt(N.^2),2);
    
    %% Influence from S1 S2 and S3    
    q = nan(len,3,3);
    q(:,:,1) = p1;
    q(:,:,2) = p2;
    q(:,:,3) = p3;

    c3 = b(:,:,3);

    delta = 0.000001;
    h = sqrt(z.^2 + delta.^2);
           
    ax2 = [];
    if flagPLOT == 1 && len == 1; temp_plt_s; end

    dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
    
    infl = nan(3,6,len,size(q,3));
    for ii = 1:size(q,3)
        n = mod(ii,3) + 1;
        infl(:,:,:,ii) = fcnSNINF(q(:,:,ii), q(:,:,n), c3, z, h, N(:,:,ii), flagPLOT, ii, pb, ax2, loc_rot);
        infl(:,:,:,ii) = reshape(fcnSTARGLOB(reshape(permute(infl(:,:,:,ii),[2 3 1]),[],3,1), loc_rot(dvenum,1), loc_rot(dvenum,2), loc_rot(dvenum,3))',3,6,[]);
    end

    % Combining S1, S2, S3 to make S
    infl_tot = sum(infl,4);

    infl_tot = fcnSTARGLOB(reshape(permute(infl_tot,[2 3 1]),[],3,1), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));
    infl_tot(isnan(infl_tot)) = 0;
        
    infl_glob(:,:,idx_chunk) = reshape(infl_tot',3,6,[]);
end

end
