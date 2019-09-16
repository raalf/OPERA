function [infl_loc] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBI, ztol)

chunk_sz = 4e5;
num_pts = length(dvenum);


infl_loc = zeros(3,6,size(dvenum,1));
for i = 1:chunk_sz:num_pts
    idx_chunk = false(size(dvenum));
    if num_pts <= chunk_sz
        idx_chunk = true(size(dvenum));
    elseif i + chunk_sz - 1 <= num_pts
        idx_chunk(i:i + chunk_sz - 1) = true;
    else
        idx_chunk(i:num_pts) = true;
    end    
    
    if ~isempty(vecBI)
        infl_loc(:,:,idx_chunk) = fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), fpg(idx_chunk,:), matPLEX, matROTANG, matCONTROL, vecBI(idx_chunk), ztol);
    else
        infl_loc(:,:,idx_chunk) = fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), fpg(idx_chunk,:), matPLEX, matROTANG, matCONTROL, [], ztol);   
    end
end

end