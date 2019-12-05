function [infl_loc] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBI, ztol)

chunk_sz = 1e7;
num_pts = length(dvenum);

GPU = false;

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
        if GPU == false
            infl_loc(:,:,idx_chunk) = fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), fpg(idx_chunk,:), matPLEX, matROTANG, matCONTROL, vecBI(idx_chunk), ztol, GPU);
        else
            infl_loc(:,:,idx_chunk) = gather(fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), gpuArray(fpg(idx_chunk,:)), gpuArray(matPLEX), gpuArray(matROTANG), gpuArray(matCONTROL), vecBI(idx_chunk), ztol, GPU)); 
        end
    else
        if GPU == false
            infl_loc(:,:,idx_chunk) = fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), fpg(idx_chunk,:), matPLEX, matROTANG, matCONTROL, [], ztol, GPU);
        else
            infl_loc(:,:,idx_chunk) = gather(fcnHDVEIND_DB(dvenum(idx_chunk), dvetype(idx_chunk), gpuArray(fpg(idx_chunk,:)), gpuArray(matPLEX), gpuArray(matROTANG), gpuArray(matCONTROL), [], ztol, GPU)); 
        end
    end
end

end