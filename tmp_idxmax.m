nph_1 = 4;
nphs = arrayfun(@(c) length(c.ph), allCells);
idx4ph = find(nphs == nph_1);
idxs = zeros(360,length(idx4ph));

for cell_i = 1:length(idx4ph)
    osi = 1;
    for oi = 1:36;
        for si = 1:10;
            cell_idx = idx4ph(cell_i);            
            tc = squeeze( allCells(cell_idx).R(oi, si, :) );            
            idxs(osi, cell_i) = indmax_tmp(tc);
            osi = osi+1;
        end
    end
end

figure(342);
hist(idxs(:), nph_1*2-1);
nn = hist(idxs(:), nph_1);
nll = ones(1, nph_1)*nnz(~isnan(idxs))/nph_1;
p = histChiSqrTest(nn, nll);
title(sprintf('p = %.2f', p))