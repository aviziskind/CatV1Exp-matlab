

nph = arrayfun(@(s) size(s.R, 3), S.allCells);
nph_use = 60;
idx_use = find(nph == nph_use);

n = length(idx_use);

ind_pref = zeros(1, n);
pref_phase = zeros(1, n);

for i = 1:n
    R = S.allCells(idx_use(i)).R;
    
    R_os = mean(R, 3);
    [mxVal, i_max] = maxElement(R_os);    
    ptc = squeeze( R(i_max(1), i_max(2), :));
    
    [mx_val, ind_pref(i)] = max(ptc);
    if nnz(mx_val == ptc) > 1
        ind_pref(i) = nan;
    end
    pref_phase(i) = (ind_pref(i)-.5)*360/nph_use;
        
end

%%
nBin_ph = 27;
binE = linspace(0, 360, nBin_ph+1);
binC = binEdge2cent(binE);
binV = histcnt(pref_phase, binE);

binC_ph = binC;
frame_ms = 1000/120;
binC_t_ms = frame_ms * [1:nBin_ph] * nph_use/nBin_ph;

figure(1);
bar(binC_t_ms, binV, 1);
% xlim([0 binC_t_ms(end)+30]);
bar(binC_ph, binV, 1);
xlim([0 360]);
ylim([0 60])
set(gca, 'xtick', 0:45:360)
xlabel('Preferred Phase (at pref ori/spf)'); ylabel('Count')
title('Drifting Gratings (2 Hz)')

3;