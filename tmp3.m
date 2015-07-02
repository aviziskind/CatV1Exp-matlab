figure(114);
% Gids = [2344        2346        2352        2354];
Gids = [1017, 1166]
m = length(Gids); n = 6;
for gi = 1:m
    [cellIds, nSpikes] = dbLookupNumSpikes(Gids(gi));
    nCells = length(cellIds);    
    for ci = 1:nCells
        nm = getName('celldata', Gids(gi), cellIds(ci));
        vr = eval(nm);
        subplot(m,n, (gi-1)*n+ci);
        OSP = vr.OSP;
        imagesc(mean(OSP.R, 3));
        set(gca, 'xtick', [], 'ytick', [])
        if ci == 1
            ylabel({sprintf('Gid = %d', Gids(gi)), sprintf('nspk = %d', sum(nSpikes))  });
        end
        title(sprintf('n = %d', nSpikes(ci)));
        rep_p = OSP.stats.rep_ori_sp_avPhase_pval;
        xlabel(sprintf('rep = %.2f', -log10(rep_p)));
%         3; vr.OSP.R
    end
end

        