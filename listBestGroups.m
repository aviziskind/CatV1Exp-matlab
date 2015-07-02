load flashedGratingCells.mat
allGids = [allCells.Gid];
[uGids, gidList] = uniqueList(allGids);
nCellsInGrp = cellfun(@length, gidList);
nSpks = [allCells.nspikes];
cc_ps = arrayfun(@(s) s.stats.allWindowStats.cc_p, allCells);
grp_ccs = cellfun(@(idx) cc_ps(idx), gidList, 'un', 0);

med_ccs = cellfun(@median, grp_ccs);
idx = ord([med_ccs .* nCellsInGrp], 'descend');
n = 60;
fprintf('Best %d groups:\n', n)
N = length(idx)
gid_idx = [1:n/2, N:-1:N-n/2];
for i = gid_idx
    j = idx(i);
    
    nCells = length(gidList{j});
    fmt_d = repmat('%d (%.2f),', 1, nCells);
    nCellSpikes = nSpks(gidList{j});
    cell_ccps = cc_ps(gidList{j});
    
    fprintf(['(%d) Group %d : %d cells [' fmt_d '] \n '], ...
        i, uGids(j), nCells, [nCellSpikes(:)'; cell_ccps(:)']);
    
end

bestGids = uGids(idx(gid_idx))'
clear allCells stats;