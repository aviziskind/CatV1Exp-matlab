%%
Gids = ratDatafilesAvailable;
for i = 60:length(Gids)
    S = load(getFileName('properties', Gids(i)));
    nspk_prop = S.nSpikes;
    S_sort = load(getFileName('clusters', Gids(i)));
    nspk_sort = length(S_sort.clusterIds);
    if nspk_prop ~= nspk_sort
        fprintf('(%d) Gid = %d. %d vs %d\n', i, Gids(i), nspk_prop, nspk_sort);
    end
    
end