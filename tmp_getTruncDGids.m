
do_gids = getAllGids('o');
S_o = load('indivClustersPruned_GLFcuw8_grating_dOr');

do_nans = zeros(size(do_gids));
progressBar('init-', length(do_gids), 30);
for i = 1:length(do_gids)    
    progressBar(i);
    Gid = do_gids(i);
    
    sd = siteDataFor('Gid', Gid, 1);
    if length(sd.cellIds) < 2
        continue;
    end
    cellId = sd.cellIds(2);
    nm = getName('celldata', Gid, cellId);
    
    v = S_o.(nm);    
    do_nans(i) = any(isnan( v.OSP.R_full(:)) );    
end
progressBar('done');
Gids_o = do_gids(find(do_nans));

%%%

ds_gids = getAllGids('s');
S_s = load('indivClustersPruned_GLFcuw8_grating_dSf');

ds_nans = zeros(size(ds_gids));
progressBar('init-', length(ds_gids), 30);
for i = 1:length(ds_gids)    
    progressBar;
    Gid = ds_gids(i);
    sd = siteDataFor('Gid', Gid, 1);
    if length(sd.cellIds) < 2
        continue;
    end
    cellId = sd.cellIds(2);
    nm = getName('celldata', Gid, cellId);
    v = S_s.(nm);
    do_nans(i) = any(isnan( v.OSP.R_full(:)) );    
end
progressBar('done');
Gids_s = ds_gids(find(ds_nans));
