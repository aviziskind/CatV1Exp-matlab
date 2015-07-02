progressBar('init-', 799)
allflds = cell(1,length(gids));
wasPresent = false(1, length(gids));
raw_fld = 'channel_all_unnorm_raw';
for i = 1:length(gids)
    gid = gids(i);
    fn = getFileName('kNN', gid, 0);
    s = who('-file', fn);
    has_raw_fld = any(strcmp(s, raw_fld));
    wasPresent(i) = has_raw_fld;
    if has_raw_fld
        S = load(fn);
        assert (isfield(S, raw_fld))
        S = rmfield(S, raw_fld);
        save(fn, '-struct', 'S');
    end
    progressBar(i);
end
   
% allflds = arrayfun(@(Gid) fieldnames( load( )), gids, 'un', 0)