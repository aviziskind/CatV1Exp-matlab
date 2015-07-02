function tmp_test

    S_grp = load('cellsGroups_grating_dOr_DB');

    S_data = load('indivCells_DB_grating_dOr');

    3;
    nCells = arrayfun(@(s) length(s.cellIds(s.cellIds>0)), S_grp.gratingGroups_dOr);
    
    
    ind5 = find(nCells == 5);
    allGids = [S_grp.gratingGroups_dOr(ind5).Gid];
    allCellIds = {S_grp.gratingGroups_dOr(ind5).cellIds};
    3;
    couldBeOK = zeros(1,length(ind5));
    for i = 1:length(allGids);        
        Gid = allGids(i);
        cellIds = allCellIds{i};
        cellIds = cellIds(cellIds>0);
        ok_tmp = zeros(1,length(cellIds));
        for ci = 1:length(cellIds)
            fn = getName('celldata', Gid, cellIds(ci));
            v = S_data.(fn);
            ok_tmp(ci) = v.OSP.stats.tuningStats.oriStats_si.cellOK;
        end

        couldBeOK(i) = isequal(ok_tmp, [0 1 1 0 1]);
        3;
        
    end
    3;

end
% Gid = 625 !!