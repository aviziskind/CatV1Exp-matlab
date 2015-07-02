function test_MID_sign_bias


    S = load('flashedGratingCells_GLFcuw8_phase.mat');
    allCells = S.allCells;
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];
    N = length(allCells);
    [uGids, gid_idx] = uniqueList(allGids);
    
    mid_maxes = nan(1, N);
    all_osis = nan(1, N);
    jackCCs = nan(1,N);    
    rsqrs = nan(1,N);    
    
%     for gi = 1:length(uGids);
%         for ci = 1:length(gid_idx{gi})
%             j = gid_idx{gi}(ci);
                        
%             if allCellIds(j) > 0
    for j = 1:N
        if allCellIds(j) > 0
                
            mid = allCells(j).MIDdata.MID;
            mid_abs_max_idx = indmax(abs(mid(:)));
            mid_max = mid(mid_abs_max_idx);                
            mid_maxes(j) = mid_max;            
            
            jackCCs(j) = allCells(j).MIDdata.jackCC;
            rsqrs(j) = allCells(j).MIDdata.rsqr;
            
            osi_j = allCells(j).stats.tuningStats.oriStats_si.OSI;
            all_osis(j) = osi_j;
            
            jackCCs = nan(1,N);            
        end                            
    end
        
    idx_ok = ~isnan(all_osis) & ~isnan(mid_maxes);
    gids = allGids(idx_ok);
    all_osis = all_osis(idx_ok);
    mid_maxes = mid_maxes(idx_ok);
    jackCCs = jackCCs(idx_ok);
    rsqrs = rsqrs(idx_ok);
    
    [uGids, gid_idx] = uniqueList(allGids);
    
    3;






end