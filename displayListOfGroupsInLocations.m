
S_mcg = load('cellsGroups_movie_fg');
mcg = S_mcg.movieGroups_fg;

onlyShowDiffFrmLengths = false;

allGids = [mcg.Gid];
locDatas = [mcg.locationData];
locIds = [locDatas.LocId];
nGroups = length(locIds);
[uLocIds, locIdsList] = uniqueList(locIds);

st     = arrayfun(@(s) [s.stimType, '_' num2str(round(s.frameLength_ms))], mcg, 'un', 0);
st_ext = arrayfun(@(s) [s.stimType, '_' num2str(round(s.frameLength_ms)) '__' num2str(length(s.cellIds)) ':', rmds(num2str(s.nSpikes))], mcg, 'un', 0);

if onlyShowDiffFrmLengths
    idx_mult_u = find( cellfun(@(idx) length(unique(st(idx)))>1, locIdsList) );
else
    idx_mult_u = find( cellfun(@(idx) length(idx)>1, locIdsList ) );
end

id = 0;

displayByNumberOfGroups = false;

if displayByNumberOfGroups
    allNGs = cellfun(@length, {locIdsList{idx_mult_u}});
    uNGs = unique(allNGs);
    for ngi = 4;
        for id = find(allNGs == ngi)        
        %     id = id+1,
            idx = locIdsList{idx_mult_u(id)},
            Gids = allGids(idx),   
            st_ext(idx)
            allSame = length(unique(st(idx)))==1;
            displayCellOSP_PSTHs(Gids, id, allSame)
            3;
        end
        
    end
    
else
   
    for id = 1:length(idx_mult_u)
    %     id = id+1,
        idx = locIdsList{idx_mult_u(id)},
        Gids = allGids(idx),   
        st_ext(idx)
        allSame = length(unique(st(idx)))==1;
        displayCellOSP_PSTHs(Gids, id, allSame)
        3;
    end

end