function generateGroupSpikeData

%     compress = true;
    
    doAllCells = false;
    
    all_txt_load = iff(doAllCells, '_all', '');    
    
    S{1} = load(['flashedGratingCells' all_txt_load '.mat']);
    S{2} = load(['driftingGratingCells' all_txt_load '.mat']);
    
    fGids = unique([S{1}.allCells.Gid]);
    dGids = unique([S{2}.allCells.Gid]);

    allGids = unique([fGids(:); dGids(:)]);    
%     allGids = 4470;
    
    nGroups = length(allGids);

%     spikeCellIds = cell(1,nGroups);
%     tetrodeAmps = cell(1,nGroups);
%     spikeTimes = cell(1,nGroups);
    
%     groupSpikeData = struct('allGids', allGids, 'cellIds', {cell(1,nGroups)}, 'tetrodeAmps', {cell(1,nGroups)});
    
    fprintf('Gathering spike data for %d groups (%d flashed, %d drifting)...\n', nGroups, length(fGids), length(dGids));
    progressBar('init-', nGroups)
    for grp_i = 1:nGroups        
        progressBar;
        Gid = allGids(grp_i);
        
        groupName = ['group_' num2str(Gid)];
        spk = getSpikes(allGids(grp_i), [], [], 1);
        
        spikeTimes_tk = int32( spk(:,1) );
        spkCellIds = int8( spk(:,2) );
        tetrodeAmps = int16 ( spk(:,3:end) );
        
        allGroupData.(groupName) = struct('spikeTimes_tk', spikeTimes_tk, 'spkCellIds', spkCellIds, 'tetrodeAmps', tetrodeAmps);
        
%         spikeCellIds{grp_i} = spkCellIds;        
%         tetrodeAmps{grp_i} = tet;
        
    end
    
%     all_txt_save = iff(doAllCells, '_all', '');        
    if length(allGids) == 1
        all_txt_save = '1'; 
    else
        all_txt_save = iff(doAllCells, '_all', '');        
    end
        
    fname = [CatV1Path 'groupSpikes' all_txt_save '.mat'];        
    allGroupData = orderfields( allGroupData ); %#ok<NASGU>
    save(fname, '-struct', 'allGroupData');

    fprintf('Saved to %s\n', fname);
    
end