function explorePvalAreaWeights

    psthWindowPrefsFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowPrefsFile.mat'];
    if ~exist(psthWindowPrefsFile, 'file')
        error('psth window preferences file does not exist');
    end
    S = load(psthWindowPrefsFile);
    allPrefs = S.allPrefs;
%            allPrefs = struct('Gid', {allCells.Gid}, 'cellId', {allCells.cellId}, 'L_bin', [], 'R_bin', [], 'nStim', []);
    idx_havePref = cellfun(@(p) ~isempty(p), {allPrefs.L_bin});
    allPrefs = allPrefs(idx_havePref);
    pref_gids = [allPrefs.Gid];
    pref_cellIds = [allPrefs.cellId];

    psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData.mat'];    
    allPsthWindowData_S = load(psthWindowDataFile);
    name_fun = @(Gid, cellId, useStim) sprintf('Gid_%04d_cell_%d_%s', Gid, cellId, iff(useStim, 'N', '0'));
    
%     idx_havePSTH = fieldnames(allPsthWindowData_S)
%     cellfun(@(p) ~isempty(p), {allPrefs.L_bin});

    S = load('flashedGratingCells_all');
    allCells = S.allCells;
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];
    nCells = length(allCells);
    idx_ok = false(nCells,1);
    
    for i = 1:nCells
        if any( pref_gids == allGids(i) & pref_cellIds == allCellIds(i) ) &&  ...
            isfield(allPsthWindowData_S, name_fun(allGids(i), allCellIds(i), 0));
            idx_ok(i) = true;
        end
    end
    
    allCells = allCells(idx_ok);
    allGids = [allCells.Gid];
    allCellIds = [allCells.cellId];
    nCells = length(allCells);
    
    binW = 25/6;
    bins = [binW/2:binW:180];

    
    [Lbins_saved, Rbins_saved, Lbins_auto, Rbins_auto] = deal( zeros(1,nCells) );
    [tstats_n, areas_n] = deal( cell(1,nCells) );

    for i = 1:nCells
        idx = find( pref_gids == allGids(i) & pref_cellIds == allCellIds(i), 1);
        Lbins_saved(i) = allPrefs(idx).L_bin;
        Rbins_saved(i) = allPrefs(idx).R_bin;      
        
        [tmp, tstats_n{i}, areas_n{i}] = getPSTHwindowData(allGids(i), allCellIds(i), 1);
    end    
    figure(155); clf; h_mags = plot(0,0, 'bo'); h_mag_tit = title('');
    figure(156); clf; h_lrdist = plot(0,0, 'bo'); 
    

    
    function updatePlot(nStimMethod, nStim, W)
        
        % go to each cell. 
        
        for ci = 1:nCells
%             Gid = allGids(i);
%             cellId = allCellIds(i);
            
            if nStim == 0
%                 [pvals, tstats, areas] = getPSTHwindowData(Gid, cellId, 0);
            else
%                 [pvals_n, tstats_n, areas_n] = getPSTHwindowData(Gid, cellId, 1);
                tstats = tstats_n{ci}(:,:,nStim);
                areas  = areas_n{ci}(:,:,nStim);
            end                
                
%             end
            [Lbins_auto(ci), Rbins_auto(ci)] = getBestLbinRbin(bins, tstats, areas, W);            
        end                
        
        dists = normV([Lbins_auto-Lbins_saved;  Rbins_auto-Rbins_saved], 1);
        set(h_mags, 'xdata', 1:nCells, 'ydata', dists);
        set(h_mag_tit, 'string', sprintf('Mean = %.2f. Median = %.2f', mean(dists), median(dists)));
        set(h_lrdist, 'xdata', Lbins_auto-Lbins_saved, 'ydata', Rbins_auto-Rbins_saved);
    end

    updatePlot('', 10, .5)

    
    nStimMethods = {'top n', 'top x%'};
    
    args = { {'nStimMethod', nStimMethods}, ...
             {'nStim', [1:50], 10}, ...
             {'wgt', [0:.01:1], .5}, ...
             };
%              {'L_margin', [1:nBinsMax], L_bin}, ...
%              {'R_margin', [1:nBinsMax], R_bin}, ...

    manipulate(@updatePlot, args, 'FigId', 250);
    



end
