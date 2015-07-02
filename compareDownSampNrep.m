function compareDownSampNrep
%     dirName = [CatV1Path 'gabor\Cells\'];
    dirName = [experimentDBPath 'MID' filesep 'Cells' filesep];
    
    cd(dirName);
    fnames = dir([dirName '*.mat']);
    
    [allGids, allCellIds, allNReps, allDsamp] = deal(zeros(1, length(fnames)));
    for i = 1:length(fnames)       
        if isempty(strfind(fnames(i).name, 'odd')) && isempty(strfind(fnames(i).name, 'even')) && isempty(strfind(fnames(i).name, '-')) 
            [allGids(i), allCellIds(i), allNReps(i), allDsamp(i) ] = dealV(sscanf(fnames(i).name, 'MID_Group_%d__cell_%d_%drep_d%d')); %#ok<AGROW>
        end
    end
    uGC = unique([allGids(:), allCellIds(:)], 'rows');
    if isequal (uGC(1,:) , [0, 0] )
        uGC = uGC(2:end,:);
    end
    allGids = uGC(:,1); allCellIds = uGC(:,2);
    nCells = length(allCellIds);
    
    doMID_reps = {1, 2, 5, 'all'}; nnreps = length(doMID_reps);
    downSampFactors = [1];   nnsamp = length(downSampFactors);
  
    showSTA = 0;
    
    allFiles_present = zeros(1,nCells);
    
    for i = 1:nCells
        Gid = allGids(i);
        cellId = allCellIds(i);
           
        allFileNames = getAllFileNames(Gid, cellId, downSampFactors, doMID_reps);
        filesPresent = cellfun(@(s) exist(s, 'file'), allFileNames);
       
        if all(filesPresent)
            3;
        end
            
        
        allFiles_present(i) = all(filesPresent);
    end        
        
    idx_cellsWithAll = find(allFiles_present);
    nCellsWithAll = length(idx_cellsWithAll);
    for j = 1:nCellsWithAll
        i = idx_cellsWithAll(j);
        Gid = allGids(i);
        cellId = allCellIds(i);
        
        if showSTA
            s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
        end
                
        figure(i); clf;
        sub_i = 1;
        gt = getGratingStimType(Gid);
        nTrials = gt.nTrials;
        IMs = cell(nnsamp, nnreps);
        ts = zeros(nnsamp, nnreps);
        nreps = zeros(nnsamp, nnreps);
        diffs_m = zeros(nnsamp, nnreps);
        diffs_s = zeros(nnsamp, nnreps);

        allFileNames = getAllFileNames(Gid, cellId, downSampFactors, doMID_reps);
        %%
        allMeanCCs = zeros(nnsamp, nnreps);
        for di = 1:nnsamp
            for ri = 1:nnreps
                nrep = doMID_reps{ri};
                if ~isnumeric(nrep),
                    nrep = nTrials;
                end
                nrep_str = sprintf('%d reps', nrep);
                nreps(di, ri) = nrep;

                S = load(allFileNames{ri, di});
                subplotGap(nnsamp+showSTA,nnreps, di, ri);
                c = S.MID;
                meanCC = getMeanMIDcc(S.v_MID{1});
                allMeanCCs(ri) = meanCC;
                
                IMs{di, ri} = c;
                t_tot = sum(S.t_calc);
                ts(di, ri) = t_tot;
                t_str = sec2hms( t_tot );
                imagesc(c);
                axis square tight;
                set(gca, 'xtick', [], 'ytick', []);
                colorbar;
                caxis([-1, 1]*max(abs(c(:))));
%                     drawnow;


                if sub_i == 1
                    ylabel(sprintf('Gid = %d, cell %d', Gid, cellId));
                end
                title(sprintf('d = %d,  %s', downSampFactors(di), nrep_str));
                meanCC_str = sprintf('cc = %.3f', meanCC);
                xlabel({t_str, meanCC_str});

                sub_i = sub_i+1;
                3;
            end
            for ri = 1:nnreps
                D = abs( (IMs{di, ri}-IMs{di, end}) ) / max(abs(IMs{di, end} )) ;
                diffs_m(di, ri) = mean( D(:) );
                diffs_s(di, ri) = std( D(:) );
            end

        end
        if showSTA
            subplotGap(nnsamp+1,nnreps, sub_i);
            STA = s.STAs.STA; STA = STA/max(abs(STA(:)));
            imagesc(STA);
            c = caxis; caxis([-1, 1]*max(abs(c)));
            axis square; colorbar;

            subplot(nnsamp+1,nnreps, sub_i+1);
            imageOSP(s.OSP, 'mean:ph');

        end
%%
        figure(100); clf;
%         subplot(2,1,1); plot(nreps', ts', 'o-'); ylabel('time')
% %         legend('d2', 'd4', 'location', 'NW');
        x = [1:nnreps];
%         subplot(2,1,2); errorbar(x', diffs_m', diffs_s'); ylabel('diff');
            plot(x, allMeanCCs, 'o-'); ylim([0, 1]);
%         legend('d2', 'd4');
        3;

            3;


        
    end

end


function allFileNames = getAllFileNames(Gid, cellId, downSampFactors, doMID_reps)

    nnreps = length(doMID_reps);
    nnsamp = length(downSampFactors);
    
    
    allFileNames = cell(nnreps, nnsamp);
        
    timeWindow = 'best';
    trialMode = 'all';
    for di = 1:nnsamp
        for ri = 1:nnreps
            nrep = doMID_reps{ri};
            if isnumeric(nrep), 
                nrep = sprintf('%drep', nrep);
            end
            
            fn = getName('MID_file', Gid, cellId, downSampFactors(di), nrep, timeWindow, trialMode);
            [~, filename, ext] = fileparts(fn);
            allFileNames{ri, di} = [filename  ext];
        end
    end             
end