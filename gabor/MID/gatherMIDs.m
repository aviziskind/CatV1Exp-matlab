function gatherMIDs(idx)

    addFittedGaborStats = 1;
    addSTAs = 0;
        flipMIDforPosCorrWithSTA = 0;
        onlyIncludeMIDsWithFits = 0;
        
    addMID_fit = 0;
        addSigPixels_gabor = 1;
        addSigPixels_2stdMID = 0;
    
    addMID_jacks = 1;
        addMID_jack_fits = 1;

%     timeWindow = curTimeWindow('');
    curResponseType('gainCorrected');
    allTimeWindows = {'best'};
%     allTimeWindows = {'best', [29, 62], [58, 91]};
    allTrialModes = {'all', 'odd', 'even'};
    
    nWindows = length(allTimeWindows);
    nTrialModes = length(allTrialModes);
    responseType = curResponseType('');
    
    S_cells = load('usedCells.mat');    
    getAllCells = 1;    
    if getAllCells         
        S.allGids = S_cells.allGids;
        S.allCellIds = S_cells.allCellIds;
    else
        S.allGids = S_cells.usedGids;
        S.allCellIds = S_cells.usedCellIds;    
    end
    
    nCells= length(S.allCellIds);
%     nCells = 4;
%     S = struct;
%     S.allMIDs = cell(nWindows,nCells);
%     S.allMIDs_odd = cell(nWindows,nCells);
%     S.allMIDs_even = cell(nWindows,nCells);
%     S.rsqr = zeros(1,nCells);
%     S.jackMeanCC = zeros(1,nCells);
%     S.gparams = cell(1,nCells);
%     S.MID_fit = cell(1,nCells);
%     S.above2std_mid = cell(1,nCells);
%     S.above2std_gabor = cell(1,nCells);
    
%     S.frac_above2std_mid = zeros(1,nCells);
%     S.frac_above2std_gabor = zeros(1,nCells);
    
    if addSTAs        
        S_indivCells = load( getFileName('indiv', 'movie_fg') );         
        S.allSTAs = cell(1,nCells);        
    end
    
    sig_pixels_gabor_flds = {}; %#ok<NASGU>
    if addSigPixels_gabor
        sig_pixels_gabor_flds = { 'above2std_gabor', [], 'frac_above2std_gabor', [] };
        if addMID_jack_fits
            sig_pixels_gabor_flds = {sig_pixels_gabor_flds{:}, 'above2std_gabor_jacks', [], 'frac_above2std_gabor_jacks', [] };
        end
    end
    
    sig_pixels_mid_flds = { };
    if addSigPixels_2stdMID 
        sig_pixels_mid_flds = { 'above2std_mid', [], 'frac_above2std_mid', [] };        
    end
    
    mid_jacks_flds = {}; %#ok<NASGU>
    if addMID_jacks
        mid_jacks_flds = {'MID_jacks', []};
        if addMID_jack_fits
            mid_jacks_flds = {mid_jacks_flds{:}, 'gparams_jacks', [], 'rsqr_jacks', [], 'MID_fit_jacks', []};
        end
    end
            
    
    emptyStruct = []; %struct('MID', [], 'jackMeanCC', [], 'gparams', [], 'rsqr', [], 'MID_fit', [], sig_pixels_gabor_flds{:}, sig_pixels_mid_flds{:}, mid_jacks_flds{:}, 't_calc', []);
    
    
    
    if nargin < 1
        idx = 1:nCells;
    end
    
    S_init = S;
    
    for wi = 1:nWindows
        
        timeWindow = allTimeWindows{wi};
        fprintf('\n\n ===================================== \n GATHERING ALL MIDS (timewindow = %s)\n =====================================\n\n', num2str(timeWindow));
        timeWindow_str = iff(strcmp(timeWindow, 'best'), '', sprintf('__%d_%d', timeWindow));           
%%        
    %     trialModes = {'all', 'odd', 'even'};
        S = S_init;
        progressBar('init-', nCells);
        for ci = 1:nCells
           
            cell_idx = idx(ci);
            Gid = S.allGids(cell_idx);
            cellId = S.allCellIds(cell_idx);
    %         Gid = 2022;
    %         cellId = 1;
            
            if addSTAs
                fn = getName('celldata', Gid, cellId);
                STA = single( S_indivCells.(fn).STAs.STA);
                S.allSTAs{cell_idx} = STA;
                3;
            end

            for tm_i = 1:nTrialModes
                trialMode = allTrialModes{tm_i};
                trialMode_str = iff(strcmp(trialMode, 'all'), '', ['_' trialMode]);
                
                mid_fileName = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode, responseType);        
                if ~exist(mid_fileName, 'file')
                    continue;
                end
                S_mid = load(mid_fileName);        
                MID = S_mid.MID;
                v_MID = S_mid.v_MID{1};
                fld_name = ['allMIDs' timeWindow_str trialMode_str];
                
                addMID_jack_fits_now = addMID_jack_fits && isequal(timeWindow, 'best'); % && any(strcmp(trialMode, {'odd', 'even'}));
                
                Si = emptyStruct;
                Si.MID = single(MID);                
                
                [nx, ny, nJacks] = size(v_MID);
                jack_MID_vec = reshape(v_MID, [nx*ny, nJacks]);

                [~, jack_ccs] = pearsonRm(jack_MID_vec);        
                meanJackCCs = mean(jack_ccs);
                Si.jackMeanCC = meanJackCCs;
                Si.t_calc = S_mid.t_calc;
                
                [gparams_jacks, rsqr_jacks, MID_fit_jacks] = deal(cell(1, nJacks));
                jackknifeIdx = [];
                if addFittedGaborStats
                    [gparams, rsqr, MID_fit] = mid_getCellGaborParams(Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx, [], addMID_fit);
                    
                    if addMID_jack_fits_now
                        for jack_idx = 1:nJacks
                            [gparams_jacks{jack_idx}, rsqr_jacks{jack_idx}, MID_fit_jacks{jack_idx}] = mid_getCellGaborParams(Gid, cellId, timeWindow, trialMode, responseType, jack_idx, [], addMID_fit);
                        end                        
                    end


                end
                if isempty(gparams)
                    if onlyIncludeMIDsWithFits
                        continue;
                    end
                    3;
%                     keyboard;
                end
                if addFittedGaborStats && ~isempty(gparams)
                    
                    if addSigPixels_2stdMID
                        MID_median = median(MID(:));
                        std_mid = std(trimmed(MID(:), 10));

                        above2std_mid = abs(MID-MID_median) > std_mid * 2;                    
                        frac_above2std_mid = nnz(above2std_mid) / numel(above2std_mid);
                        
                        Si.above2std_mid = above2std_mid;
                        Si.frac_above2std_mid = frac_above2std_mid;
                    end
                    
                    if addSigPixels_gabor
                        Si.above2std_gabor = sig_gabor_pixels(gparams, Gid);
                        Si.frac_above2std_gabor = nnz(Si.above2std_gabor) / numel(Si.above2std_gabor);
                        
                        if addMID_jack_fits_now
                            for jack_idx = 1:nJacks
                                 Si.above2std_gabor_jacks{jack_idx} = sig_gabor_pixels(gparams_jacks{jack_idx}, Gid);
                                 Si.frac_above2std_gabor_jacks{jack_idx} = nnz(Si.above2std_gabor) / numel(Si.above2std_gabor_jacks{jack_idx});
                            end
                        end                        
                    end

                else
                    [gparams,   rsqr, MID_fit] = deal([]);
%                     [,    above2std_gabor,    frac_above2std_gabor] = deal([]);
                end
                [Si.gparams, Si.rsqr, Si.MID_fit, Si.gparams_jacks, Si.rsqr_jacks, Si.MID_fit_jacks] = deal(...
                    gparams,    rsqr,    MID_fit, gparams_jacks,    rsqr_jacks,    MID_fit_jacks);                                   
%                 [Si.above2std_mid, Si.frac_above2std_mid, Si.above2std_gabor, Si.frac_above2std_gabor] = deal(...
%                     above2std_mid,    frac_above2std_mid,    above2std_gabor,    frac_above2std_gabor);                   
                    
                if addMID_jacks
                    assert(length(S_mid.v_MID) == 1);
                    Si.MID_jacks = single( S_mid.v_MID{1} );
                    
                end
                
                show = 0;
                if show
                    [xs, ys] = getStimulusXY(Gid);
                    figure(656); clf;
                    subplot(1,4,1); imagesc(xs, ys, MID); axis equal tight; title(sprintf('jack_{cc} = %.2f', meanJackCCs))
                    subplot(1,4,2); imagesc(xs, ys, MID_fit); axis equal tight; title(sprintf('rsqr = %.2f', rsqr))
                    subplot(1,4,3); imagesc(xs, ys, above2std_mid); axis equal  tight; hold on;  caxis([0 1]);
                    subplot(1,4,4); imagesc(xs, ys, above2std_gabor); axis equal  tight; hold on;  caxis([0 1]);
                end
                3;

                if addSTAs && flipMIDforPosCorrWithSTA
                    cc_sm = normDotProd(MID, STA);
                    if (cc_sm < 0) && flipMIDforPosCorrWithSTA

                        cc = zeros(1, 4);
                        for jj = 1:4
                            cc(jj) = normDotProd(STA, v_MID(:,:,jj));
                        end                        
                        
                        figure(67); clf;
                        subplot(2,4,1:2); imagesc(STA); axis square; colorbar
                        h(5) = subplot(2,4,3:4); imagesc(MID); axis square; colorbar;
                        for jj = 1:4
                            h(jj) = subplot(2,4,4+jj); imagesc(v_MID(:,:,jj)); title(sprintf('dot = %.2f', cc(jj)))
                        end
                        matchAxes('C', h)
                        3;
                    end
                end
                
                if isempty(emptyStruct)
                    emptyStruct = blankStruct(Si);
                end
                
                if ~isfield(S, fld_name)
                    S.(fld_name)(nCells) = emptyStruct;
                end
                S.(fld_name)(cell_idx) = Si;
                
                
            end
            progressBar(ci);
        end
        
       
        
        all_mids_filename = ['allMIDs' timeWindow_str '.mat'];

        save([CatV1Path all_mids_filename], '-struct', 'S', '-v6');
            %%
        mid_types = fieldnames(S); mid_types = mid_types(strncmp(mid_types, 'allMID', 6));
        for type_i = 1:length(mid_types);
            nMIDs(type_i) = nnz( arrayfun(@(cell_i) ~isempty(cell_i.MID), S.(mid_types{type_i})) );
            nMID_fits(type_i) = nnz( arrayfun(@(cell_i) ~isempty(cell_i.gparams), S.(mid_types{type_i}) ) ) ;
            all_t_calc =      arrayfun(@(cell_i) cell_i.t_calc, S.(mid_types{type_i}), 'un', 0);
            t_calc(type_i) =  sum([all_t_calc{:}]);
        end

        fprintf('Saved the following to file %s : \n', all_mids_filename);
        for i = 1:length(mid_types)        
            fprintf('  %s : %d MIDs, %d fits,   [ took %s]\n', mid_types{i}, nMIDs(i), nMID_fits(i), sec2hms(t_calc(i)) );
        end
        fprintf('Total time for calculation: %s \n', sec2hms( sum(t_calc) ));
        fprintf('Total time for calculation (divided by 8): %s \n', sec2hms( sum(t_calc)/8 ));
        
        
    end
   
    
        %%
    3;

%     nMIDs = cellfun(@(mid_type) nnz( arrayfun(@(st) ~isempty(st.MID), S.(mid_type))),  mid_types);
%     nMID_fits = cellfun(@(mid_type) nnz( arrayfun(@(st) ~isempty(st.t_calc), S.(mid_type))),  mid_types);
%     t_calc = cellfun(@mid_type)  arrayfun(@(st
    
    
    mid_getCellGaborParams('save');
%     cellfun(@(mid_type) fprintf('  %s : %d\n', mid_type, nnz(~cellfun(@isempty, S.(mid_type)))), mid_types, 'un', 0)
    %%
%     mid_types = fieldnames(S); mid_types = mid_types(strncmp(mid_types, 'allMID', 6));
    if 0  %this crashes matlab
        SS.a = S.allMIDs;
        mid_types = {'a'};

        %%
        cellfun(@(mid_type) fprintf('  %s : %d\n', mid_type, nnz(~cellfun(@isempty, SS.(mid_type)))), mid_types, 'un', 0)
    
        
%         cellfun(@(s) fprintf('  %s : %d\n', mid_type, nnz(~cellfun(@isempty, S.(mid_type)))), mid_types, 'un', 0);
    end    
    %%
    
    

    %%
%     S.a = zeros(1, 10); S.b = zeros(1, 10); S.c = zeros(1, 10);
%     fn = fieldnames(S);
%     cellfun(@(fld_name) fprintf('%s : %d \n', fld_name, nnz(~cellfun(@ isempty, S.(fld_name)))), fn);
    
    
    %%
    
    
    
    %%
%     cellfun(@(s) fprintf('  %s : %d\n', mid_types, nnz(~cellfun(@isempty, S.(mid_type))), mid_types, 'un', 0)
%     nAll_MIDs = ;
%     nOdd_MIDs = nnz(~cellfun(@isempty, S.allMIDs_odd));
%     nEven_MIDs = nnz(~cellfun(@isempty, S.allMIDs_even));

    
    
    
    
%     (%d all-trial MIDs, %d odd-trial MIDs, %d even-trial MIDs to file %s\n', nAll_MIDs, nOdd_MIDs, nEven_MIDs, fn);
      
    3;
%     
%     mid_dir = [experimentDBPath 'MID' filesep 'Cells' filesep];
%     mid_files = dir([mid_dir 'MID*.mat']);
%         
%     nCells = length(mid_files);
%     
%     allMIDs = cell(1,nCells);
%     allGids = zeros(1,nCells);
%     allCellIds = zeros(1,nCells);
%     for i = 1:nCells
%         gc = sscanf(mid_files(i).name, 'MID_Group_%d__cell_%d_');
%         
%         S = load([mid_dir mid_files(i).name]);        
%         allMIDs{i} = S.MID;        
%         allGids(i) = gc(1);
%         allCellIds(i) = gc(2);
%     end
%     save([CatV1Path 'allMIDs.mat'], 'allMIDs', 'allGids', 'allCellIds');

    
end