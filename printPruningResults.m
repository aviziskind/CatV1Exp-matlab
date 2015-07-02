function printPruningResults(splitArg)

    Gids = ratDatafilesAvailable(1);

    pruningFeatures = curPruningFeatures('');
    nGids = length(Gids);
%     nGids = 10;
            
    showSortingStats = 0;
    showPruningStats = 1;

    showWithAutoRefrPeriods = 1;
    showWithICRefrPeriods = 1;
    
    doAutoPruning = 1;
%     doCrossPruning_N = 0;
    doCrossPruning_D = 0;
    doMergePruning = 0;
    doPruneMergePruning = 0;

    curJustICclusters(1);
    
    if nargin == 1
        assert(ibetween(splitArg, 1, 8));
        pruningType_id_active = mod(splitArg-1, 4)+1;        
        allPruningIds = num2cell(false(1,4)); allPruningIds{pruningType_id_active} = 1;
        
        refrType_id_active = double(splitArg>4) + 1;
        allRefrIds = num2cell(false(1,2)); allRefrIds{refrType_id_active} = 1;
        [doAutoPruning, doCrossPruning, doMergePruning, doPruneMergePruning] = allPruningIds{:};
        [showWithAutoRefrPeriods, showWithICRefrPeriods] = allRefrIds{:};
        
        [doCrossPruning_N, doCrossPruning_D] = deal( doCrossPruning );
    end
    
%     pruningMode = 'autoPruning';
%     pruningMode = 'crossPruning';
%     pruningMode = 'mergedCrossPruning';
    
    cellSorting = 'clusters';
    if showSortingStats
        [fracFalseNeg_detect, fracFalseNeg_sorting, fracFalsePos_sorting, ...
         nCorrect_detect, nFalseNeg_detect, nFalseNeg_sorting, nFalsePos_sorting] = deal(...
            zeros(1, nGids));
        
        progressBar('init-', nGids);
        for i = 1:nGids
            [stats, idxs] = identifyCellfromIC(Gids(i), cellSorting, 0);                

            if isfield(stats, 'fracFalseNegDetections')
                3;
                
                fracFalseNeg_detect(i) = stats.fracFalseNegDetections;
                fracCorrect_detect(i) = stats.fracCorrectDetections;
%                 fracFalseNeg_sorting(i) = stats.fracFalseNegSorting;
%                 fracFalsePos_sorting(i) = stats.fracFalsePosSorting;
                
%                 nCorrect_detect(i) = stats.nCorrect_sorting;
%                 nFalseNeg_detect(i) = stats.nFalseNeg_detection;
%                 nFalseNeg_sorting(i) = stats.nFalseNeg_sorting;
%                 nFalsePos_sorting(i) = stats.nFalsePos_sorting;
                3;
            else
                3;
                fprintf('[Gid = %d]', Gids(i))
            end
                
%             progressBar;
        end

        %%
        figure(900); plot(nCorrect_detect, nFalseNeg_detect, '.'); xlabel('nCorrect'); ylabel('nMissed')
        3;
        figure(901); hist(fracFalseNeg_detect, 50);  xlabel({'Fraction of False-negative (detection) IC spikes', 'IC spike, but no EC spike detected'})
        figure(902); hist(fracFalseNeg_sorting, 50); xlabel({'Fraction of False-negative (sorting) IC spikes', 'IC spike, but EC assigned a to different cell'})
        figure(903); hist(fracFalsePos_sorting, 50); xlabel({'Fraction of False-positive (sorting) IC spikes', 'EC assigned to cell, but no IC spike'})
        3;
    end

    sortingFeatureSets_representative1 = { ...
                'Neg', 'PCAcuw4', 'PCAcuw8', 'PCAcuw12', 'PCAcuw16', ...
                'PCAcur4', 'PCAcur8', 'PCAcur12', 'PCAcur16', ...
                'PCAsur1', 'PCAsur2', 'PCAsur3',  'PCAsur4'};

    sortingFeatureSets_representative48 = { ...
         'Neg', 'PCAcuw4', 'PCAcuw8', ...
                'PCAcur4', 'PCAcur8', ...
                'PCAsur1', 'PCAsur2'};
            
    sortingFeatureSets_representative4 = { ...
         'Neg', 'PCAcuw4', 'PCAcur4', 'PCAsur1'};
            
%     features_smoothN = { ...   removed 'Neg' from beginning
%                   'PCAcuw1', 'PCAcuw2', 'PCAcuw3', 'PCAcuw4', 'PCAcuw5', 'PCAcuw6', 'PCAcuw8', 'PCAcuw10', 'PCAcuw12', 'PCAcuw16', ...
%                   'PCAcur1', 'PCAcur2', 'PCAcur3', 'PCAcur4', 'PCAcur5', 'PCAcur6', 'PCAcur8', 'PCAcur10', 'PCAcur12', 'PCAcur16', ...
%                   'GLFcuw1', 'GLFcuw2', 'GLFcuw3', 'GLFcuw4', 'GLFcuw5', 'GLFcuw6', 'GLFcuw8', 'GLFcuw10', 'GLFcuw12', 'GLFcuw16', ...
%                   'GLFcur1', 'GLFcur2', 'GLFcur3', 'GLFcur4', 'GLFcur5', 'GLFcur6', 'GLFcur8', 'GLFcur10', 'GLFcur12', 'GLFcur16'};    

        sortingFeatureSets_representative2_short = { ...
         'Neg', 'PCAcuw4', 'PCAcuw6', 'PCAcuw8', 'PCAcuw10', 'PCAcuw12', ...
                'PCAcur4', 'PCAcur6', 'PCAcur8', 'PCAcur10', 'PCAcur12', ...
        };

        sortingFeatureSets_full = { ...
         'Neg',  'PCAcuw2', 'PCAcuw3', 'PCAcuw4', 'PCAcuw5', 'PCAcuw6', 'PCAcuw8', 'PCAcuw12', ...
                 'PCAcur2', 'PCAcur3', 'PCAcur4', 'PCAcur5', 'PCAcur6', 'PCAcur8', 'PCAcur12' ...
                 'GLFcuw2', 'GLFcuw3', 'GLFcuw4', 'GLFcuw5', 'GLFcuw6', 'GLFcuw8', 'GLFcuw12', ...
                 'GLFcur2', 'GLFcur3', 'GLFcur4', 'GLFcur5', 'GLFcur6', 'GLFcur8', 'GLFcur12' ...
        };
    
        features_smoothN = { ...
              'PCAsur1', 'PCAsur2', 'PCAsur3', ...
              'GLFsur1', 'GLFsur2', 'GLFsur3', ...            
...               'PCAsur1', 'PCAsur2', 'PCAsur3', ...
...               'PCAsuw1', 'PCAsuw2', 'PCAsuw3', ...
...               'GLFsur1', 'GLFsur2', 'GLFsur3', ...
...               'GLFsuw1', 'GLFsuw2', 'GLFsuw3', ...
              'PCAcuw1', 'PCAcuw2', 'PCAcuw3', 'PCAcuw4', 'PCAcuw5', 'PCAcuw6', 'PCAcuw8', 'PCAcuw10', 'PCAcuw12', 'PCAcuw16', ...
              'PCAcur1', 'PCAcur2', 'PCAcur3', 'PCAcur4', 'PCAcur5', 'PCAcur6', 'PCAcur8', 'PCAcur10', 'PCAcur12', 'PCAcur16', ...
              'GLFcuw1', 'GLFcuw2', 'GLFcuw3', 'GLFcuw4', 'GLFcuw5', 'GLFcuw6', 'GLFcuw8', 'GLFcuw10', 'GLFcuw12', 'GLFcuw16', ...
              'GLFcur1', 'GLFcur2', 'GLFcur3', 'GLFcur4', 'GLFcur5', 'GLFcur6', 'GLFcur8', 'GLFcur10', 'GLFcur12', 'GLFcur16', ...
        };
              
              
              %         {{'Egy', 'PCAsnr1'}},     
%     sortingFeatureSets = { ...
%         'GLFcuw12'    'GLFcuw16'    'GLFcuw3'    'GLFcuw32' ... 
%         'GLFcuw4'    'GLFcuw48'    'GLFcuw5'    'GLFcuw6'    'GLFcuw8'    'GLFsnr3'    'Neg' ...
%         'PCAcur1'    'PCAcur10'    'PCAcur12'    'PCAcur2'    'PCAcur3'    'PCAcur4' ... 
%         'PCAcur5'    'PCAcur6'    'PCAcur7'    'PCAcur8'    'PCAcuw1'    'PCAcuw12' ... 
%         'PCAcuw2'    'PCAcuw3'    'PCAcuw4'    'PCAcuw5'    'PCAcuw6'    'PCAcuw7' ...
%         'PCAcuw8'    'PCAsnr3'    'PCAsur1'    'PCAsur2'    'PCAsur3'    'PCAsur4' ...
%         'PCAsuw3'};

%     sortingFeatureSets = {'PCAcuw8'}; 
        
%     sortingFeatureSets = {'PCAcur1'    'PCAcur10'    'PCAcur12'    'PCAcur2'    'PCAcur3'    'PCAcur4' ... 
%         'PCAcur5'    'PCAcur6'    'PCAcur7'    'PCAcur8'    'PCAcuw1'    'PCAcuw12'};

%     
    sortingFeatureSets = { 'Neg',  'PCAcur4', 'GLFcuw4', 'GLFcuw8', 'GLFcuw12', 'GLFcuw16'  ...
         ...'GLFcuw32' ... 
         ...'GLFcuw48'    'GLFsnr3'    ...
         'PCAcur12'   ... 
         'PCAcur8'    'PCAcuw12' ...         
         'PCAcuw8'    ...'PCAsnr3'   'PCAsur3'  ...
         ...'PCAsuw3'
         };
%         sortingFeatureSets = {'GLFcuw8'};

        sortingFeatureSets_representative = { ...
        'GLFcuw4'    'GLFcuw8', 'GLFcuw12'    'GLFcuw16'   'GLFcuw32', 'GLFcuw48'  ... 
        'GLFsnr3'    'Neg', ...
         'PCAcur4', 'PCAcur8', 'PCAcur12', ...
         'PCAcuw4', 'PCAcuw8', 'PCAcuw12' ...         
         'PCAsnr3',  'PCAsur3', 'PCAsuw3', ...
        };

        sortingFeatureSets_representative2_have = { ...
         'Neg', 'PCAcuw3', 'PCAcuw4', 'PCAcuw5', 'PCAcuw6', 'PCAcuw8', 'PCAcuw10', 'PCAcuw12', 'PCAcuw16', ...
                'PCAcur3', 'PCAcur4', 'PCAcur5', 'PCAcur6', 'PCAcur8', 'PCAcur10', 'PCAcur12', 'PCAcur16' ...
        };

    
        sortingFeatureSets_representative5 = { ...
         'Neg', 'PCAcuw4', 'PCAcuw8', 'PCAcuw12', 'PCAcuw16', ...
                'PCAcur4', 'PCAcur8', 'PCAcur12', 'PCAcur16', ...
                'PCAsur1', 'PCAsur2', 'PCAsur3', 'PCAsur4'};


%     sortingFeatureSets = {'Neg', 'GLFcuw4', 'GLFcuw8', 'PCAcuw4', 'PCAcuw8', 'GLFcur4', 'PCAcur4'};
%     sortingFeatureSets = {'Neg', 'GLFcuw4', 'GLFcuw8', 'PCAcuw4'};
%     sortingFeatureSets = {'GLFcuw4', 'PCAcuw4', 'Neg', 'GLFcur4', 'PCAcur4'};
    
%     sortingFeatureSets = sortingFeatureSets_representative;
    
%     sortingFeatureSets = { 'GLFcuw4', 'GLFcuw8'};
    sortingFeatureSets = { 'Neg', 'PCAcuw4', 'PCAcuw8'};
    
%     pruningFeatureSets = { 'Neg', 'PCAcur4', 'GLFcuw4', 'GLFcuw8', 'GLFcuw12',    'GLFcuw16' };   
    pruningFeatureSets = { 'PCAcuw4', 'PCAsur1', 'PCAcur4'}; %, 'PCAcur4', 'GLFcuw4', 'GLFcuw8', 'GLFcuw12',    'GLFcuw16' };   
    
%     pruningFeatureSets = { 'Neg', 'GLFcuw4', 'GLFcuw8', 'PCAcuw4', 'PCAcuw8'};    
%     pruningFeatureSets = { 'Neg', 'GLFcuw4', 'GLFcuw8', 'PCAcuw4', 'PCAcuw8'};    
    prune_availNow = {'Neg', 'PCAcuw4'};
    
    if strcmp(getenv('computername'), 'AVI-PC') 
    
%         sortingFeatureSets = sortingFeatureSets_representative2_short;
        sortingFeatureSets = sortingFeatureSets_full;
        pruningFeatureSets = features_smoothN;
%         sortingFeatureSets = sortingFeatureSets_representative5;
%         pruningFeatureSets = prune_availNow; %sortingFeatureSets_representative1;
    
    else
        sortingFeatureSets = sortingFeatureSets_representative2_have;
        pruningFeatureSets = features_smoothN;
        
    end
        
    
%    pruningFeatureSets = { ...
%         'Neg', 'GLFcuw4', 'GLFcuw8', 'GLFcuw12',    'GLFcuw16', ...
%         'GLFsnr3', ...
%          'PCAcur4', 'PCAcur8', 'PCAcur12', ...
%          'PCAcuw4', 'PCAcuw8', 'PCAcuw12' ...         
%          'PCAsnr3',  'PCAsur3', 'PCAsuw3', ...
%         };    
    
    doCheckForFiles = 0;
    if doCheckForFiles
        check_opt = struct('showWithAutoRefrPeriods', showWithAutoRefrPeriods, 'showWithICRefrPeriods', showWithICRefrPeriods);    
        avail_tf = checkFilesAvailable(sortingFeatures, pruningFeatures, Gids, pruningMode, check_opt);
    end
    
%     if any(~hasData)
%         fprintf('Sorting features excluded (data not complete) : %s \n', cellstr2csslist(sortingFeatureSets(~hasData)));
%     end
    
%     sortingFeatureSets = sortingFeatureSets(hasData);
            
    overlapPlots = 0;

    nSortingSets = length(sortingFeatureSets);
    nPruningSets = length(pruningFeatureSets);
            
    refrPeriodTypes = [showWithAutoRefrPeriods, showWithICRefrPeriods];
    
    if showPruningStats
    
        allStats_SPG_refr_C = cell(1, 2);
        refrPeriodTypes_idxs = find(refrPeriodTypes);

        file_opt.redo = 0;
        file_opt.redoOldFiles = 1;
        file_opt.redoBeforeDate = 735235.010349; %%;  % sprintf('%.6f', now)
        file_opt.redoPrimaryDatafile = 0;
        
        % Step 1: Calculate / Load data for all sorting features, all pruning features, all Groups.
        [allStats_SPG_auto_C, allStats_SPG_crossD_C, allStats_SPG_crossN_C, allStats_SPG_merge_C, allStats_SPG_pruneMerge_C] = deal(cell(1, 2));
        
        
        for refr_type_i = refrPeriodTypes_idxs                                    
            useICcellRefrPeriods = (refr_type_i == 2);                
            if doAutoPruning
                allStats_SPG_auto_C{refr_type_i} = loadAllPruningGroupStats(Gids, sortingFeatureSets, pruningFeatureSets, 'autoPruning', useICcellRefrPeriods, file_opt);                
            end
            if doCrossPruning_D
                allStats_SPG_crossD_C{refr_type_i} = loadAllPruningGroupStats(Gids, sortingFeatureSets, pruningFeatureSets, 'crossPruning', useICcellRefrPeriods, file_opt, 'D');
            end
%             if doCrossPruning_N
%                 allStats_SPG_crossN_C{refr_type_i} = loadAllPruningGroupStats(Gids, sortingFeatureSets, pruningFeatureSets, 'crossPruning', useICcellRefrPeriods, file_opt, 'N');
%             end            
            if doMergePruning
                allStats_SPG_merge_C{refr_type_i} = loadAllPruningGroupStats(Gids, sortingFeatureSets, pruningFeatureSets, 'mergePruning', useICcellRefrPeriods, file_opt);
            end
            if doPruneMergePruning
                allStats_SPG_pruneMerge_C{refr_type_i} = loadAllPruningGroupStats(Gids, sortingFeatureSets, pruningFeatureSets, 'pruneMergePruning', useICcellRefrPeriods, file_opt);
            end            
            
        end
        
        % Step 2: Decide if want to limit only one set for each group, or if can have multiple sorting features per group.
        %%
%         maximizeWhenChoosingGid = 'nPruned';
%             maximizeWhenChoosingGid = 'nClusts';
        if nargin == 1
            return;
        end
        displayGidStats = 0;
        if displayGidStats
            printGroupStats(allStats_SPG_auto_C{1})
        end
        
%         return;
        maximizeWhenChoosingGid = 'all';
        
        plot_opt.nPrunedTotal_min = 2;
        plot_opt.nPrunedCluster_min = 2;
        plot_opt.minFracIC = 0.5;
        plot_opt.figId = 550;
        plot_opt.maximizeWhenChoosingGid = maximizeWhenChoosingGid;
        
        
                plot_opt.maximizeWhenChoosingGid = 'nClusts';
%         compareFracPruned_vs_fracFalsePos(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, plot_opt);
        doFigs = 3;
        
        %%
        if any(doFigs == 3)
            allPruningTypes = {'autoPruning', 'crossPruning', 'mergePruning', 'pruneMergePruning'}; % figures: 3, 5
            for i = 1:1; %length(allPruningTypes)
                %%
%                 for j = 1:length(sortingFeatureSets)
                    pruningType = allPruningTypes{1};

                    plot_opt.sortFet = 'PCAcur5'; 
                    plot_opt.pruneFet = 'GLFcuw3';
                    [plot_opt1, plot_opt2] = deal(plot_opt);        
                    plot_opt1.figId = 100 + i*100; 
                    plot_opt2.figId = 110 + i*100; 
%                     [medP_auto(j), meanE_auto(j)] = plot_frTP_vs_frFP(allStats_SPG_auto_C{1}, pruningType, 0, sortingFeatureSets, pruningFeatureSets, plot_opt);
                    plot_frTP_vs_frFP(allStats_SPG_auto_C{1}, pruningType, 0, sortingFeatureSets, pruningFeatureSets, plot_opt);
                    plot_frTP_vs_frFP(allStats_SPG_auto_C{2}, pruningType, 1, sortingFeatureSets, pruningFeatureSets, plot_opt2);
%                 end
                3;
            end
        end
        3;
                

%         plot_frTP_vs_frFP(allStats_SPG_cross_C{1}, 'crossPruning', 0, sortingFeatureSets, pruningFeatureSets, plot_opt1)
%         plot_frTP_vs_frFP(allStats_SPG_cross_C{2}, 'crossPruning', 1, sortingFeatureSets, pruningFeatureSets, plot_opt2)
%         plot_frTP_vs_frFP(allStats_SPG_merge_C{1}, 'mergePruning', 0, sortingFeatureSets, pruningFeatureSets, plot_opt1)
%         plot_frTP_vs_frFP(allStats_SPG_merge_C{2}, 'mergePruning', 1, sortingFeatureSets, pruningFeatureSets, plot_opt2)

        
%         showAverageFracRemoved_vs_refrPeriod(allStats_SPG_auto_C{1}, struct('figId', 800));
        
%         showPerformance_Vs_Sorting_Pruning_Features(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt);
        if any(doFigs == 5)
              plot_opt.figId = 100; 
              showPerformance_Vs_variables(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt);
%               
%               plot_opt.figId = 200; 
%               showPerformance_Vs_X(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt, 'PCA_vs_GLF');
%               
%               plot_opt.figId = 300; 
%               showPerformance_Vs_X(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt, 'ccw_vs_raw');
%               
%               plot_opt.figId = 400; 
%               showPerformance_Vs_X(allStats_SPG_auto_C{1}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt, 'cat_vs_sep');
        
        end
%         showPerformance_Vs_Sorting_Pruning_Features
        
%         plot_opt.figId = 551;
%         showmeanRatio_Vs_Sorting_Pruning_Features(allStats_SPG_auto_C{2}, sortingFeatureSets, pruningFeatureSets, 'autoPruning', plot_opt);
        
        3;
%         compareTwoPruningMethods(allStats_SPG_cross_C{1}, allStats_SPG_merge_C{1}, sortingFeatureSets, pruningFeatureSets, plot_opt);

%         
%         plot_opt.maximizeWhenChoosingGid = 'nClusts';
%         compareTwoPruningMethods(allStats_SPG_cross_C{1}, allStats_SPG_pruneMerge_C{1}, 'crossPruning', 'mergePruning', plot_opt);
%         compareTwoPruningMethods(allStats_SPG_merge_C{1}, allStats_SPG_pruneMerge_C{1}, 'Cross Pruning', 'Merge Pruning', plot_opt);

        compareCrossPruning_D_vs_N = 0;
        if compareCrossPruning_D_vs_N
            compareTwoPruningMethods(allStats_SPG_crossD_C{1}, allStats_SPG_crossN_C{1}, 'crossPruning_D', 'crossPruning_N', plot_opt);
        end

%         compareMergePruning_and_PruneMergePrune(allStats_SPG_merge_C{1}, allStats_SPG_pruneMerge_C{1}, plot_opt);
        
        3;
%         comparePruningMethods(allStats_SPG_cross_C{1}, allStats_SPG_pruneMerge_C{1}, 'Cross Pruning', 'Merge Pruning', plot_opt);
        
        3;
        
        % *** Sorting False negative: IC spike, but EC was assigned a different clustId/cellId
        % *** Sorting False positive: spike was assigned to IC cell (by spike-sorting algorithm), but no IC spike
        
%         if overlapPlots
%             plotOffset = 0; % first do ideal = 0, then ideal = 1;            
%             makeNewPlot = ~useICcellRefrPeriods;
%         else
%             
%             makeNewPlot = 1;
%         end
% 
%         if makeNewPlot
%             markerStyle = 'b.';
%         else
%             markerStyle = 'go';
%         end
%         
%         cols = jet(nSortingSets);       
%                 
%         
%         for refr_type_i = find(refrPeriodTypes)
%             
%         end
        
    3;

    end

end



function avail_tf = checkFilesAvailable(sortingFeatures, pruningFeatures, Gids, pruningMode, opt)
    
    nSortingSets = length(sortingFeatures);
    nPruningSets = length(pruningFeatures);
    
    avail_tf = zeros(nSortingSets, nPruningSets);
    
    %     sortingFeatureSets = sortingFeatureSets(1:4);
    switch pruningMode
        case 'autoPruning',  fileType = 'clusterPrunings';
        case 'crossPruning', fileType = 'clusterData2';
        case 'mergePruning', fileType = 'clusterData3';
    end
    %     fileType_auto = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'clusterPrunings', 'clusterData2'});
    %     fileType_IC = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'clusterPrunings_ideal', 'clusterData2'});
    
    haveFile = @(filename) exist(filename, 'file');
    hasData = false(nSortingSets, nPruningSets);
    for i = 1:nSortingSets
        for j = 1:nPruningSets
            file_opt_auto = struct('sortingFeatures', sortingFeatures{i}, 'pruningFeatures', pruningFeatures{j}, 'ideal', 0);
            file_opt_IC   = struct('sortingFeatures', sortingFeatures{i}, 'pruningFeatures', pruningFeatures{j}, 'ideal', 1);
            %         clustPruning_fileType = iff(useIdealClusterPrunings, '', 'clusterPrunings');
            
            fileName_auto = getFileName(fileType, Gids(1), [], file_opt_auto);
            allFileNames_auto = arrayfun(@(G) strrep(fileName_auto, num2str(Gids(1)), num2str(G)), Gids, 'un', 0);
            haveFiles_auto = cellfun(haveFile, allFileNames_auto);
            
            fileName_IC = getFileName(fileType, Gids(1), [], file_opt_IC);
            allFileNames_IC = arrayfun(@(G) strrep(fileName_IC, num2str(Gids(1)), num2str(G)), Gids, 'un', 0);
            haveFiles_IC = cellfun(haveFile, allFileNames_IC);
            
            hasData(i,j) = (all(haveFiles_auto) || ~opt.showWithAutoRefrPeriods) && ...
                           (all(haveFiles_IC)   || ~opt.showWithICRefrPeriods);
        end
    end
end
    


function showPerformance_Vs_Sorting_Pruning_Features(allStats_SPG_C, sortingFeatureSets, pruningFeatureSets, pruningMode, opt)

    anyData_all = ~cellfun(@isempty, allStats_SPG_C);
    idx_gids_withData = find(any(any(anyData_all, 1),2));
    nGids_usable = length(idx_gids_withData);

    allStats_SPG_C = allStats_SPG_C(:, :, idx_gids_withData);
    
%     anyData_all = anyData_all(:,:,idx_gids_withData);

%     allStats_used_C = cell(1,nGids_usable);
%     idx_use = zeros(1, nGids_usable);

%             nPrunedTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpikesPruned'});
%             nSpkTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpkInClusts_Tot'});

    %%
    
    nSortingSets = length(sortingFeatureSets);
    nPruningSets = length(pruningFeatureSets);
    
    meanRatios = nan(nSortingSets, nPruningSets);
    ismeanRatioInf = true(nSortingSets, nPruningSets);
    median_falsePosRate_changes = zeros(nSortingSets, nPruningSets);
    
    for sort_fet_idx = 1:nSortingSets
        for prune_fet_idx = 1:nPruningSets
        
            allStats_C = [allStats_SPG_C{sort_fet_idx, prune_fet_idx, :}];    
            if isempty(allStats_C)
                continue;
            end
            S = getResultsStruct(allStats_C, opt, pruningMode);
            allP = S.falsePos_fracRemoved ./ S.correctICspikes_fracRemoved;
            meanRatio = median(allP);
            
            ismeanRatioInf = isinf(meanRatio);
            if ismeanRatioInf
                meanRatio = max(allP(isfinite(allP)));
            end
            
            all_fp_changes = S.sorting_falsePosRate_change;
%             fp_pct_before = [S.sorting_falsePosRate_before];
%             fp_pct_after = [S.sorting_falsePosRate_after];
%             %         sorting_falsePosRate_change = -[allStats.sorting_falsePosRate_change];
%             all_fp_changes = (fp_pct_before - fp_pct_after) ./ fp_pct_before;
            median_falsePosRate_change = mean(all_fp_changes);
            
            meanRatios(sort_fet_idx, prune_fet_idx) = meanRatio;
            
            
            median_falsePosRate_changes(sort_fet_idx, prune_fet_idx) = median_falsePosRate_change;
            ismeanRatioInf(sort_fet_idx, prune_fet_idx) = ismeanRatioInf;

            
        end
    end
%%
    fmtSetNames = @(s) upper(strrep(strrep(s, 'PCA', ''), 'u', ''));
    figure(opt.figId); clf;    
%     subplot(1,100, [10:95])
    
    imagesc(meanRatios);
    sortingFeatureSetLabels = cellfun(fmtSetNames, sortingFeatureSets, 'un', 0);
    pruningFeatureSetLabels = cellfun(fmtSetNames, pruningFeatureSets, 'un', 0);
    
    set(gca, 'ytick', 1:nSortingSets, 'xtick', 1:nPruningSets, 'yticklabel', sortingFeatureSetLabels, 'xticklabel', pruningFeatureSetLabels, 'fontsize', 10);   
    xticklabel_rotate;
%     ylabel('Sorting Features'); xlabel('Pruning Features');
    title('Mean Precision (P)');    
    colorbar;
    axis square;
    3;
    
    
    figure(opt.figId+10); clf;
    imagesc(median_falsePosRate_changes*100);
    set(gca, 'ytick', 1:nSortingSets, 'xtick', 1:nPruningSets, 'yticklabel', sortingFeatureSetLabels, 'xticklabel', pruningFeatureSetLabels, 'fontsize', 10);
%     ylabel('Sorting Features'); xlabel('Pruning Features');
    xticklabel_rotate;
    title('Mean Effectiveness (E), (%) ');
    axis square
    colorbar;
                
    3;
end



% function n = getNfromFet(fet_str)
%     if any(strncmpi(fet_str, {'PCA', 'GLF'}, 3))
%         n = str2double(fet_str(7:end));
%         cat_char = fet_str(6);
%         if strcmp(cat_char, 's');
%             n = n*4;
%         end        
%         
%     elseif strcmpi(fet_str, 'Neg')
%         n = 4;
%     end
% 
% end

function [proj_type, ccw_tf, concat_tf, nFet] = getFeatureDetails(fet_str)
    % PCAcur / cuw
    if any(strncmpi(fet_str, {'PCA', 'GLF'}, 3))
        proj_type = double( strcmp(fet_str(1), 'P') );
                
        nFet = str2double(fet_str(7:end));
        cat_char = fet_str(4);
        norm_char = fet_str(5); %#ok<NASGU>
        ccw_char = fet_str(6);
        
        if strcmp(cat_char, 's');
            nFet = nFet*4;
        end        
                
        
        ccw_tf = strcmp(ccw_char, 'w');
        concat_tf = strcmp(cat_char, 'c');        
        
    elseif strcmpi(fet_str, 'Neg')        
        
        proj_type = 2;
        ccw_tf = false;
        concat_tf = false;
        nFet = 4;
        
    end


end



function [allP, allE, allfFPremain, allfTPremain, uCellIds] = gatherP_E_vs_allVariables(...
    allStats_SPG_C_allSortFet, sortingFeatureSets_all, pruningFeatureSets_all, pruningMode, opt, pruning_filters)


    anyData_all = ~cellfun(@isempty, allStats_SPG_C_allSortFet);
    idx_gids_withData = find(any(any(anyData_all, 1),2));
    nGids = length(idx_gids_withData);

    allStats_SPG_C_allSortFet = allStats_SPG_C_allSortFet(:, :, idx_gids_withData);

    allStats_SPG_C = allStats_SPG_C_allSortFet;
    
    %%
%     sortingFeature_use = 'PCAcuw4';
%     sf_idx_use = strcmp(sortingFeatureSets, sortingFeature_use);
    sf_idx_use = 1:length(sortingFeatureSets_all); %  compute all now, but can select later
%     sf_idx_use = 6;  % [1, 3, 6, 12];
    allStats_SPG_C = allStats_SPG_C(sf_idx_use, :,:);
    
    sortingFeatureSets = sortingFeatureSets_all;
%     anyData_all = anyData_all(:,:,idx_gids_withData);

%     allStats_used_C = cell(1,nGids_usable);
%     idx_use = zeros(1, nGids_usable);

%             nPrunedTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpikesPruned'});
%             nSpkTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpkInClusts_Tot'});

    %%
        
    
%     pruningFeatureSets = pruningFeatureSets_all(idx_pruning_fet_use);
    pruningFeatureSets = pruningFeatureSets_all; 
    [pruning_PCAproj, pruning_ccw, pruning_concat, pruning_nFet] = cellfun(@getFeatureDetails, pruningFeatureSets);
    allStats_SPG_C = allStats_SPG_C_allSortFet; %(:,idx_pruning_fet_use,:);
%     [proj_type, ccw_tf, concat_tf, nFet] = getFeatureDetails(fet_str);
    
    
    [uProj, projTypeIdxs] = uniqueList(pruning_PCAproj);
    [uConcat, concatIdxs] = uniqueList(pruning_concat);
    [uCCW, ccwTypeIdxs] = uniqueList(pruning_ccw);
    [uNPruneFet, pruneFetIdxs] = uniqueList(pruning_nFet);    
    
    
%     maxN = max(uNPruneFet);
%     nMaxPerPruneFet = max( cellfun(@length, pruneFetIdxs) );
    
    nProjTypes = length(uProj);    
    nConcatTypes = length(uConcat);
    nCCWtypes = length(uCCW);
    nNPruneFets = length(uNPruneFet);
    
    
    %%
    sortingFeatureSets_use = sortingFeatureSets(sf_idx_use); %(1:2);
    nSortUse = length(sortingFeatureSets_use);
    
%     sort_fet_idx = 2;

    allS = [allStats_SPG_C{:}];
    uCellIds = unique( [allS.GC] );
    uCellIds_Gid = round(uCellIds/1000);
    uGids = unique(uCellIds_Gid);
    nCells = length(uCellIds);

    %%
    allNan  = nan(  nCells, nNPruneFets, nProjTypes, nConcatTypes, nCCWtypes, nSortUse, 'single');
%     allBlnkCell = cell( nCells, nNPruneFets, nProjTypes, nCCWtypes, nSortUse);
    
%     [avRatios, avFPchanges, avTPchanges, avfFPremain, avfTPremain] = deal( allNan );        
    [allP, allE, allTPchanges, allfFPremain, allfTPremain]  = deal( allNan );

%     [avRatios, avFPchanges, avTPchanges, avfFPremain, avfTPremain] = deal( allNan );        
%     [allP, allE, allTPchanges, allfFPremain, allfTPremain]  = deal( allNan );
    
    
    %%
%     allE   = cell( nNPruneFets, nCells);
%     allYields      = cell( nNPruneFets, nCells);    
%     sf_i = 1;
    
%     [avRatios, avYields, avFPchanges] = deal( nan( nNPruneFets, nCells) );    
%     [allP, allE, allYields]      = deal( cell( nNPruneFets, nCells) );
%%
%     intersect_3 = @(a,b,c) intersect(intersect(a,b),c);
%     intersect_4 = @(a,b,c,d) intersect(intersect(intersect(a,b),c),d);

    for sf_i = 1:nSortUse
            
       for g_i = 1:nGids
        
           for ccw_i = 1:nCCWtypes
%                ccw_idxs = ccwTypeIdxs{ccw_i};
               
               ccw_tf = (pruning_ccw == uCCW(ccw_i));
               

                % prune_ccw_idxs = intersect(prune_n_idxs, ccw_idxs);
                for proj_i = 1:nProjTypes
%                     proj_idxs = projTypeIdxs{proj_i};
                    
                    proj_tf = (pruning_PCAproj == uProj(proj_i));
                    

                    for cat_i = 1:nConcatTypes
%                         concat_idxs = concatIdxs{cat_i};
                        
                        concat_tf = (pruning_concat == uConcat(cat_i));

                        for ni = 1:nNPruneFets
%                             n_idxs = pruneFetIdxs{ni};
                            
                            n_tf = (pruning_nFet == uNPruneFet(ni));

%                             pruneFet_idxs = intersect_4(n_idxs, ccw_idxs, proj_idxs, concat_idxs);
                            pruneFet_idxs = find(n_tf & ccw_tf & proj_tf & concat_tf);                            
%                             assert(isequal(pruneFet_idxs, pruneFet_idxs2));
                            
                            assert(length(pruneFet_idxs) <= 1);
                            if isempty(pruneFet_idxs)
                                continue;
                            end


                            allStats_C = [allStats_SPG_C{sf_i, pruneFet_idxs, g_i}];                                                                       
                            if isempty(allStats_C)
                                continue;
                            end                        
                            resultStruct = getResultsStruct(allStats_C, opt, pruningMode);
                            idx_cell = binarySearch(uCellIds, resultStruct.GC);
                            assert(all(diff(idx_cell) > 0))
%                             [uIdx_cell, idx_cell_idxs] = uniqueList(idx_cell);


                            [cells_p, cells_e, cells_frTP,  cells_fFPrem, cells_fTPrem] = getRatios_FP_TP_changes(resultStruct);      

                            for ci = 1:length(idx_cell)
    %                         allNan  = nan(  nCells, nNPruneFets, nProjTypes, nCCWtypes, nSortUse);

                                allP(idx_cell(ci), ni, proj_i, cat_i, ccw_i, sf_i) = cells_p(ci);
                                allE(idx_cell(ci), ni, proj_i, cat_i, ccw_i, sf_i) = cells_e(ci); 

    %                             allTPchanges{ni, uIdx_cell(ci), sf_i, proj_i} = all_frTP(idx_cell_idxs{ci}); 
                                allfFPremain(idx_cell(ci), ni, proj_i, cat_i, ccw_i, sf_i) = cells_fFPrem(ci); 
                                allfTPremain(idx_cell(ci), ni, proj_i, cat_i, ccw_i, sf_i) = cells_fTPrem(ci);                     

    %                             avRatios    (ni, uIdx_cell(ci), sf_i, proj_i) = nanmedian(allP(idx_cell_idxs{ci}));
    %                             avFPchanges (ni, uIdx_cell(ci), sf_i, proj_i) = nanmean(all_frFP(idx_cell_idxs{ci})); %medianFPchange;
    %                             avTPchanges (ni, uIdx_cell(ci), sf_i, proj_i) = nanmean(all_frTP(idx_cell_idxs{ci})); %medianFPchange;
    %                             avfFPremain (ni, uIdx_cell(ci), sf_i, proj_i) = nanmean(allfFPrem(idx_cell_idxs{ci})); %medianFPchange;
    %                             avfTPremain (ni, uIdx_cell(ci), sf_i, proj_i) = nanmean(allfTPrem(idx_cell_idxs{ci})); %medianFPchange;

                            end

                        end    
                    end
                end
            end
        end
    end
    %%
%     if nSortUse > 1
%         avRatios_vN = reshape(avRatios, [nNPruneFets, maxN*nSortUse, ]);
%         avFPchanges_vN = reshape(avFPchanges, [nNPruneFets, maxN*nSortUse]);
%         avYields_vN = reshape(avYields, [nNPruneFets, maxN*nSortUse]);
%                 
%     end
    %%
%     
%     sf_idxs_use = 3:5;
%     
%     if nSortUse > 1
%         allP_use_vN = reshape(avRatios, [nNPruneFets, maxN*nSortUse, ]);
%         avFPchanges_vN = reshape(avFPchanges, [nNPruneFets, maxN*nSortUse]);        
%                 
%     end
    
    %%    
%     [sorting_PCAproj, sorting_ccw, sorting_cat, sorting_nFet] = cellfun(@getFeatureDetails, sortingFeatureSets_use);

    
end

function showPerformance_Vs_variables(allStats_SPG_C_all, sortingFeatureSets_all, pruningFeatureSets_all, pruningMode, opt)
    

%     pruning_filters = struct('ccw', 1, 'cat', 1);
    
    [pruning_PCAproj_all, pruning_ccw_all, pruning_cat_all, pruning_nFet_all] = cellfun(@getFeatureDetails, pruningFeatureSets_all);
    
%     idx_pruning_fet_use = pruning_ccw_all; % & pruning_cat_all;
    idx_pruning_fet_use = pruning_cat_all;
%     idx_pruning_fet_use = true(size(pruning_cat_all));
    

    allStats_SPG_C = allStats_SPG_C_all(:,idx_pruning_fet_use,:);
    
    %     if isfield(pruning_filters, 'ccw') && ~isempty(pruning_filters.ccw);
%         idx_pruning_fet_use = idx_pruning_fet_use & (pruning_ccw == pruning_filters.ccw);
%     end
%     if isfield(pruning_filters, 'cat') && ~isempty(pruning_filters.cat);
%         idx_pruning_fet_use = idx_pruning_fet_use & (pruning_concat == pruning_filters.cat);
%     end

    pruningFeatureSets = pruningFeatureSets_all(idx_pruning_fet_use);
    [pruning_PCAproj, pruning_ccw, pruning_concat, pruning_nFet] = cellfun(@getFeatureDetails, pruningFeatureSets);

    uProj= unique(pruning_PCAproj);     %nProjTypes = length(uProj);      
    uCCW = uniqueList(pruning_ccw);     %nCCWtypes = length(uCCW);
    uNPruneFet = unique(pruning_nFet);  %nNPruneFets = length(uNPruneFet);
    uConcat = unique(pruning_concat);
    %%    
%     profile on;
    [allP, allE, allfFPremain, allfTPremain, uCellIds] = gatherP_E_vs_allVariables(...
        allStats_SPG_C, sortingFeatureSets_all, pruningFeatureSets, pruningMode, opt);
%     profile viewer;
    %%    
    [nCells, nNPruneFets, nProjTypes, nConcatTypes, nCCWtypes, nSortingFet] = size(allP);        
    assert(isequal([nProjTypes, nConcatTypes, nCCWtypes, nNPruneFets], [length(uProj), length(uConcat), length(uCCW), length(uNPruneFet) ]));
    
    
    

    
    %%    
    [sorting_PCAproj, sorting_ccw, sorting_cat, sorting_nFet] = cellfun(@getFeatureDetails, sortingFeatureSets_all);
    
    %%
%     sortFet_idxs = find(strcmp(sortingFeatureSets_all, sf_use), 1);
%     iii = iii+1;
    sortFet_idxs = 4; %1:length(sortingFeatureSets_use);
%     sortFet_idxs = find(ibetween(sorting_nFet, ));
%     sf_use = 'PCAcuw4';
    
%     sortFet_idxs = [5:20];
    sortFet_idxs = 1:length(sortingFeatureSets_all);
    
    sortFet_idxs = 6; %[1:10];

%     sortFet_idxs = sortFet_idxs + 10;
%     sortFet_idxs = sortFet_idxs(sortFet_idxs <= length(sortingFeatureSets_use));
        
    nSortFetPlot = length(sortFet_idxs);
    addToPlots = 1;
    makePlots = 0;
    
    
    addToPlots = 1;
    
    plotInGrayScale = 0;
    
    doPlotN = 1;
    doPlotBinaryVariables = 1;
    
    
    
    if doPlotN

        
        separateGLFandPCA = 1;
        plotAvAllOnTop = 0;
        
        if plotAvAllOnTop
%             sortFet_idxs = 1:length(sortingFeatureSets_all);
        end
        nMinHaveN = 5;
%         idx_fromHave_useN = 4;

        measures = {allP, allE};
%         measure_label_short = {'P', 'E'};
        measure_ylabels = {'Precision (P)', 'Effectiveness (E)'};
        
        av_op = {@nanmedian, @nanmean};  av_std_op = {@getMedianPrctle, @getMeanStd};        
        
%         av_op = {@nanmedian, @nanmedian}; av_std_op = {@getMedianPrctle, @getMedianPrctle};        

        nMeasures = length(measures);

        figure(903); clf;

        for mi = 1:nMeasures

            allM = measures{mi};
            mnStdFunc = av_std_op{mi};
            meanFunc = av_op{mi};
            justUseConcat = 1;
            
            
            % take average over ccw,(& pca)
            % average over cells
            % separate out for each n

            allM_c_n_oth = reshape(allM(:,:,:, :,:, sortFet_idxs), nCells, nNPruneFets, [nProjTypes * nCCWtypes * nConcatTypes * nSortFetPlot], 1);
            idx_cell_use_i = find( sum(sum( ~isnan(allM_c_n_oth),2),3) >= nProjTypes*nMinHaveN )';
            
            if ~separateGLFandPCA

                allM_c_n_oth = allM(idx_cell_use_i,:,:,:,:, sortFet_idxs);                

                [ms_vs_N_med, ms_vs_N_mad] = mnStdFunc(apply2dims(meanFunc, allM_c_n_oth, [3 4 5 6]), 1);

%                 allM_c_n_oth = reshape(allM(:,:,:, :,:, sortFet_idxs), nCells, nNPruneFets, [nProjTypes * nCCWtypes * nConcatTypes * nSortFetPlot], 1);
%                 idx_cell_use_i = find( sum(sum( ~isnan(allM_c_n_oth),2),3) >= nProjTypes*nMinHaveN )';
            %         idx_cell_use_P = idx_cell_use_P(idx_fromHave_useN);


%                 allM_c_n_oth = allM_c_n_oth(idx_cell_use_i,:,:);
%                 [Ms_vs_N_av, Ms_vs_N_std] = mnStdFunc(meanFunc(allM_c_n_oth, 3), 1);                
                
                subplot(nMeasures,1,mi); hold on; box on;
                errorbar(uNPruneFet, ms_vs_N_med, ms_vs_N_mad, ['rs-']);

            elseif separateGLFandPCA
                allM_c_n_oth_PCA = allM(idx_cell_use_i,:,uProj==1,:,:, sortFet_idxs);
                allM_c_n_oth_GLF = allM(idx_cell_use_i,:,uProj==0,:,:, sortFet_idxs);

                [ms_vs_N_med_PCA, ms_vs_N_mad_PCA] = mnStdFunc(apply2dims(meanFunc, allM_c_n_oth_PCA, [4 5 6]), 1);
                [ms_vs_N_med_GLF, ms_vs_N_mad_GLF] = mnStdFunc(apply2dims(meanFunc, allM_c_n_oth_GLF, [4 5 6]), 1);

                subplot(nMeasures,1,mi); hold on; box on;
                errorbar(uNPruneFet, ms_vs_N_med_PCA, ms_vs_N_mad_PCA, ['bs-']);
                errorbar(uNPruneFet, ms_vs_N_med_GLF, ms_vs_N_mad_GLF, ['rs-']);
            end
                
            %         P_use = allP(:,:,:,:,sortFet_idxs(sf_i));

            %         plot(uNPruneFet,
            
            %         errorbar(uNPruneFet, nanmean(allP(:,:,sortFet_idxs(sf_i), 1), dim), nanstderr(avRatios(:,:,sortFet_idxs(sf_i), 1), dim), [color_s(sf_i) 's-']);

            %         errorbar(uNPruneFet, nanmean(avRatios(:,:,sortFet_idxs(sf_i), 2), dim), nanstderr(avRatios(:,:,sortFet_idxs(sf_i), 2), dim), [color_s(sf_i+2) 'o:']);

            xticks = [1, 4:4:16];
            xlims = [0, 17];
            L1 = ms_vs_N_med_GLF-ms_vs_N_mad_GLF; U1 = ms_vs_N_med_GLF+ms_vs_N_mad_GLF;  L2 = ms_vs_N_med_PCA-ms_vs_N_mad_PCA; U2 = ms_vs_N_med_PCA+ms_vs_N_mad_PCA; 
            
            ylims_P = lims([0; L1(:); L2(:); U1(:); U2(:)], .05);
            set(gca, 'xtick', xticks, 'xlim', xlims, 'ylim', ylims_P);
            xlabel('# Pruning Features'); ylabel(measure_ylabels{mi});

            if mi == 2
               legend({'PCA', 'GLF'}, 'location', 'NE');
            end

        end

        %     axis tight



    %     for sf_i = 1:1;%nSortFetPlot
    % 
    %         allE_c_n_oth = reshape(allE(:,:,:,:, sortFet_idxs), nCells, nNPruneFets, [nProjTypes * nCCWtypes * nSortFetPlot], 1);
    %         idx_cell_use_E = find( sum(sum( ~isnan(allE_c_n_oth),2),3) >= nProjTypes*nMinHaveN )';
    %         %         idx_cell_use_E = idx_cell_use_E(idx_fromHave_useN);
    %         [P_vs_N_med, P_vs_N_mad] = getMedianMAD(nanmedian(allE_c_n_oth, 3), 1);
    % 
    % 
    %         allE_c_n_oth_PCA = allE(idx_cell_use_E,:,uProj==1,:, sortFet_idxs);
    %         allE_c_n_oth_GLF = allE(idx_cell_use_E,:,uProj==0,:, sortFet_idxs);
    % 
    %         [E_vs_N_med_PCA, E_vs_N_mad_PCA] = getMeanStd(apply2dims(@nanmedian, allE_c_n_oth_PCA, [4 5]), 1);
    %         [E_vs_N_med_GLF, E_vs_N_mad_GLF] = getMeanStd(apply2dims(@nanmedian, allE_c_n_oth_GLF, [4 5]), 1);
    % 
    % 
    %         errorbar(uNPruneFet, E_vs_N_med_PCA, E_vs_N_mad_PCA, ['bs-']);
    %         errorbar(uNPruneFet, E_vs_N_med_GLF, E_vs_N_mad_GLF, ['rs-']);
    % 
    %         %         errorbar(uNPruneFet, nanmean(avFPchanges(:,:,sortFet_idxs(sf_i), 1), dim), nanstderr(avFPchanges(:,:,sortFet_idxs(sf_i), 1), dim), [color_s(sf_i) 's-']);
    %         %         errorbar(uNPruneFet, nanmean(avFPchanges(:,:,sortFet_idxs(sf_i), 2), dim), nanstderr(avFPchanges(:,:,sortFet_idxs(sf_i), 2), dim), [color_s(sf_i+2) 'o:']);
    %     end


        figure(904); clf; hold on;
        %     sortFet_idxs = 1:nSortFetPlot; %7; %nSortUse;
        %     nSortFetPlot = length(sortFet_idxs);
        %     sortFet_idxs = 3;

        for sf_i = 1:1; %1:nSortFetPlot; %nSortFetPlot
            %average
            allFP_PCA = allfFPremain(:,:,uProj==1,uConcat==1,:, sortFet_idxs);
            allTP_PCA = allfTPremain(:,:,uProj==1,uConcat==1,:, sortFet_idxs);

            xx_PCA = apply2dims(@nanmean, allFP_PCA, [1 3 4 5 6]);
            yy_PCA = apply2dims(@nanmean, allTP_PCA, [1 3 4 5 6]);


            allFP_GLF = allfFPremain(:,:,uProj==0,uConcat==1,:, sortFet_idxs);
            allTP_GLF = allfTPremain(:,:,uProj==0,uConcat==1,:, sortFet_idxs);

            xx_GLF = apply2dims(@nanmean, allFP_GLF, [1 3 4 5 6]);
            yy_GLF = apply2dims(@nanmean, allTP_GLF, [1 3 4 5 6]);


            %         xx1 = nanmean(avfFPremain(:,:,sortFet_idxs(sf_i), 1), dim); xe1 = nanstderr(avfFPremain(:,:,sortFet_idxs(sf_i), 1), dim);
            %         yy1 = nanmean(avfTPremain(:,:,sortFet_idxs(sf_i), 1), dim); ye1 = nanstderr(avfTPremain(:,:,sortFet_idxs(sf_i), 1), dim);

            %         xx2 = nanmean(avfFPremain(:,:,sortFet_idxs(sf_i), 2 ), dim); xe2 = nanstderr(avfFPremain(:,:,sortFet_idxs(sf_i), 2), dim);
            %         yy2 = nanmean(avfTPremain(:,:,sortFet_idxs(sf_i), 2 ), dim); ye2 = nanstderr(avfTPremain(:,:,sortFet_idxs(sf_i), 2), dim);

            %         ploterr(xx1, yy1, xe1, ye1, [color_s(sf_i) ':']);
            %         ploterr(xx1, yy1, xe1, ye1, [color_s(sf_i) ':']);

            hnd(1) = plot(xx_PCA, yy_PCA, [color_s(sf_i) 's-']);
            plot(xx_PCA(1), yy_PCA(1), [color_s(sf_i) 's-'], 'markerfacecolor', color_s(sf_i));

            hnd(2) = plot(xx_GLF, yy_GLF, [color_s(sf_i+2) 'o:']);
            plot(xx_GLF(1), yy_GLF(1), [color_s(sf_i+2) 'o:'], 'markerfacecolor', color_s(sf_i+2));


            %         errorbar(xx, yy, ye, ['g']);
            %         plot(, nanmean(avfTPremain(:,:,sortFet_idxs(sf_i)), dim), );
            %         plot(nanmean(avFPchanges(:,:,sortFet_idxs(sf_i)), dim), nanmean(avTPchanges(:,:,sortFet_idxs(sf_i)), dim), [color_s(sf_i) 'o:']);
        end
        %     set(gca, 'xtick', [4:4:16]);
        xlabel('Fraction of FP remaining'); ylabel('Fraction of TP remaining');
        L = lims([xx_PCA(:); yy_PCA(:); xx_GLF(:); yy_GLF(:)], .05);
        axis([L; L]);
        %     axis([0 1 0 1]);
        fplot(@(x) x, L, 'k:');
        axis square; box on;
        legend(hnd, {'PCA', 'GLF'}, 'location', 'SE');
        3;

    end
%     testAllCombos(a,b,c);
    
    if doPlotBinaryVariables
        makePlots = 1;
        addStatsStrToXlabel = 0;
        plotErrorBars = 0;

    %     binaryVariablesToPlot = {'PCA_vs_GLF', 'ccw_vs_raw', 'cat_vs_sep'};
        binaryVariablesToPlot = {'PCA_vs_GLF', 'ccw_vs_raw'};
    %     binaryVariablesToPlot = {'PCA_vs_GLF'};

        [pval_t, pval_w] = deal( nan(length(sortingFeatureSets_all), length(binaryVariablesToPlot), 2));

        if plotInGrayScale
            mk_col = 'k';
            unity_line_col = 'k';
        else
            mk_col = 'b';
            unity_line_col = 'r';
        end

        for sf_i = 1; %1:length(sortingFeatureSets_all)

    %         sortFet_idxs = sf_i;
            for vi = 1:length(binaryVariablesToPlot)
                bvariable = binaryVariablesToPlot{vi};


                switch bvariable
                    case 'PCA_vs_GLF', label1 = 'PCA'; label2 = 'GLF';
                        fig_id = 100;
                        uVals = uProj;
                    case 'ccw_vs_raw', label1 = 'Raw'; label2 = 'CCW';
                        fig_id = 200;   
                        uVals = uCCW;
                    case 'cat_vs_sep', label1 = 'concatenated'; label2 = 'separate';
                        fig_id = 300;
                        uVals = uConcat;                
                end



                nMinHaveN = 5;
        %         idx_fromHave_useN = 4;

                measures = {allP, allE};
        %         measure_label_short = {'P', 'E'};
                measure_ylabels = {'Precision', 'Effectiveness'};
                measure_abbrev = {'P', 'E'};

                av_op_names = {'Median', 'Mean'};

    %             av_op = {@nanmedian, @nanmedian}; av_std_op = {@getMedianMAD, @getMedianMAD};        

                nMeasures = length(measures);
    %             nMeasures = 1;

                if makePlots
                    figure(fig_id); clf;
                end

                if length(uVals) == 1
                    fprintf('Can''t plot %s if only have uVal = %d ...\n', bvariable, uVals);
                    continue;
                end        

        %         splitGLF_PCA = 0;

                for mi = 1:nMeasures

                    allM = measures{mi};
                    switch av_op_names{mi}
                        case 'Median',                                                                 
    %                         meanFunc = @nanmedian;
                            mnStdFunc = @getMedianPrctle;
                            pm_str = '\pm p25';
                        case 'Mean'
    %                         meanFunc = @nanmean;
                            mnStdFunc = @getMeanStd;
                            pm_str = '\pm stdev';
                    end

                    % take average over ccw,(& pca)
                    % average over cells
                    % separate out for each n


        %             allM_c_n_oth = reshape(allM(:,:,:,:, sortFet_idxs), nCells, nNPruneFets, [nProjTypes * nCCWtypes * nSortFetPlot], 1);
        %             idx_cell_use_i = find( sum(sum( ~isnan(allM_c_n_oth),2),3) >= nProjTypes*nMinHaveN )';
                    %         idx_cell_use_P = idx_cell_use_P(idx_fromHave_useN);
        %             

        %             allM_c_n_oth = allM_c_n_oth(idx_cell_use_i,:,:);
        %             [Ms_vs_N_av, Ms_vs_N_std] = mnStdFunc(meanFunc(allM_c_n_oth, 3), 1);

                    switch bvariable
                        case 'PCA_vs_GLF',
                            allM1 = allM(:,:,uProj==1,:,:, sortFet_idxs);
                            allM2 = allM(:,:,uProj==0,:,:, sortFet_idxs);
                        case 'ccw_vs_raw',
                            allM1 = allM(:,:,:,:,uCCW==0, sortFet_idxs);
                            allM2 = allM(:,:,:,:,uCCW==1, sortFet_idxs);                    
                        case 'cat_vs_sep',
                            allM1 = allM(:,:,:,uConcat==1,:, sortFet_idxs);
                            allM2 = allM(:,:,:,uConcat==0,:, sortFet_idxs);                    
                    end

                    idx_both_ok = ~isnan(allM1) & ~isnan(allM2);
        %             allM1 = allM1(idx_both_ok);
        %             allM2 = allM2(idx_both_ok);
                    allM1(~idx_both_ok) = nan;
                    allM2(~idx_both_ok) = nan;

                    nDim2 = numel(allM1)/nCells;
                    allM1_v_cell = reshape(allM1, [nCells, nDim2]);
                    allM2_v_cell = reshape(allM2, [nCells, nDim2]);

                    [ms1_med, ms1_L, ms1_U] = mnStdFunc(allM1_v_cell, 2, 1);
                    [ms2_med, ms2_L, ms2_U] = mnStdFunc(allM2_v_cell, 2, 1);
                    idx_nonnan = ~isnan(ms1_med) & ~isnan(ms2_med);
                    ms1_med = ms1_med(idx_nonnan); ms1_L = ms1_L(idx_nonnan); ms1_U = ms1_U(idx_nonnan);
                    ms2_med = ms2_med(idx_nonnan); ms2_L = ms2_L(idx_nonnan); ms2_U = ms2_U(idx_nonnan);

                    [stats_str, p_w(sf_i, vi, mi), p_t(sf_i, vi, mi)] = getStatsStrs(ms1_med, ms2_med, label1, label2);

                    if makePlots
                        subplotGap(1, nMeasures, 1, mi); hold on; box on;
                        if ~plotErrorBars
                            plot(ms1_med, ms2_med, 'o', 'color', mk_col);            
                        else
                            ploterr(ms1_med, ms2_med, {ms1_L, ms1_U}, {ms2_L, ms2_U}, 'bo', 'hh', 1.5);                                    
                        end



            %             xticks = [1, 4:4:16];
            %             xlims = [0, 17];
            %             ylims_P = [0, max([ms_vs_N_med_GLF+ms_vs_N_mad_GLF, ms_vs_N_med_PCA+ms_vs_N_mad_PCA])*1.05];
            %             set(gca, 'xtick', xticks, 'xlim', xlims, 'ylim', ylims_P);
    %                     xlabel(label1); 
    %                     ylabel(label2);
                        title(sprintf('%s (%s) for each cluster (%s %s)', measure_ylabels{mi}, measure_abbrev{mi}, av_op_names{mi}, pm_str));

                        if ~plotErrorBars
                            L = lims([ms1_med; ms2_med], .05 );
                        else
                            L = lims([ms1_L; ms1_U; ms2_L; ms2_U], .05);
                        end
                        axis([L; L]);
                        % errorbar(x,y,distance from y)  % not (x,y,abs lower y, abs upper y)


                        hold on;
                        [xx, yy] = fplot(@(x) x, xlim);
                        plot(xx, yy, 'color', unity_line_col, 'linestyle', ':');
            %             title(sprintf('Mean %s %s for each cluster', measure_labels{mi}, measure_label_short{mi}));        
                        axis square;    

                        xlabel_str = sprintf('%s (%s)', measure_abbrev{mi}, label1);
                        if addStatsStrToXlabel
                            xlabel_str = [xlabel_str, stats_str]; %#ok<AGROW>
                        end
                        xlabel(xlabel_str);

                        ylabel_str = sprintf('%s (%s)', measure_abbrev{mi}, label2);
                        ylabel(ylabel_str);
            %     xlabel(label1);     
                    end



        % 
        %              figure(50+10*(mi-1)   + 1); clf;
        %             barPlotTwoSamples(m1_forClusts, m2_forClusts, label1, label2);
        %             xlabel(sortingFeatureSets_use{1});            

                end        

            end
        end

        if 0
            %%
            figure(610); clf;
            plot(-log10(p_t(:,:,1)), 'o-'); hold on;
            plot(-log10(p_t(:,:,2)), 's:'); 

            xlabel('sorting feature idx'); ylabel('-log_{10} p_T'); title('-log pval for mean'); 
            leg_str = [cellfun(@(s) [s ' (P)'], binaryVariablesToPlot, 'un', 0), cellfun(@(s) [s ' (E)'], binaryVariablesToPlot, 'un', 0)]
            legend(leg_str, 'interpreter', 'none', 'location', 'bestoutside');
            drawHorizontalLine([-log10(.05) 2 3], 'linestyle', ':', 'color', 'k');

            figure(611);
            clf;
            plot(-log10(p_w(:,:,1)), 'o-'); hold on;
            plot(-log10(p_w(:,:,2)), 's:'); 

            xlabel('sorting feature idx'); ylabel('-log_{10} p_W'); title('-log pval for median'); 
            legend(leg_str, 'interpreter', 'none', 'location', 'bestoutside')


        end
    3;
    end
3;
    
end


% function showPerformance_Vs_NFeatures(allP, allE, uNPruneFet, sortFet_idxs)
% 
% end


function showPerformance_Vs_X(allStats_SPG_C, sortingFeatureSets, pruningFeatureSets, pruningMode, opt, X)

    anyData_all = ~cellfun(@isempty, allStats_SPG_C);
    idx_gids_withData = find(any(any(anyData_all, 1),2));
    nGids_usable = length(idx_gids_withData);

    allStats_SPG_C = allStats_SPG_C(:, :, idx_gids_withData);
    
%     anyData_all = anyData_all(:,:,idx_gids_withData);

%     allStats_used_C = cell(1,nGids_usable);
%     idx_use = zeros(1, nGids_usable);

%             nPrunedTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpikesPruned'});
%             nSpkTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpkInClusts_Tot'});

        
    [pruning_projTypes, pruning_ccw_tf, pruning_cat_tf, pruning_nFet] = cellfun(@getFeatureDetails, pruningFeatureSets);        
    
        
%     idx_4 = find( cellfun(@(s) ~isempty(regexp(s, 'PCAc..4')) || ~isempty(regexp(s, 'PCAs..1')) || strcmpi(s, 'Neg'), pruningFeatureSets) ); %#ok<RGXP1>
%     idx_8 = find( cellfun(@(s) ~isempty(regexp(s, 'PCAc..8')) || ~isempty(regexp(s, 'PCAs..2')), pruningFeatureSets) ); %#ok<RGXP1>
%     idx_12 = find( cellfun(@(s) ~isempty(regexp(s, 'PCAc..12')) || ~isempty(regexp(s, 'PCAs..3')), pruningFeatureSets) ); %#ok<RGXP1>
%     idx_16 = find( cellfun(@(s) ~isempty(regexp(s, 'PCAc..16')) || ~isempty(regexp(s, 'PCAs..4')), pruningFeatureSets) ); %#ok<RGXP1>
%     all_Ns = [4 8 12 16];
%     pruneFetIdxs = {idx_4, idx_8, idx_12, idx_16};
%   maxN = max(cellfun(@length, pruneFetIdxs));        
  %%
    [all_Ns, pruneFetIdxs] = uniqueList(pruning_nFet);
        
    all_pruning_features = [pruning_projTypes(:), pruning_ccw_tf(:), pruning_cat_tf(:), pruning_nFet(:)];
%     X = 'PCA_vs_GLF';
%     X = 'ccw_vs_raw';
    X = 'cat_vs_sep';
    
    switch X
        case 'PCA_vs_GLF', idx_1 = find(pruning_projTypes); label1 = 'PCA'; label2 = 'GLF';
            col_idx = 1;             
        case 'ccw_vs_raw', idx_1 = find(pruning_ccw_tf);  label1 = 'CCW'; label2 = 'Raw';
            col_idx = 2;            
        case 'cat_vs_sep', idx_1 = find(pruning_cat_tf);  label1 = 'concatenated'; label2 = 'separate';
            col_idx = 3;            
            
    end
    idx_2 = zeros(1, length(idx_1));
    other_cols = setdiff(1:4, col_idx);    
    for i = 1:length(idx_1)        
        features_2_orig = all_pruning_features(idx_1(i), :);
        features_2_opp = features_2_orig;
        features_2_opp(col_idx) = ~features_2_opp(col_idx); 
        
        idx_2_i = findRows(features_2_opp, all_pruning_features);
        
        if ~isempty(idx_2_i)            
            idx_2(i) = idx_2_i;
        else
            3;
        end
    end
    
    idx_remove = idx_2 == 0;
    idx_1(idx_remove) = [];
    idx_2(idx_remove) = [];
    
    
                   
    %%
%     maxN = max(all_Ns);
    
%     nNPruneFets = length(all_Ns);
    %%
    sf_idxs = 7;
    sortingFeatureSets_use = sortingFeatureSets(sf_idxs); %(1:2);
    nSortUse = length(sortingFeatureSets_use);
    
%     sort_fet_idx = 2;
    
    

    
    allStats_SPG_C_1 = allStats_SPG_C(sf_idxs, idx_1, :);
    allStats_SPG_C_2 = allStats_SPG_C(sf_idxs, idx_2, :);
    opt.maximizeWhenChoosingGid = 'all';
    
    S1 = getResultsStruct(allStats_SPG_C_1, opt, pruningMode);
    S2 = getResultsStruct(allStats_SPG_C_2, opt, pruningMode);
    
    % 1 - plot all vs all
%     plot_all_v_all = 0;
%     if plot_all_v_all
%     end

    
%     allStats_SPG_C_1_v = allStats_SPG_C_1(:);
%     allStats_SPG_C_2_v = allStats_SPG_C_2(:);
%     idx_use = ~cellfun(@isempty, allStats_SPG_C_1_v) & ~cellfun(@isempty, allStats_SPG_C_2_v);
%     allStats_SPG_C_1_v = [allStats_SPG_C_1_v{idx_use}];
%     allStats_SPG_C_2_v = [allStats_SPG_C_2_v{idx_use}];
%     assert(isequal([allStats_SPG_C_1_v.Gid], [allStats_SPG_C_2_v.Gid]))
%     assert(isequal({allStats_SPG_C_1_v.clustIds}, {allStats_SPG_C_1_v.clustIds}))
%     N = length(allStats_SPG_C_1_v);

    allStats_SPG_C_1_v = S1(:);
    allStats_SPG_C_2_v = S1(:);
    N = length(allStats_SPG_C_1_v);
    
            
    
    
    p1_C = cell(1, N); p2_C = cell(1, N);
    e1_C = cell(1, N); e2_C = cell(1, N);
%     Gids_C = cell(1, N);
%     prune_Proj = cell(1,N);
    for i = 1:N
        [~, ia, ib] = intersect(allStats_SPG_C_1_v(i).GC, allStats_SPG_C_2_v(i).GC);
%         p_C{i} = allStats_SPG_C_1_v(i)./
        p1_C{i} = allStats_SPG_C_1_v(i).precision(ia);
        p2_C{i} = allStats_SPG_C_2_v(i).precision(ib);
        
        e1_C{i} = allStats_SPG_C_1_v(i).sorting_falsePosRate_change(ia);
        e2_C{i} = allStats_SPG_C_2_v(i).sorting_falsePosRate_change(ib);        
    end
    N_each = cellfun(@length, e2_C);
    idx_each = NtoIdx(N_each);

    e1 = [e1_C{:}]; e2 = [e2_C{:}];
    idx_ok = ~isnan(e1) & ~isnan(e2);
    e1 = e1(idx_ok);
    e2 = e2(idx_ok);
    
    p1 = [p1_C{:}]; p2 = [p2_C{:}];
    idx_ok = ~isnan(p1) & ~isnan(p2);
    p1 = p1(idx_ok);
    p2 = p2(idx_ok);
    
    

    statGCs = [allStats_SPG_C_1_v.GC];        
    GCs_each = statGCs(idx_each);
    GCs_each = GCs_each(idx_ok);    
    [uGCs, GC_idxs] = uniqueList(GCs_each);    
    
    
    
%     [uPruneFet, pruneFetIdxs] = uniqueList({allStats_SPG_C_1_v.pruningFeatures});
    
%     [proj_type_1, ccw_tf_1, concat_tf_1, nFet_1] = cellfun(@getFeatureDetails, {allStats_SPG_C_1_v.pruningFeatures});
%     [proj_type_2, ccw_tf_2, concat_tf_2, nFet_2] = cellfun(@getFeatureDetails, {allStats_SPG_C_2_v.pruningFeatures});
    
    
    % 1. simple scatter plot, all vs all.
%     figure(45); plot(e1, e2, 'o', 'markersize', 3);
%     title('Individual clusters (no averaging)');
%     
%     figure(46); clf; 
%     barPlotTwoSamples(e1, e2, label1, label2);
    
    % plot, but color code
    
    % do for [1] precison [2] effectiviness:
    XX = {{p1, p2}, {e1, e2}};
    measure_labels = {'Precision', 'Effectiveness'};
    measure_label_short = {'(P)', '(E)'};
    nMeasures = length(XX);
    for mi = 1:nMeasures
        m1 = XX{mi}{1};
        m2 = XX{mi}{2};
    
        figure(50+10*(mi-1))
    
        % plot average for each group
        m1_forClusts = cellfun(@(list) mean(m1(list)), GC_idxs);
        m2_forClusts = cellfun(@(list) mean(m2(list)), GC_idxs);

        plot(m1_forClusts, m2_forClusts, 'o');    
        L = lims([m1_forClusts, m2_forClusts], .05 );
        axis([L; L]);

        hold on;
        fplot(@(x) x, xlim, 'r:')
        title(sprintf('Mean %s %s for each cluster', measure_labels{mi}, measure_label_short{mi}));        
        axis square;    

        ylabel(label2);
    %     xlabel(label1);     

        stats_str = getStatsStrs(m1_forClusts, m2_forClusts, label1, label2);
        xlabel([label1, stats_str]);
        
         figure(50+10*(mi-1)   + 1); clf;
        barPlotTwoSamples(m1_forClusts, m2_forClusts, label1, label2);
        xlabel(sortingFeatureSets_use{1});
    end
    

    
    
    
   

    
    
    %%
    if col_idx == 1,   assert(all(proj_type_1 ~= proj_type_2));  else    assert(all(proj_type_1 == proj_type_2));   end
    if col_idx == 2,   assert(all(ccw_tf_1    ~= ccw_tf_2));     else    assert(all(ccw_tf_1    == ccw_tf_2));   end
    if col_idx == 3,   assert(all(concat_tf_1 ~= concat_tf_2));  else    assert(all(concat_tf_1 == concat_tf_2));   end
    assert(all(nFet_1 == nFet_2));
    
    
    %%
    
    
%     e1 = [allStats_SPG_C_1_v.sorting_falsePosRate_change];
%     e2 = [allStats_SPG_C_2_v.sorting_falsePosRate_change];
    
%     L1 = cellfun(@length, {allStats_SPG_C_1_v.sorting_falsePosRate_change});
%     L2 = cellfun(@length, {allStats_SPG_C_2_v.sorting_falsePosRate_change});
    %%
    
    % color
    
    
    
%     gscatter(e1, e2, 
    
    %%
    
    % GLF/PCA - symbol : square/ circle
    % ccw vs raw - 
    % cat vs indiv - filled vs not.
    % n - color;

    
    
    figure(56); 
    
    
    
    %%
% 
%     allmeanRatios_1 = nan( nNPruneFets, maxN, nSortUse);
%     allIsmeanRatiosInf = nan( nNPruneFets, maxN, nSortUse);
%     allMedianFPchanges = nan( nNPruneFets, maxN, nSortUse);
%         
%     for ni = 1:nNPruneFets
%         for pf_i = 1:length(pruneFetIdxs{ni})
%             
%             for sf_i = 1:nSortUse
%                 pf_idx = pruneFetIdxs{ni}(pf_i);
%                 pruneFet_str = pruningFeatureSets{pf_idx};
% 
%                 allStats_C = [allStats_SPG_C{sf_i, pf_idx, :}];    
%                 if isempty(allStats_C)
%                     continue;
%                 end            
%                 resultStruct = getResultsStruct(allStats_C, opt, pruningMode);            
%                 [meanRatio, medianFPchange, yield, ismeanRatiosInf] = getmeanRatio_FPChange(resultStruct);
%                 allmeanRatios(ni, pf_i, sf_i) = meanRatio;
%                 allIsmeanRatiosInf(ni, pf_i, sf_i) = ismeanRatiosInf;
%                 allMedianFPchanges(ni, pf_i, sf_i) = medianFPchange;
%             end
%         end
%     end
%     
%     if nSortUse > 1
%         allmeanRatios = reshape(allmeanRatios, [nNPruneFets, maxN*nSortUse, ]);
%         allMedianFPchanges = reshape(allMedianFPchanges, [nNPruneFets, maxN*nSortUse]);
%     end
%     
%     dim = 2;
%     figure(903); subplot(2,1,1);
%     errorbar(all_Ns, nanmean(allmeanRatios, dim), nanstd(allmeanRatios, [], dim), 'b.-');
%     set(gca, 'xtick', [4:4:16]);
%     xlabel('# Pruning Features'); ylabel('Median Precision (P)');
%     
%     subplot(2,1,2);
%     errorbar(all_Ns, nanmean(allMedianFPchanges, dim), nanstd(allMedianFPchanges, [], dim), 'b.-');
%     set(gca, 'xtick', [4:4:16]);
%     xlabel('# Pruning Features'); ylabel('Mean Effectiveness (E)');
    3;

end


function [stats_str, p_w, p_t] = getStatsStrs(x1, x2, label1, label2)
    x1 = nonnans(x1);
    x2 = nonnans(x2);
    [~, p_t] = ttest(x1, x2);
    p_w = ranksum(x1, x2);
    str1 = sprintf('%s : mean : %.2g. median: %.2g ', label1, mean(x1), median(x1));
    str2 = sprintf('%s : mean : %.2g. median: %.2g ', label2, mean(x2), median(x2));
    str3 = sprintf('p_t = %.2g. p_W = %.2g', p_t, p_w);

    stats_str = {str1, str2, str3};        
end

function showAverageFracRemoved_vs_refrPeriod(allStats_SPG_C, fig_opt)

    anyData_all = ~cellfun(@isempty, allStats_SPG_C);
    idx_gids_withData = find(any(any(anyData_all, 1),2));
    nGids_usable = length(idx_gids_withData);

    allStats_SPG_C = allStats_SPG_C(:, :, idx_gids_withData);
    
    allStats_SPG = [allStats_SPG_C{:}];
    allFracsRemoved = [allStats_SPG.allFracsRemoved{:}];
    
    
    S = getResultsStruct(allStats_C, opt, 'autoPruning');
    3;


end

function compareFracPruned_vs_fracFalsePos(allStats_SPG_auto_C, sortingFeatureSets, pruningFeatureSets, plot_opt)
    S_auto = getResultsStruct(allStats_SPG_auto_C, plot_opt, 'autoPruning');

    allStats = gatherStatsForGroups(allStats_SPG_auto_C, 'autoPruning', plot_opt.maximizeWhenChoosingGid);
        
    S_auto = getResultsStruct({allStats}, plot_opt, 'autoPruning');
%         allFracsRemoved_ci_sm = gaussSmooth(allFracsRemoved_ci, 2);
%         allFNull_t{ci} = getFNull_t(refrPeriod_range_ms, allFracsRemoved_ci);
%         maxDiffs_Fremoved_Null(ci) = max( allFracsRemoved_ci - fNull_t );
%         maxDiffs_Fremoved_Null_sm(ci) = max( allFracsRemoved_ci_sm - fNull_t );

    
%     for i = 1:length
    3;
    %%
%     pruned_lim = [.0 1];
    
%     'meanDiff'
%     'meanDiff_norm'
%     'medianDiff'
%     'medianDiff_norm'
%     'diffAtTrefrac'
%     'diffAtTrefrac_norm'
%     'diffAroundTrefrac_1ms'
%     'diffAroundTrefrac_1ms_norm'

%     plotNames = {'FracPruned', 'FalsePos'};
    plotNames = {'FracPruned', 'Change'};
%     plotNames = {'pctRefr', 'FalsePos'};
%     plotNames = {'diffAroundTrefrac_1ms_norm', 'FalsePos'};
    [data, labels] = deal(cell(1,2));
    
%     MU_fields = fieldnames(S_auto.MU_stats(1));
    for i = 1:2
        switch plotNames{i}
            case 'Change',     data{i} = -S_auto.sorting_falsePosRate_change; labels{i} = 'Effectiveness (Reduction in False Positives)';
            case 'FracPruned', data{i} = S_auto.fracOfClusterPruned; labels{i} = 'Fraction Pruned';
            case 'FalsePos',   data{i} = S_auto.fracOfClusterFalsePos; labels{i} = 'Fraction False Positives';
            case 'pctRefr',   data{i} = S_auto.pctRefrSpikePairs; labels{i} = '% refractory spikes';
            case 'nRefr',   data{i} = S_auto.nRefrSpikePairs; labels{i} = '# refractory spikes';
            case MU_fields, data{i} = [S_auto.MU_stats.(plotNames{i})]; labels{i} = plotNames{i};                
            otherwise, error('Unknown variable');
        end
    end
    
    
    
%     fracPruned_all = S_auto.fracOfClusterPruned;
%     fracFalsePos_all = S_auto.fracOfClusterFalsePos;
%     
%     idx_use = ibetween(fracPruned_all, pruned_lim);
%     fracPruned = fracPruned_all(idx_use);
%     fracFalsePos = fracFalsePos_all(idx_use);
%     %%
    differentMarkersForGroups = 0;
        
    figure(34); clf; hold on; box on;
    if ~differentMarkersForGroups
        plot(data{1}, data{2}, '.');    
    else
        Gids = S_auto.Gid;
        [uGids, Gid_idx] = uniqueList(Gids);
        cols = jet(length(uGids));
        for i = 1:length(uGids)
            plot(data{1}(Gid_idx{i}), data{2}(Gid_idx{i}), ['.'], 'color', cols(i,:), 'marker', marker(i));        
        end   
    end
    hold on;
%     plot(fracPruned, fracFalsePos, 'mo');
    p = polyfit(data{1}, data{2}, 1);
    hold on;
    fplot(@(x) polyval(p, x), xlim, 'r:');
    [cc, pval] = corr(data{1}(:), data{2}(:));
%     hold on; 
%     fplot(@(x) x, xlim, 'r:');
    xlabel(labels{1}); ylabel(labels{2});
%     title(sprintf('cc = %.2f. p = %.2g', cc, pval));

    %%
%     figure(35); clf; hold on; box
%     plot(S_auto.pctRefrSpikePairs, S_auto.fracOfClusterFalsePos,'.');
    
    
    
    3;
    %%
     
    3;



end


                
function compareTwoPruningMethods(allStats_SPG_C1, allStats_SPG_C2, pruneType1, pruneType2, opt)
    
    opt.mergeCrossPruningsFromSameGroup = 1;
%     opt.maximizeWhenChoosingGid = 'all';
    opt.maximizeWhenChoosingGid = 'nClusts';
    allStats1 = gatherStatsForGroups(allStats_SPG_C1, pruneType1, opt);
    allStats2 = gatherStatsForGroups(allStats_SPG_C2, pruneType2, opt);
        
    S1 = getResultsStruct({allStats1}, opt, pruneType1);
    S2 = getResultsStruct({allStats2}, opt, pruneType2);
    
    [common, idx_1, idx_2] = intersect(S1.Gid_sort_prune, S2.Gid_sort_prune);
    %%
    figure(88); clf; hold on; box on;
    e1 = S1.sorting_falsePosRate_change(idx_1);
    e2 = S2.sorting_falsePosRate_change(idx_2);
%     plot(e1, e2, '.');
    
    label_1 = pruneType2Label(pruneType1);
    label_2 = pruneType2Label(pruneType2);
        
    xlabel(label_1); 
    ylabel(label_2); 
    3;
    Gids = S1.Gid(idx_1);
    [uGids, Gid_idx] = uniqueList(Gids);
    cols = jet(length(uGids));
    for i = 1:length(uGids)
        plot(e1(Gid_idx{i}), e2(Gid_idx{i}), ['o'], 'color', cols(i,:), 'marker', marker(i));        
    end    
    hold on; fplot(@(x) x, xlim, 'r:')    
    axis square

    maxOverlaps = S1.clustPairOverlaps(idx_1);
    
    
    
    %%
%     S_1.Gid_sort_prune
    
    figure(2); clf;
    xx = [1, 2]; yy = [mean(e1), mean(e2)]; ee = [stderr(e1), stderr(e2)];
    bar(xx, yy); hold on;
    errorbar(xx, yy, ee, 'r.');
    set(gca, 'xticklabel', {label_1, label_2});
    xlim([.5, 2.5]);
    [~, p_t] = ttest(e1, e2);
    p_w = ranksum(e1, e2);
    title(sprintf('p_t = %.2g. p_W = %.2g', p_t, p_w))
    
    
    3;



end
    
function s = pruneType2Label(pruneType)
    switch pruneType
        case 'autoPruning',        s = 'Pruning';
        case 'crossPruning',       s= 'Cross-Pruning';
        case 'crossPruning_D',     s= 'Cross-Pruning (D)';
        case 'crossPruning_N',     s= 'Cross-Pruning (N)';
        case 'mergePruning',       s= 'Post-merge Pruning';
        case 'pruneMergePruning',  s= 'Prune-merge Pruning';
    end
end

function compareMergePruning_and_PruneMergePrune(allStats_SPG_merge_C, allStats_SPG_pruneMerge_C, plot_opt)

    plot_opt.useCombinedPruningStats = 1;
    S_mergePrune      = getResultsStruct(allStats_SPG_merge_C, plot_opt, 'mergePruning');
    S_pruneMergePrune = getResultsStruct(allStats_SPG_pruneMerge_C, plot_opt, 'pruneMergePruning');
    
    %%
    [common, idx_mergePrune, idx_pruneMergePrune] = intersect(S_mergePrune.Gid_sort_prune, S_pruneMergePrune.Gid_sort_prune);
    
    figure(90); clf; hold on; box on;
    change_mergePrune = -S_mergePrune.sorting_falsePosRate_change(idx_mergePrune);
    change_pruneMergePrune = -S_pruneMergePrune.sorting_falsePosRate_change(idx_pruneMergePrune);
%     plot(change_mergePrune, change_pruneMergePrune, '.');
    xlabel('Change in fraction of FP (Merge Pruning)'); ylabel('Change in fraction of FP (Prune Merge Pruning)');
    3;
    Gids = S_mergePrune.Gid(idx_mergePrune);
    [uGids, Gid_idx] = uniqueList(Gids);
    cols = jet(length(uGids));
    for i = 1:length(uGids)
        plot(change_mergePrune(Gid_idx{i}), change_pruneMergePrune(Gid_idx{i}), ['o'], 'color', cols(i,:), 'marker', marker(i));        
    end    
    hold on; fplot(@(x) x, xlim, 'r:')    
    axis square
    3;
    %%
    diffs = change_mergePrune(Gid_idx{i}) - change_pruneMergePrune(Gid_idx{i});
    figure(91); hist( diffs, 50);
    title(sprintf('Median : %.2f', median(diffs)));
    xlabel( 'Merge Prune - Prune-Merge-Prune');

end



function S = getResultsStruct(allStats, opt, pruningMode)
    %%
    if iscell(allStats)
        allStats = [allStats{:}];
    end               

    idx_group_use = true(1, length(allStats));

    nPruned_allClusts = [allStats.nPruned_allClusts];
    if isfield(opt, 'nPrunedTotal_min')                    
        idx_group_use = idx_group_use & (nPruned_allClusts > opt.nPrunedTotal_min);
    end
    S_use = allStats(idx_group_use);
    nGroups = length(S_use);
        
    mergePruningsFromSameGroup = isfield(opt, 'mergeCrossPruningsFromSameGroup') && opt.mergeCrossPruningsFromSameGroup && strcmp(pruningMode, 'crossPruning');
    if mergePruningsFromSameGroup
%         S_use = arrayfun(@mergeCrossPruningStatPairs, S_use);
        for i = 1:length(S_use)
            S_use(i) = mergeCrossPruningStatPairs(S_use(i), opt);
        end
        
    end
        
    %%
    if strcmp(pruningMode, 'pruneMergePruning') && isfield(opt, 'useCombinedPruningStats') && opt.useCombinedPruningStats
        fnames = fieldnames(S_use(1).autoAndMergePruningStats);
        for gi = 1:length(S_use)
            for fld_i = 1:length(fnames)
                S_use(gi).(fnames{fld_i}) = S_use(gi).autoAndMergePruningStats.(fnames{fld_i}); 
            end
        end
        3;
    end
    %%
    addGroupSortPrune_field = isfield(opt, 'doGrpSortPrune') && opt.doGrpSortPrune == 1;
    if addGroupSortPrune_field
        Group_Sort_Prune            = arrayfun(@(s) sprintf('%d_%s_%s', s.Gid, s.sortingFeatures, s.pruningFeatures), S_use, 'un', 0);
    end
    Gids                        = [S_use.Gid];
    GC                          = [S_use.GC];   
    
    regularizeFracRemoved = 1;
    if regularizeFracRemoved
        div = @(x,y) (x+1)./(y+1);
    else
        div = @(x,y) x./y;
    end
    
    falsePos_fracRemoved        = div([S_use.nFalsePosPruned], [S_use.nSpkInClust_FalsePos]);
    correctICspikes_fracRemoved = div([S_use.nICspikesPruned], [S_use.nSpkInClust_IC]);
    precision                   = falsePos_fracRemoved ./ correctICspikes_fracRemoved;
    fracPruned_falsePos         = [S_use.fracPruned_falsePos];
    fracPruned_ICspikes         = [S_use.fracPruned_ICspikes];                                
    sorting_falsePosRate_before = [S_use.sorting_falsePosRate_before];
    sorting_falsePosRate_after = [S_use.sorting_falsePosRate_after];
%     sorting_falsePosRate_change = [S_use.sorting_falsePosRate_change];
    sorting_falsePosRate_change = (sorting_falsePosRate_before - sorting_falsePosRate_after)./sorting_falsePosRate_before;
    if strcmp(pruningMode, 'autoPruning')
        fracOfClusterFalsePos       = 1-[S_use.fracOfCluster_IC];
        fracOfClusterPruned         = [S_use.nSpikesPruned]./[S_use.nSpkInClust_Tot];
        nSpkInClust_tot             = [S_use.nSpkInClust_Tot];
        nRefrSpikePairs             = [S_use.nRefrSpikePairs];
        pctRefrSpikePairs           = [S_use.pctRefrSpikePairs];
        allRefrPeriods              = [S_use.refrPeriod_ms];
        allRefrPeriod_ranges        = {S_use.refrPeriod_range_ms};
        allFracRemoved              = [S_use.allFracsRemoved];
%         allFNull                    = [S_use.allFNull_t];        
    end
    if isfield(S_use, 'clustPairOverlaps');
       clustPairOverlaps = [S_use.clustPairOverlaps];
    end
    3;

    if ~mergePruningsFromSameGroup && length(Gids) < length(falsePos_fracRemoved);
        %%
        nPerGroup = arrayfun(@(s) length(s.clustIds), S_use);
        idxs_expand_C = arrayfun(@(i, n) ones(1,n)*i, 1:nGroups, nPerGroup, 'un', 0);
        idxs_expand = [idxs_expand_C{:}];
        
        Gids = Gids(idxs_expand);
        if addGroupSortPrune_field
            Group_Sort_Prune = Group_Sort_Prune(idxs_expand);        
        end
        if strcmp(pruningMode, 'autoPruning')
            allRefrPeriod_ranges = allRefrPeriod_ranges(idxs_expand);
%             MU_stats = cellfun(@(rng, fr, fn, rp) getMUstats(rng, fr, fn, rp), allRefrPeriod_ranges, allFracRemoved, allFNull, num2cell(allRefrPeriods));            
        end
        pruningFeatures = {S_use(idxs_expand).pruningFeatures};
    end
    
    %%
    idx_keep_clust = true(1, length(falsePos_fracRemoved));
    if isfield(opt, 'nPrunedCluster_min')
        nPrunedEachCluster = [S_use.nSpikesPruned];
        idx_keep_clust = nPrunedEachCluster > opt.nPrunedCluster_min;
    end


    switch pruningMode
        case {'autoPruning', 'mergePruning'}, fracIC = [S_use.fracOfCluster_IC];
        case 'crossPruning', 
            if opt.mergeCrossPruningsFromSameGroup
                fracIC = [S_use.fracOfClusts_IC];
            else
                fracIC = min( cat(2, S_use.fracOfClustPair_IC), [],  1);
            end
    end           
        
    if strcmp(pruningMode, 'crossPruning')
        if isfield(opt, 'minFracIC')
            idx_keep_clust = idx_keep_clust & (fracIC > opt.minFracIC) & ~isnan(fracIC);
        end
    end    
    

    
    
%     S.pruningFeatures             = pruningFeatures(idx_keep_clust);
    if addGroupSortPrune_field
        S.Gid_sort_prune              = Group_Sort_Prune(idx_keep_clust);
    end       
    S.Gid                         = Gids(idx_keep_clust);
    S.GC                          = GC(idx_keep_clust);
    S.falsePos_fracRemoved        = falsePos_fracRemoved(idx_keep_clust);
    S.correctICspikes_fracRemoved = correctICspikes_fracRemoved(idx_keep_clust);
    S.fracPruned_falsePos         = fracPruned_falsePos(idx_keep_clust);
    S.fracPruned_ICspikes         = fracPruned_ICspikes(idx_keep_clust);
    S.sorting_falsePosRate_change =  sorting_falsePosRate_change(idx_keep_clust);
    if strcmp(pruningMode, 'autoPruning')
        S.fracOfClusterPruned         = fracOfClusterPruned(idx_keep_clust);
        S.fracOfClusterFalsePos       = fracOfClusterFalsePos(idx_keep_clust);  
        S.nSpkInClust_tot             = nSpkInClust_tot(idx_keep_clust);
        S.nRefrSpikePairs             = nRefrSpikePairs(idx_keep_clust);
        S.pctRefrSpikePairs           = pctRefrSpikePairs(idx_keep_clust);     
        S.precision                   = precision(idx_keep_clust);
%         S.MU_stats                    = MU_stats(idx_keep_clust);
    end
    if exist('clustPairOverlaps', 'var')
        S.clustPairOverlaps = clustPairOverlaps(idx_keep_clust);
    end
    
    
end


function MU_stats = getMUstats(refrPeriod_ranges, fracRemoved, fNull, t_refrac)
    %%
    if ~isnan(fracRemoved(1)) && length(fracRemoved) > 1
        idx_allowed = refrPeriod_ranges >= 0.8;
            fracRemoved = fracRemoved(idx_allowed);
            refrPeriod_ranges = refrPeriod_ranges(idx_allowed);
            fNull = fNull(idx_allowed);

        idx_trefrac = indmin(abs(refrPeriod_ranges-t_refrac));

        idx_tefrac_1ms = find(abs(refrPeriod_ranges-t_refrac) < 1);
    else
        [idx_trefrac, idx_tefrac_1ms] = deal(1);        
    end
    fracRemoved_norm = fracRemoved/ fracRemoved(end);
    fNull_norm = fNull / fNull(end);
    
    MU_stats.meanDiff = mean(fracRemoved - fNull);
    MU_stats.meanDiff_norm = mean(fracRemoved_norm - fNull_norm);
    MU_stats.medianDiff = median(fracRemoved - fNull);
    MU_stats.medianDiff_norm = median(fracRemoved_norm - fNull_norm);
    MU_stats.diffAtTrefrac = fracRemoved(idx_trefrac) - fNull(idx_trefrac);
    MU_stats.diffAtTrefrac_norm = fracRemoved_norm(idx_trefrac) - fNull_norm(idx_trefrac);
    MU_stats.diffAroundTrefrac_1ms = mean( fracRemoved(idx_tefrac_1ms) - fNull(idx_tefrac_1ms) );
    MU_stats.diffAroundTrefrac_1ms_norm = mean( fracRemoved_norm(idx_tefrac_1ms) - fNull_norm(idx_tefrac_1ms) );


end

function idx = search(data, items)
    idx = zeros(size(items));
    for i = 1:length(items)
        idx(i) = find(data == items(i), 1);                
    end
end

function s = mergeCrossPruningStatPairs(s_orig, opt)
%     search = @(data, x) arrayfun(@(xi) find(data == xi, 1), x);
    %%
    s = s_orig;
            
    pair_ok = true(1, length(s.GC));
    if isfield(opt, 'minFracIC')
%         if isfield(s, 'fracOfClustPair_IC')
            pair_ok = pair_ok & min(s.fracOfClustPair_IC, [], 1) > opt.minFracIC;
%         elseif isfield(s, 'fracOfCluster_IC')
%             pair_ok = pair_ok & s.fracOfCluster_IC > opt.minFracIC;
%         end        
    end
    if isfield(opt, 'nPrunedCluster_min')
        pair_ok = pair_ok & s.nSpikesPruned >= opt.nPrunedCluster_min;
    end
    idx_pair_use = find( pair_ok );

    [uClustId_use] = unique(s.clustIdPairs(:,idx_pair_use) );
    first_idx_clust_use = search(s.clustIdPairs(:), uClustId_use);

    s.clustIds = uClustId_use(:)';    
    s.clustIdPairs = s.clustIdPairs(:,idx_pair_use);
    s.clustPairOverlaps = max(s.clustPairOverlaps(:,idx_pair_use) );
        if isempty(s.clustPairOverlaps), s.clustPairOverlaps = nan; end    
    s.GC = s.GC(idx_pair_use);

    s.nSpkInClusts_Tot = sum( s.nSpkInClustPair_Tot(first_idx_clust_use) );
    s.nSpkInClusts_IC = sum( s.nSpkInClustPair_IC(first_idx_clust_use) );
    s.nSpkInClusts_falsePos = sum( s.nSpkInClustPair_falsePos(first_idx_clust_use) );

    s.nSpkInClustPair_Tot = s.nSpkInClustPair_Tot(:, idx_pair_use);
    s.nSpkInClustPair_IC = s.nSpkInClustPair_IC(:, idx_pair_use);
    s.nSpkInClustPair_falsePos = s.nSpkInClustPair_falsePos(:, idx_pair_use);
    s.fracOfClustPair_IC = s.fracOfClustPair_IC(:, idx_pair_use);
    s.fracOfClustPair_falsePos = s.fracOfClustPair_falsePos(:, idx_pair_use);

    s.fracOfClusts_IC = s.nSpkInClusts_IC / s.nSpkInClusts_Tot;
    s.fracOfClusts_falsePos = s.nSpkInClusts_falsePos / s.nSpkInClusts_Tot;

    s.nSpikesPruned = sum(s.nSpikesPruned(idx_pair_use));
    s.nFalsePosPruned = sum(s.nFalsePosPruned(idx_pair_use));
    s.nICspikesPruned = sum(s.nICspikesPruned(idx_pair_use));
    s.nPruned_allClusts = s.nSpikesPruned;

    s.fracPruned_falsePos = s.nFalsePosPruned  / s.nSpikesPruned;
    s.fracPruned_ICspikes = s.nICspikesPruned  / s.nSpikesPruned;
    s.falsePos_fracRemoved = s.nFalsePosPruned / s.nSpkInClusts_falsePos;
    s.correctICspikes_fracRemoved = s.nICspikesPruned / s.nSpkInClusts_IC;
    s.falsePos_fracRemaining = 1 - s.falsePos_fracRemoved;
    s.correctICspikes_fracRemaining = 1 - s.correctICspikes_fracRemoved;
    s.sorting_falsePosRate_before = s.nSpkInClusts_falsePos / s.nSpkInClusts_Tot;
    s.sorting_falsePosRate_after = (s.nSpkInClusts_falsePos - s.nFalsePosPruned) / (s.nSpkInClusts_Tot - s.nSpikesPruned);
    s.sorting_falsePosRate_change = s.sorting_falsePosRate_after - s.sorting_falsePosRate_before;

    s.refrPeriod_ms = min(s.refrPeriod_ms(idx_pair_use) );
end


function allStats = gatherStatsForGroups(allStats_SPG_C, pruningMode, groupSelectionMode)
                        
    % take allStats_SPG_C, which is nSortFet x nPruneFet x nGids and gather all the data for each
    % group into a cell array 1 x nGids.
    % For a particular group, either collect all data or select based on nClusts / nPruned.

    anyData_all = ~cellfun(@isempty, allStats_SPG_C);
    idx_gids_withData = find(any(any(anyData_all, 1),2));
    nGids_usable = length(idx_gids_withData);
    %         Gids_usable = Gids(idx_gids_withData);

%     allStats_SPG_C = allStats_SPG_C(:, idx_gids_withData);
%     anyDataForGidsAndSets = anyData_all(:,idx_gids_withData);

    allStats_used_C = cell(1,nGids_usable);
    

%             nPrunedTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpikesPruned'});
%             nSpkTot_field = switchh(pruningMode, {'autoPruning', 'crossPruning'}, {'nPrunedTotal', 'nSpkInClusts_Tot'});
            
            %%
    for gi = 1:nGids_usable
        gid_idx = idx_gids_withData(gi);
        [idx_sort_fet, idx_prune_fet] = find(anyData_all(:,:, gid_idx));
        idx_usable_gi = sub2indV(size(anyData_all), [idx_sort_fet, idx_prune_fet, gid_idx*ones(size(idx_sort_fet))]);
        allStats_usable_gi = allStats_SPG_C(idx_usable_gi)';
        if length(idx_usable_gi) > 1
            3;
            if any(strcmp(pruningMode, {'autoPruning', 'mergePruning'}))
                nClusts = cellfun(@(s) length(s), allStats_usable_gi);
                nPrunedTotal = cellfun(@(s) s.nPruned_allClusts, allStats_usable_gi);
                nInClustsTotal = cellfun(@(s) sum(s.nSpkInClust_Tot), allStats_usable_gi);
                fracPrunedTotal = nPrunedTotal ./ nInClustsTotal;                                        

            elseif strncmp(pruningMode, 'crossPruning', 10)

                nClusts = cellfun(@(s) length(s.GC), allStats_usable_gi);                    
                nPrunedTotal = cellfun(@(s) sum(s.nSpikesPruned), allStats_usable_gi);
                nInClustsTotal = cellfun(@(s) sum(s.nSpkInClusts_Tot), allStats_usable_gi);
                fracPrunedTotal = nPrunedTotal ./ nInClustsTotal;                                        

            end

            A = [nClusts(:), fracPrunedTotal(:), nPrunedTotal(:)];                

            if isstruct(groupSelectionMode)
                groupSelectionMode = groupSelectionMode.maximizeWhenChoosingGid;
            end
            
            if ~strcmp(groupSelectionMode, 'all')
                switch groupSelectionMode
                    case 'nClusts', column_order = [1 2 3];
                    case 'nPruned', column_order = [2 1 3];
                end                                
                [A_sorted, sort_idx] = sortrows(A, column_order);

                idx_use_now = sort_idx(1); % pick best match
            else
                idx_use_now = 1:length(allStats_usable_gi);                    
            end
            allStats_this_gid  = [allStats_usable_gi{idx_use_now}];
        else
            
            allStats_this_gid = [allStats_usable_gi{:}];
        end
%         idx_use(gi) = idx_use_now(1); %            
        
        
        allStats_used_C{gi} = {allStats_this_gid}; %[allStats_usable_gi(idx_use_now)];
        assert( all([allStats_this_gid.Gid] == allStats_this_gid(1).Gid))
    end
    
    allStats_C = [allStats_used_C{:}];
    allStats = [allStats_C{:}];
        
end

function [medianRatio, meanE] = plot_frTP_vs_frFP(allStats_SPG_C_all, pruningMode, useICcellRefrPeriods, sortingFeatureSets_all, pruningFeatureSets_all,  plot_opt)
%         if strcmp(pruningMode
        3;        
            
%         allStats_SPG_C = allStats_SPG_refr_C{refr_type_i};
%         groupSelectionMode = plot_opt.maximizeWhenChoosingGid;
%         groupSelectionMode = 'all';
        groupSelectionMode = 'nClusts';        
        
%         sortingFeatureSets_all = sortingFeatureSets;
%         pruningFeatureSets_all = pruningFeatureSets;
%         allStats_SPG_C_all = allStats_SPG_C;
        
        %%
        doPlots = nargout == 0;
        
        allStats_SPG_C = allStats_SPG_C_all;
        if isfield(plot_opt, 'sortFet')
            sort_fet_idx = find(strcmp(plot_opt.sortFet, sortingFeatureSets_all));
            assert(length(sort_fet_idx) == 1);
            
            allStats_SPG_C = allStats_SPG_C(sort_fet_idx, :, :);
            sortingFeatureSets = sortingFeatureSets_all(sort_fet_idx);
        end
        if isfield(plot_opt, 'pruneFet')
            prune_fet_idx = find(strcmp(plot_opt.pruneFet, pruningFeatureSets_all));
            assert(length(prune_fet_idx) == 1);
            
            allStats_SPG_C = allStats_SPG_C(:, prune_fet_idx, :);
            pruningFeatureSets = pruningFeatureSets_all(prune_fet_idx);
        end
        
%         allStats = gatherStatsForGroups(allStats_SPG_C, pruningMode, 'all');
        allStats = gatherStatsForGroups(allStats_SPG_C, pruningMode, groupSelectionMode);
        
        S = getResultsStruct(allStats_SPG_C, plot_opt, pruningMode);
%         S2 = getResultsStruct(allStats, plot_opt, pruningMode);
        
%         allStats = gatherStatsForGroups(allStats_SPG_C, pruningMode, 'nPruned');

        %%  PLOTTING
        
%         usedGids = cell(1, nSortingSets);
        nSortingSets_here = 1; % nSortingSets
%         nRefrTypes = nnz(refrPeriodTypes);
%         overlapPlots = nRefrTypes == 2;
        
        markerStyle = 'b.';
        
        cols = [0 0 1]; %jet(nSortingSets);
        
        doScatterPlot = 1;
        doHistFracPrunedFP = 0 && 1;
        doHistChangeFP = 1;
        
        
        figId_scatter = plot_opt.figId+1;
        figId_fracFP = plot_opt.figId+2;
        figId_changeFP = plot_opt.figId+3;
        
%         firstPlot = ~overlapPlots ||  ((refr_type_i == 1) || (nRefrTypes == 1));
        if doPlots
            if doScatterPlot, figure(figId_scatter);  hold on; box on; end
            if doHistFracPrunedFP, figure(figId_fracFP);  hold on; box on; end
            if doHistChangeFP, figure(figId_changeFP);  hold on; box on; end

            firstPlot = 1;

            if firstPlot && doScatterPlot, clf(figId_scatter); end
            if firstPlot && doHistFracPrunedFP, clf(figId_fracFP); end
            if firstPlot && doHistChangeFP, clf(figId_changeFP); end
        end
%         else
%             markerStyle = 'go';
%         end            

                  

        nPrunedTotal_min = 2;

        idx_group_use = [allStats.nPruned_allClusts] >= nPrunedTotal_min;
        allStats = allStats(idx_group_use);

        if isempty(allStats)
            3;
        end
%             usedGids{fet_idx} = [allStats.Gid];
%%

        falsePos_fracRemoved        = S.falsePos_fracRemoved;
        correctICspikes_fracRemoved = S.correctICspikes_fracRemoved;

        fracPruned_falsePos = S.fracPruned_falsePos;
        fracPruned_ICspikes = S.fracPruned_ICspikes;


        sorting_falsePosRate_change = S.sorting_falsePosRate_change;

        if nSortingSets_here > 1
            col = cols(fet_idx, :);
            mk = marker(fet_idx);
            markerStyle = [mk col];
        else
            markerStyle = 'b.';
        end
        %%
        

        allP = falsePos_fracRemoved ./ correctICspikes_fracRemoved;
        medianRatio = nanmedian(allP);
        if doScatterPlot && doPlots
            %%
            
            figure(figId_scatter); clf; hold on; box on;
            plot(correctICspikes_fracRemoved, falsePos_fracRemoved, markerStyle);
            L = max([correctICspikes_fracRemoved, falsePos_fracRemoved])*1.05;            

            h1 = gca;
            
            
%             ismedianRatioInf = isinf(medianRatio);
%             if ismedianRatioInf
%                 medianRatio = max(allP(isfinite(allP)));
%             end
            %%            
%             L = 0.4;
                axis([0 L, 0 L]); 
                fplot(@(x) x, [0 L], 'r:')
                xlabel('frac reduction of true pos (frTP)');
                ylabel('frac reduction of false pos (frFP)');
                pruningMode_str = switchh(pruningMode, {'autoPruning', 'crossPruning', 'mergePruning', 'pruneMergePruning'}, {'Cluster Pruning', 'Cross-pruning', 'Post-merge Pruning', 'Post-merge Pruning (net)'});
                ideal_str = iff(useICcellRefrPeriods, '  (using t_{true})', '  (using estimated t_{refrac})');
    %             pruningFet_str = sprintf('Pruning Features: %s', pruningFeatures);
                pruningFet_str = '';
                str1 = [pruningMode_str ideal_str];
    %             s_brack = iff(ismedianRatioInf, '[]', '  ');
                title({str1, sprintf('Median Precision (P): %.1f', medianRatio )});            
        end
        
        
        meanE = nanmean(sorting_falsePosRate_change);        
        if  doHistChangeFP && doPlots
            doPercentage = 1;
%             f = iff(doPercentage, 100, 1);
%             pct_txt = iff(doPercentage, 'Percentage', 'Fractional');
            if doPercentage
                f = 100;
                fmt = 'Mean: %.1f%%. Range : [%.1f%% -- %.1f%%]';
                pct_txt = 'Percentage';
            else
                f = 1;
                fmt = 'Mean: %.3f. Range : [%.2f -- %.2f]';
                pct_txt = 'Fractional';
            end
            
            figure(figId_changeFP); clf;            
            nBins = 30;
            L_change = lims(sorting_falsePosRate_change*f, .01); db = diff(L_change)/nBins;            
            binE = roundToNearest(L_change(1), db, 'down') : db : roundToNearest(L_change(2), db, 'up');
            
            hist(sorting_falsePosRate_change*f, binEdge2cent( binE ));
            title({str1, sprintf(fmt, ...
                nanmean(sorting_falsePosRate_change)*f, prctile( sorting_falsePosRate_change*f, [5, 95]))});
            3;
            xlims = lims(sorting_falsePosRate_change*f, .1);
            xlim(xlims);
            xlabel( sprintf('%s Decrease in FP (E)', pct_txt));
            ylabel(' Number of clusters ');
            h3 = gca;
        
        end
        
        3;
%         figure(654); clf;
%         plot(fracPruned_ICspikes, fracPruned_falsePos, 'o');
%         L = max([fracPruned_ICspikes, fracPruned_falsePos])*1.05;
%         axis equal;
%         axis([0 L 0 L]);
%         hold on;
%         fplot(@(x) x, [0 L], 'r:')
%         xlabel('frac removed - correct');
%         ylabel('frac removed - false pos');
        3;

        if doHistFracPrunedFP
            figure(figId_fracFP); clf;
%             
%             plot(fracPruned_ICspikes, fracPruned_falsePos, mk, 'color', col);
%             L = max([fracPruned_ICspikes, fracPruned_falsePos])*1.05;
% %             axis equal;
%             axis([0 inf 0 inf]);
% 
% %             hold on;
%             fplot(@(x) x, [0 L], 'r:')
%             xlabel('frac removed - correct');
%             ylabel('frac removed - false pos');
            hist(fracPruned_falsePos, 20);
            title({str1, sprintf('Pruned spikes: False Pos/Total (mean: %.2f)', nanmean(fracPruned_falsePos) )})
            h2(refr_type_i) = gca;

%         figure(655); clf;
%         hist([allStats.sorting_falsePosRate_change], 50);

        end

        %%
%         figure(201); hist(fracPruned_falsePos, 50); xlabel('Fraction of Pruned spikes that were false positives')
%         figure(202); hist(fracFalsePos_pruned, 50); xlabel('Fraction of False-positives that were pruned');

        
                3;
           
        
        
    
    
end
    

function compareRefrPeriod_auto_vs_IC
    
    if nRefrTypes == 2
        matchAxes('xy', h1);        
        matchAxes('xy', h2);
        matchAxes('xy', h3);
        
        %%
        allStats_auto = [allStatsForFet_C_refr{1}{:}];
        allStats_IC = [allStatsForFet_C_refr{2}{:}];
        
        GC_auto = [allStats_auto.GC];
        GC_IC = [allStats_IC.GC];
        [GC_common, idx_auto, idx_ic] = intersect(GC_auto, GC_IC);

        
        %%
        minFracIC_use = 0.6;
        figure(680); clf; hold on; box on;
        fracRemoved = [allStats_auto.correctICspikes_fracRemoved] + [allStats_auto.falsePos_fracRemoved];
        fracRemoved = fracRemoved(idx_auto);
        refrPeriod_auto = [allStats_auto.refrPeriod_ms];  refrPeriod_auto = refrPeriod_auto(idx_auto);
        if strcmp(pruningMode, 'autoPruning')
            refrPeriod_IC = [allStats_IC.refrPeriod_ms_ideal_allEC];  refrPeriod_IC = refrPeriod_IC(idx_ic);
        else
            refrPeriod_IC = [allStats_IC.refrPeriod_ms];  refrPeriod_IC = refrPeriod_IC(idx_ic);
        end
            
        if strcmp(pruningMode, 'autoPruning')
            fracIC = [allStats_auto.fracOfCluster_IC]; fracIC = fracIC(idx_auto);
        else
            fracIC = min([allStats_auto.fracOfClustPair_IC], [], 1); fracIC = fracIC(idx_auto);
        end
%         fracIC_IC = [allStats_IC.fracOfCluster_IC]; fracIC_IC = fracIC_IC(idx_ic);
        
        idx_use = fracIC > minFracIC_use & fracRemoved > .2;
        fracRemoved = fracRemoved(idx_use);
        refrPeriod_auto = refrPeriod_auto(idx_use);
        refrPeriod_IC = refrPeriod_IC(idx_use);

        L = max([refrPeriod_auto, refrPeriod_IC])*1.05;
%         L = 8;
        plot(refrPeriod_IC, refrPeriod_auto, '.');
        
        cc = spearmanRho(refrPeriod_IC, refrPeriod_auto);
        xlabel('Refr period (IC)');
        ylabel('Refr period (auto)'); 
        axis([0 L 0 L]);
        hold on;
        fplot(@(x) x, [0 L], 'r:');
        title(sprintf('cc = %.2f', cc))
        
%         figure(
        3;
    end
end

function idx_each = NtoIdx(N_each)
    %%
    idx_each = zeros(1, sum(N_each));        
    cum_idx = [0, cumsum(N_each)];
    for i = 1:length(N_each)
       idx_each(cum_idx(i)+1 : cum_idx(i+1)) = i;        
    end


end


function printGroupStats(allStats_SPG_auto_C)
    3;


end







%                         if (prunedStats.nPruned > 0)
%                     falsePos_fracRemoved(i) = prunedStats.falsePos_fracRemoved;
%                     correctICspikes_fracRemoved(i) = prunedStats.correctICspikes_fracRemoved;
%                     
%                     fracPruned_falsePos(i) = prunedStats.fracPruned_falsePos;
%                     fracPruned_ICspikes(i) = prunedStats.fracPruned_ICspikes;
%                 end                
%                 3;


%      'Egy_PCAsnr1'    'GLFcuw12'    'GLFcuw16'    'GLFcuw3'    'GLFcuw32'
%     'GLFcuw4'    'GLFcuw48'    'GLFcuw5'    'GLFcuw6'    'GLFcuw8'    'GLFsnr3'    'Neg'
%     'PCAcur1'    'PCAcur10'    'PCAcur12'    'PCAcur2'    'PCAcur3'    'PCAcur4'
%     'PCAcur5'    'PCAcur6'    'PCAcur7'    'PCAcur8'    'PCAcuw1'    'PCAcuw12'
%     'PCAcuw2'    'PCAcuw3'    'PCAcuw4'    'PCAcuw5'    'PCAcuw6'    'PCAcuw7'
%     'PCAcuw8'    'PCAsnr3'    'PCAsur1'    'PCAsur2'    'PCAsur3'    'PCAsur4'
%     'PCAsuw3'    '~old cat sorting'



%{

    allmeanRatios      = nan( nNPruneFets, nMaxPerPruneFet, nSortUse);
    allIsmeanRatiosInf = nan( nNPruneFets, nMaxPerPruneFet, nSortUse);
    allMedianYields            = nan( nNPruneFets, nMaxPerPruneFet, nSortUse);
    allMedianFPchanges         = nan( nNPruneFets, nMaxPerPruneFet, nSortUse);

    allP      = cell( nNPruneFets, nMaxPerPruneFet, nSortUse);
    allE   = cell( nNPruneFets, nMaxPerPruneFet, nSortUse);
    allYields      = cell( nNPruneFets, nMaxPerPruneFet, nSortUse);
    
        
    for ni = 1:nNPruneFets
        for pf_i = 1:length(pruneFetIdxs{ni})

            pf_idx = pruneFetIdxs{ni}(pf_i);
            pruneFet_str = pruningFeatureSets{pf_idx};
            
            for sf_i = 1:nSortUse
                % averaging over Gids (weird, b/c this is the independent dimension).
                allStats_C = [allStats_SPG_C{sf_i, pf_idx, :}];    
                if isempty(allStats_C)
                    continue;
                end            
                resultStruct = getResultsStruct(allStats_C, opt, pruningMode);            
                [meanRatio, medianFPchange, yield, ismeanRatiosInf, allP, allE, allY] = getRatios_FP_TP_changes(resultStruct);
                allmeanRatios     (ni, pf_i, sf_i) = meanRatio;
                allMedianYields           (ni, pf_i, sf_i) = yield;
                allIsmeanRatiosInf(ni, pf_i, sf_i) = ismeanRatiosInf;
                allMedianFPchanges  (ni, pf_i, sf_i) = medianFPchange;
                
                allP   {ni, pf_i, sf_i} = allP;
                allE{ni, pf_i, sf_i} = allE;
                allYields   {ni, pf_i, sf_i} = allY;                
                
            end
        end
    end
    %%
    if nSortUse > 1
        allmeanRatios_vN = reshape(allmeanRatios, [nNPruneFets, maxN*nSortUse, ]);
        allMedianFPchanges_vN = reshape(allMedianFPchanges, [nNPruneFets, maxN*nSortUse]);
        allMedianYields_vN = reshape(allMedianYields, [nNPruneFets, maxN*nSortUse]);
                
    end
    %%
    dim = 2;
    figure(903); clf;
    subplot(2,1,1);
    errorbar(all_Ns, nanmean(allmeanRatios, dim), nanstd(allmeanRatios, [], dim), 'b.-');
    set(gca, 'xtick', [4:4:16]);
    xlabel('# Pruning Features'); ylabel('Median Precision (P)');
    
    subplot(2,1,2);
    errorbar(all_Ns, nanmean(allMedianFPchanges, dim), nanstd(allMedianFPchanges, [], dim), 'b.-');
    set(gca, 'xtick', [4:4:16]);
    xlabel('# Pruning Features'); ylabel('Mean Effectiveness (E)');
    3;
%     a = nanmean(allmeanRatios, dim);
%     b = nanmean(allMedianFPchanges, dim);
%     c = nanmean(allMedianYields, dim);
    
    figure(904); clf;
    errorbar(all_Ns, nanmean(allMedianYields, dim), nanstd(allMedianYields, [], dim), 'b.-');
    set(gca, 'xtick', [4:4:16]);
    xlabel('# Pruning Features'); ylabel('Mean Yield (Y)');
    
%}