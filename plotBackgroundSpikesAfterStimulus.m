function plotBackgroundSpikesAfterStimulus(Gid_sets, Gid_set_labels)

    if nargin == 0

    %     sd = siteDataFor( getAllGids );
        sd = siteDataFor( getAllGids('fd', 'Gid') );
        locData = [sd.locationData];
        locIds = [locData.LocId];

        idx_f = find(strncmp({sd.stimType}, 'Movie:Flashed', 10));
        idx_o = find(strncmp({sd.stimType}, 'Grating:Orientation', 15));
        idx_s = find(strncmp({sd.stimType}, 'Grating:Spatial', 15));

        f_locIds = locIds(idx_f);
        o_locIds = locIds(idx_o);
        s_locIds = locIds(idx_s);
        allTypes = 'fos';
        allTypeLabels = {'Flashed', 'Orientation', 'Spatial Freq'};
        allLocIds = {f_locIds, o_locIds, s_locIds};
        allIdxs = {idx_f, idx_o, idx_s};

        matchTypes = 0;
        doTypes = 'fos';
        
        if length(doTypes) == 3        
            if matchTypes
                [cLocIds, ci_f, ci_o, ci_s] = intersect3(f_locIds, o_locIds, s_locIds);
                idx_f = idx_f(ci_f);
                idx_o = idx_o(ci_o);
                idx_s = idx_s(ci_s);
            end

            Gid_sets = {[sd(idx_f).Gid], [sd(idx_o).Gid], [sd(idx_s).Gid]};
            Gid_set_labels = allTypeLabels;
        elseif length(doTypes) == 2
            id1 = find(doTypes(1) == allTypes, 1);
            id2 = find(doTypes(2) == allTypes, 1);
            if matchTypes
                [cLocIds, ci_1, ci_2] = intersect(allLocIds{id1}, allLocIds{id2});
                idx1 = allIdxs{id1}(ci_1);
                idx2 = allIdxs{id2}(ci_2);
            else
                idx1 = allIdxs{id1};
                idx2 = allIdxs{id2};                
            end
            Gid_sets = {[sd(idx1).Gid], [sd(idx2).Gid]};
            Gid_set_labels = allTypeLabels([id1, id2]);

        end

            
    end
    
    if ~iscell(Gid_sets)
        Gid_sets = {Gid_sets};
    end
    if ~iscell(Gid_set_labels)
        Gid_set_labels = {Gid_set_labels};
    end    
    
    t_max_sec = 2;
    binSize_ms = 4 + 1/6;
%     binSize_sec = binSize_ms/1000;
    nAv = 10;
    bckgSkip_ms = 250;
    
    doRawPlot = 0;
    
    doAvPlot = 1;
    doAvCountPlot = 1;
    
    allOnSamePlot = 1;
    if allOnSamePlot
        fig_offsets = zeros(1, length(Gid_sets));
    else
        fig_offsets = 10*[0:length(Gid_sets)-1];
    end
    
    fig2_tit = cell(1, length(Gid_sets));
    figure(2); clf;
    figure(3); clf;
    indivMeanRates = cell(1, length(Gid_sets));
        
    for j = 1:length(Gid_sets)
        Gids = Gid_sets{j};
        Set_label = Gid_set_labels{j};
    
        [binSpks, binCount, meanRate_all, meanRate_cut, indivMeanRates{j}] = getBinnedData(Gids, t_max_sec, bckgSkip_ms);
        nBins = length(binSpks);
%         meanRate = (sum(binSpks) / sum(binCount))/binSize_sec;
        
        % 1. Just raw plots
        binMeanSpikes = binSpks ./ binCount;        
        binSize_sec = binSize_ms/1000;

        if doRawPlot
            figure(1+fig_offsets(j));
            plot([0:nBins-1]*binSize_ms, binMeanSpikes/binSize_sec, color_s(j));
            hold on;
            title(sprintf('%s,  %.2f (all), (%.2f in view)', Set_label, meanRate_all, meanRate_cut));
        end
        
        if doAvPlot
%             binSize_sec_av = binSize_sec*nAv;
            idxs = arrayfun(@(i) (i-1)*nAv+1 : i*nAv, 1:floor(nBins/nAv), 'un', 0);
            nBins_av = length(idxs);
            binMeanSpikes_av = cellfun(@(i) mean(binMeanSpikes(i)), idxs);
            binCount_av = cellfun(@(i) mean(binCount(i)), idxs)/binSize_sec;

            figure(2+fig_offsets(j));            
            binSize_av_ms = binSize_ms*nAv;
            plot([0:nBins_av-1]*binSize_av_ms, binMeanSpikes_av/binSize_sec, [color_s(j) '.-'])
            tit_str = sprintf('%s, %.2f (all), (%.2f in view)', Set_label, meanRate_all, meanRate_cut);
            if allOnSamePlot
                fig2_tit{j} = tit_str;
                title(fig2_tit);
                hold on;
            else
                title(tit_str);
            end
            
            

            3;
            if doAvCountPlot
                figure(3+fig_offsets(j));
                plot([0:nBins_av-1]*binSize_av_ms, binCount_av, [color_s(j) '.-']);
                hold on;
                set(gca, 'yscale', 'log');
                title(Set_label);
            end

        end
        
    end
    figure(2);
    legend(Gid_set_labels, 'location', 'bestOutside');
    
    figure(5); clf;
    m = cellfun(@mean, indivMeanRates);
    s = cellfun(@stderr, indivMeanRates);
    M = diag(m);
    x = 1:length(Gid_sets);
    hb = bar(x, M, 'stacked');
    hold on;
    errorbar(x, m, s, 'ks');
    set(gca, 'xtick', x, 'xticklabel', Gid_set_labels)
    s2 = [];
    for i = 1:length(Gid_sets)
        set(hb(i), 'facecolor', switchh(doTypes(i), {'f', 'o', 's'}, {'b', 'g', 'r'}));
        s2 = [s2, Gid_set_labels{i} ': ' num2str(m(i), '%.2f') ' Hz; '];
    end   
    x_e = binCent2edge(x);
    xlim(x_e([1, end]));    
    
    allNs = cellfun(@length, Gid_sets);
    if length(unique(allNs)) == 1
        s1 = sprintf('Paired %s experiments (N = %d)', cellstr2csslist(Gid_set_labels), allNs(1)) ;
    else
        s = [Gid_set_labels; num2cell(allNs)]';
        s1 = ['All ' sprintf('%s (N = %d); ', s{:}) ];
    end
        
    title({s1, s2}) ;
    3;
    
    
end


function [binSpks, binCount, meanRate_all, meanRate_cut, indivMeanRates] = getBinnedData(Gids, t_max_sec, bckgSkip_ms)
    binSize_ms = 4 + 1/6;
    binSize_sec = binSize_ms / 1000;
    
    t_max_nbins = round(t_max_sec/binSize_sec);
    bckgSkip_nbins = round(bckgSkip_ms/binSize_ms);
    
    inclExpStart = 1;   % perhaps ignore first interval (long, and doesn't follow stimulus.
    inclPostTrials = 0;
    
    nG = length(Gids);
    allGrpBckgData = cell(1,nG);
    nMaxIntervalBins = zeros(1,nG);
    nCellsPerGrp = zeros(1,nG);

    useAllSpikeMU = 0;
    includeSmallSpikes = 0;
    
    progressBar('init-', length(Gids), 40);
    for gi = 1:length(Gids)
        Gid = Gids(gi);                
                
        sd = siteDataFor('Gid', Gid, 1);
        cellIds = sd.cellIds;
        if nnz(cellIds >0) == 0
            continue;
        end
        
        
%         if useAllSpikeMU
%             grpBckgData_MU = getBackgroundSpikes(Gid, 100);            
%             grpBckgData_MU = cellfun(@minus, grpBckgData_MU, getBackgroundSpikes(Gid, -1), 'un', 0);
%             if ~includeSmallSpikes
%                 grpBckgData0 = getBackgroundSpikes(Gid, 0);                
%                 grpBckgData_MU = cellfun(@minus, grpBckgData_MU, grpBckgData0, 'un', 0);                
%             end                        
%         else
            if ~includeSmallSpikes
                cellIds = cellIds(cellIds>0);
            end
            
            grpBckgData_INDIV = getBackgroundSpikes(Gid, cellIds(1));
            for i = 2:length(cellIds)
                grpBckgData_i = getBackgroundSpikes(Gid, cellIds(i));                
                grpBckgData_INDIV = cellfun(@plus, grpBckgData_INDIV, grpBckgData_i, 'un', 0);
            end
            
%             assert(isequal(grpBckgData_MU, grpBckgData_INDIV))
            grpBckgData = grpBckgData_INDIV;
%         end
    
        nCellsPerGrp(gi) = length(cellIds);
        interval_ok = true(1, length(grpBckgData));
        nBinsEachInt = cellfun(@length, grpBckgData);
        
        if ~inclExpStart
            interval_ok(1) = 0;
        end
        if ~inclPostTrials
            interval_ok(2:end) = 0;
        end
            
        if ~inclPostTrials && inclExpStart
            nMaxInterStim = nBinsEachInt(1);
        elseif inclPostTrials
%             if length(nBinsEachInt) > 2 % have at least 1 inter-trial interval
%                 nMaxInterStim = max(nBinsEachInt(2:end-1));                        
%             else 
%                 nMaxInterStim = 0; % don't want to drag out end of experiment.
%             end
            nMaxInterStim = max(nBinsEachInt);                        
        end
        nMaxIntervalBins(gi) = nMaxInterStim;
            
        grpBckgData = grpBckgData(interval_ok);  
        allGrpBckgData{gi} = grpBckgData;
        
        progressBar(gi);
        
%         nBinsEachInt = cellfun(@length, grpBckgData);
        
        
%         progressBar(gi);
    end
    progressBar('done');
    
    
    nGlobMaxBins = max(nMaxIntervalBins);
    binSpks  = zeros(1, nGlobMaxBins);
    binCount = zeros(1, nGlobMaxBins);
    
    
    indivMeanRates  = zeros(1,nG);
    
    for gi = 1:nG
        grpBckgData = allGrpBckgData{gi};
        nInt = length( grpBckgData );
        grpBinSpks = zeros(1, nGlobMaxBins);
        grpBinCount = zeros(1, nGlobMaxBins);
        
        for i = 1:nInt
            nBins = min(length(grpBckgData{i}), nGlobMaxBins)-bckgSkip_nbins;
            grpBinSpks(1:nBins) = grpBinSpks(1:nBins) + grpBckgData{i}([1:nBins]+bckgSkip_nbins);
            grpBinCount(1:nBins) = grpBinCount(1:nBins) + nCellsPerGrp(gi);
        end               
        indivMeanRates(gi) = (sum(grpBinSpks)/sum(grpBinCount))/binSize_sec;
        
        binSpks = binSpks + grpBinSpks;
        binCount = binCount + grpBinCount;        
    end
            
    indivMeanRates = nonnans(indivMeanRates);
    meanRate_all = (sum(binSpks)/sum(binCount))/binSize_sec;
    
    if t_max_nbins < nGlobMaxBins
        binSpks = binSpks(1:t_max_nbins);
        binCount = binCount(1:t_max_nbins);
    end
    
    meanRate_cut = (sum(binSpks)/sum(binCount))/binSize_sec;

end



    
%         sd = siteDataFor('Gid', Gid, 1);
%         cellIds = sd.cellIds;
%         nC = length(cellIds);
%                 
%         grpBckgData = getBackgroundSpikes(Gid, cellIds(1));
%         for ci = 2:nC
%             grpBckgData_i  = getBackgroundSpikes(Gid, cellIds(ci));            
%             grpBckgData = cellfun(@plus, grpBckgData, grpBckgData_i, 'un', 0);                        
%         end
