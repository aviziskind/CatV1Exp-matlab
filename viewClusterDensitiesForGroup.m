function viewClusterDensitiesForGroup(Gid)
    sortingFeatureSets = {'Neg', 'GLFcuw4', 'GLFcur4'};
    nSets = length(sortingFeatureSets);

    tfToOnOff = @(tf) iff(tf, 'On', 'Off');    
    clustGrouping = 'clusters';
    [uClustIds, clustSpkIdxs] = getCellSorting(Gid, clustGrouping);
    minFracIC = 0.5;
    
    [features, featureLabels] = getGroupFeatures(Gid, sortingFeatureSets, 1, 1);
            
    IC_stats = identifyCellfromIC(Gid, struct('clustGrouping', clustGrouping));
    idx_clusts_use = (IC_stats.fracOfCluster_IC > minFracIC) & (IC_stats.EC_clustIds > 0);
                
    uClustIds_use = uClustIds(idx_clusts_use);
    clustSpkIdxs = clustSpkIdxs(idx_clusts_use);
    nClust = nnz(idx_clusts_use);
    clustNSpikes = cellfun(@length, clustSpkIdxs);            
    
    if nClust < 2
        fprintf('Only %d clusters.\n', nClust);
        return;
    end
    
    

    cellColors =  jet(nClust);    
                
    densityTypes = {'spike density', 'cluster density'};
    scaleTypes = {'linear', 'log'};
        
    nDensityBins = 64;
        
    densityType_prev = '';
    densityType = densityTypes{2};
    nSets = length(sortingFeatureSets); 
%     fet_idx = 1;
                
    nMult = 3;
    nFetsPerBasis = cellfun(@(x) size(x,2), features);        
    combos2D = arrayfun(@(n) getNMultPairs(n, nMult), nFetsPerBasis, 'un', 0);        
        
    [spkDensities, axesLims2D, axesTickLims2D] = deal(cell(1,nSets));
    for bi = 1:nSets
        [spkDensities{bi}, axesLims2D{bi}, axesTickLims2D{bi}] = ...
            getDensityHists(features{bi}, combos2D{bi}, clustSpkIdxs, nDensityBins, 0);
    end            

    cmap_spikes = [0 0 0; jet(64)];
    [cmap_clusters, cmapSclFactor_cell, cmapSclFactor_grp] = getStackedCmap(cellColors, 1);

    doColorBar = 1;
%     featureSet_idx0 = 1;
    featureSet_id = 1;
    
    densityType  = densityTypes{1};
    scaleType = scaleTypes{1};
    
    [h_2D_ax, h_2D_im, h_2D_xlab, h_2D_ylab] = deal(cell(1, nSets)); % deal( zeros(1,nSubplots) );
    
    for fet_idx = 1:nSets
        figId = 100+fet_idx;                
        figure(figId); clf; set(figId, 'Name', sortingFeatureSets{fet_idx});%, 'NumberTitle', 'off');
        %%
        if fet_idx == 1
            h_featureSetSelect = uicontrol('style', 'popupmenu', 'string', sortingFeatureSets, 'parent', figId, 'value', fet_idx, ...
                                  'units', 'normalized', 'position', [.05 .95, .20, .04], 'callback', @updateSpikeDensities);
            h_densityTypeSelect = uicontrol('style', 'popupmenu', 'string', densityTypes, 'parent', figId, 'value', 2, ...
                                  'units', 'normalized', 'position', [.30 .95, .20, .04], 'callback', @updateSpikeDensities);
            h_scaleSelect = uicontrol('style', 'popupmenu', 'string', scaleTypes, 'parent', figId, 'value', 1, ...
                                  'units', 'normalized', 'position', [.55 .95, .20, .04], 'callback', @updateSpikeDensities);
            h_showLabels_chk = uicontrol('style', 'checkbox', 'string', 'Labels', 'parent', figId, 'value', 1, ...
                                         'units', 'normalized', 'position', [.80 .95, .1, .04], 'callback', @updateSpikeDensityLabels);
            h_showTicks_chk = uicontrol('style', 'checkbox', 'string', 'Ticks', 'parent', figId, 'value', 0, ...
                                         'units', 'normalized', 'position', [.90 .95, .1, .04], 'callback', @updateSpikeDensityTicks);
        end
                          
        subM = 3; subN = 2;
        nSubplots = subM*subN;        
        
        densFigMargin = [0 0 .05 .05];
        for plot_i = 1:nSubplots
            h_2D_ax{fet_idx}(plot_i) = mySubplot(subM,subN,plot_i, [], 0, densFigMargin); 
            tks = axesTickLims2D{fet_idx}(plot_i,:);
            h_2D_im{fet_idx}(plot_i) = imagesc(tks(1:2), tks(3:4), zeros(nDensityBins, nDensityBins));
            
            [t1, t2] = deal(combos2D{fet_idx}(plot_i,1), combos2D{fet_idx}(plot_i,2));
            axis(h_2D_ax{fet_idx}(plot_i), axesLims2D{fet_idx}(plot_i,:) );
            axis(h_2D_ax{fet_idx}(plot_i), 'xy');
            h_2D_xlab{fet_idx}(plot_i) = xlabel(featureLabels{fet_idx}{t1}); 
            h_2D_ylab{fet_idx}(plot_i) = ylabel(featureLabels{fet_idx}{t2});
            set(h_2D_ax{fet_idx}(plot_i), 'xtick', [], 'ytick', []);
        end
        
        if doColorBar
            p6 = get(h_2D_ax{fet_idx}(end), 'position');
            dx = .001;
            dummy_pos = [p6(1:2) + [p6(3), 0], dx, p6(4)];
            h_dummy_ax(fet_idx) = axes('position', dummy_pos); 
            h_dummy_im(fet_idx) = imagesc(zeros(1,2)); 
            hColorBar(fet_idx) = colorbar;
            p_colbr = get(hColorBar(fet_idx), 'position'); 
            p_colbr(1) = p6(1)+p6(3)+dx;
            p_colbr(3) = p_colbr(3)*2;
            set(hColorBar(fet_idx), 'position', p_colbr);

            set(h_dummy_ax(fet_idx), 'xtick', [], 'ytick', [], 'position', p_colbr, 'visible', 'off');
        end

        curClims = [nan, nan];
        

        curCellGrouping = {1:nClust};
        curCellZoom = 1;
        
%         for plot_j = 1:nSubplots
%             ax = axesLims2D{fet_idx}(plot_j,:);     
%             ax_tck = axesTickLims2D{fet_idx}(plot_j,:);                 
%             
%             [cdata, curClims] = getSpikeDensityCData(spkDensities{fet_idx}, plot_j, 1, ...
%                 densityType, scaleType, curClims, cmapSclFactor_cell, cmapSclFactor_grp, clustNSpikes, curCellGrouping, curCellZoom);
% 
%             set(h_2D_im(plot_j), 'xdata', ax_tck(1:2), 'ydata', ax_tck(3:4), 'cdata', cdata' );                    
%             set(h_2D_ax(plot_j), 'xlim', ax(1:2), 'ylim', ax(3:4));
% 
%             [i1, i2] = deal(combos2D{fet_idx}(plot_j,1), combos2D{fet_idx}(plot_j,2));
%             set(h_2D_xlab(plot_j), 'string', featureLabels{fet_idx}{i1}); 
%             set(h_2D_ylab(plot_j), 'string', featureLabels{fet_idx}{i2});             
%         end
%         
        
    end    


    updateSpikeDensities;
    
    function updateSpikeDensities(~,~)                
        
%         featureSet_id = get(h_featureSetSelect, 'value');        
        densityType = densityTypes{get(h_densityTypeSelect, 'value')};
        scaleType = scaleTypes{get(h_scaleSelect, 'value')};        
        
        for fet_id = 1:nSets
        
            if ~strcmp(densityType_prev, densityType)
                switch densityType
                    case 'spike density',  colormap(h_2D_ax{fet_id}(1), cmap_spikes);
                    case 'cluster density', colormap(h_2D_ax{fet_id}(1), cmap_clusters);
                end                
            end


            showCells = true(1, nClust); % clustsAvailable & absUseClusters & useClusters(curNCells,:);
            nCells = nClust;
%             curClustIds = nan;            
            if nCells <= 6
                cellColors =  lines(nCells);
            else
                cellColors =  jet(nCells);
            end

    %         else % only show clusters in current cell
    %             showCells = false(1,nClust);
    %             curClustIds = curCellGrouping{curCellZoom};
    %             nCells = length(curClustIds);
    %             showCells(curClustIds) = clustsAvailable(curClustIds) & absUseClusters(curClustIds) & useClusters(curNCells, curClustIds);
    %                         
    %             cellColors = cat(1, clustColors{curClustIds});
    %         end
            curClims = [nan, nan];                       

            [cmap_clusters, cmapSclFactor_cell, cmapSclFactor_grp] = getStackedCmap(cellColors(1:nCells,:), 1);

            for plot_j = 1:nSubplots
                ax = axesLims2D{fet_id}(plot_j,:);     
                ax_tck = axesTickLims2D{fet_id}(plot_j,:);                 

                [cdata, curClims] = getSpikeDensityCData(spkDensities{fet_id}, plot_j, showCells, ...
                    densityType, scaleType, curClims, cmapSclFactor_cell, cmapSclFactor_grp, clustNSpikes, curCellGrouping, curCellZoom);

                set(h_2D_im{fet_id}(plot_j), 'xdata', ax_tck(1:2), 'ydata', ax_tck(3:4), 'cdata', cdata' );                    
                set(h_2D_ax{fet_id}(plot_j), 'xlim', ax(1:2), 'ylim', ax(3:4));

                [i1, i2] = deal(combos2D{fet_id}(plot_j,1), combos2D{fet_id}(plot_j,2));
                set(h_2D_xlab{fet_id}(plot_j), 'string', featureLabels{fet_id}{i1}); 
                set(h_2D_ylab{fet_id}(plot_j), 'string', featureLabels{fet_id}{i2});             
            end

            if any(isnan(curClims)),
                curClims = [0 1];
            end                

            switch densityType
                case 'spike density', 
                    set(h_dummy_im(fet_id), 'cdata', [curClims]);
                    set([h_2D_ax{fet_id} h_dummy_ax(fet_id)], 'climmode', 'auto');
                    set(hColorBar(fet_id), 'yTickMode', 'auto', 'yTickLabelMode', 'auto');
                case 'cluster density', 

                    colormap(h_2D_ax{fet_id}(1), cmap_clusters );
                    set(h_dummy_im(fet_id), 'cdata', [0, nCells]);
                    set(h_dummy_ax(fet_id), 'climMode', 'auto');
                    set(h_2D_ax{fet_id}, 'clim', [0 nCells]);
                    set(hColorBar(fet_id), 'yTick', []);
    %                         set(hColorBar, 'yTick', [0:nCells-1]+5, 'yTickLabel', cellfun(@num2str, num2cell(0:nCells-1), 'un', 0));
                    skp = ceil(nCells/10);%iff(nCells<10, nCells, round( nCells/10 );
                    cellIdIdx = 1:skp:nCells;
                    set(hColorBar(fet_id), 'yTickMode', 'auto', 'yTickLabelMode', 'auto');
                    yticks = uClustIds;
    %                 if ~isempty(curCellZoom)                     
    %                 else
    %                     yticks = cellIdIdx;
    %                 end

                    drawnow;
                    set(hColorBar(fet_id), 'yTick', [cellIdIdx]-.5, 'yTickLabel', arrayfun(@num2str, yticks, 'un', 0));                
    %                 set(hColorBar, 'yTick', [cellIdIdx]-.5, 'yTickLabel', arrayfun(@num2str, yticks, 'un', 0));                
                    3;
    %                 set(hColorBar, 'yTickMode', 'auto', 'ytickLabelMode', 'auto');
                    3;
            end          
        
        end        
        
        
        densityType_prev = densityType;        
    end
    

    function updateSpikeDensityLabels(~,~)
        
        showLabels = get(h_showLabels_chk, 'value');
        for fet_id = 1:nSets
            set([h_2D_xlab{fet_id}, h_2D_ylab{fet_id}], 'visible', tfToOnOff(showLabels));
        end
        
        
    end

    function updateSpikeDensityTicks(~,~)
        showTicks = get(h_showTicks_chk, 'value');
        
        for fet_id = 1:nSets
            if showTicks
                set(h_2D_ax{fet_id}, 'xtickmode', 'auto', 'ytickmode', 'auto');
            else
                set(h_2D_ax{fet_id}, 'xtick', [], 'ytick', []);
            end
        end
        
    end

    
    3;
end