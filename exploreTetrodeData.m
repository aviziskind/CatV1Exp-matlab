function exploreTetrodeData(Gid, autoPruningFlag, pruningMode)

    % autoPruningFlag: 0: no pruning
    %                  1: standard pruning of clusters, using automatic refr period.                      2: prune using ground truth IC refr period
    %                  11: standard pruning of IC cell (combined clusters), using automatic refr period. 12: prune IC cell (combined clusters) using ground truth IC refr period

    %%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    autoClassOpt.cutRemainTh = 0; % 0.66;
    autoClassOpt.assignToMU_mahalaSmallOrWide = 0;
    autoClassOpt.assignToMU_waveformMultPeaks = 0;
    autoClassOpt.meanWvfmSecondPeakRatio = 0.5;
    autoClassOpt.nFirstBins_mode = 1;  % maybe first *2* bins?
    [refrPeriodRange_ms, maxRefrPeriod_ms, minRefrPeriod_ms] = getGlobals(...
        'refrPeriod_ms_range', 'maxRefrPeriod_ms', 'minRefrPeriod_ms');    
    autoClassOpt.refrRange_ms = refrPeriodRange_ms;
    autoClassOpt.minRefrPeriod_ms = minRefrPeriod_ms;  
    autoClassOpt.maxRefrPeriod_ms = maxRefrPeriod_ms;  % had "5" here?
    autoClassOpt.refr_leeway_ms = .15;
    
    % 3. Spike waveforms (& Haar wavelet coefficients)    
%     showWaveformPCA = false;
%     showHad_and_PCA = true;
%     showWaveformGLF = false;

%     allFeaturesLabels    = {'Neg amps (raw)', 'Hadamard Transform', 'Neg amps (ccw)', 'Channel PCA', 'Waveform PCA', 'Channel GLF', 'Waveform GLF'};
%     allFeaturesAvailable = {'Neg_raw',        'Neg_had',            'Neg',            'PCAs',        'PCAc',         'GLFs',        'GLFc'        };
%     featuresToShow = {'Neg', 'Neg_had', 'PCA_sep.2', 'GLF_sep.2', 'GLF_cat.4'};
%     featuresToShow = {'Neg', 'PCAc16', 'GLFc8'};
    if (nargin < 2) || isempty(autoPruningFlag)
        autoPruningFlag = 0;
    end

    if nargin < 3 || isempty(pruningMode)
        pruningMode = 'clusterPruning';
    end
    
    doAutoPruning = (autoPruningFlag ~= 0);    
    useICrefrPeriod = curIdealPruning && (Gid > 6000);         
    
    switch pruningMode
        case 'clusterPruning',    clustGrouping = 'clusters';           cluster_data_idx = 1;
        case 'mergePruning',      clustGrouping = 'cells_onlyIC';       cluster_data_idx = 3;
        case 'pruneMergePruning', clustGrouping = 'cells_onlyPrunedIC'; cluster_data_idx = 4;        
        otherwise, error('Invalid pruningMode: %s', pruningMode)
    end
%     useCombinedICclusters = any(autoPruningFlag == [11, 12, 21, 22]);    
    
    useICrefrPeriod_useIndivClustRefrPeriod = 1;
    if useICrefrPeriod % || 1
        IC_stats = identifyCellfromIC(Gid, struct('clustGrouping', clustGrouping));
    end
    
    if useICrefrPeriod
        assert(curJustICclusters);
    end
    
    showAllFigures = ~doAutoPruning;
    pruningFeatures = curPruningFeatures('');
%     featuresToShow = {'Neg'};
    featuresToShow = {'Neg', 'Neg_raw', 'PCAcuw8', pruningFeatures}; %, 'GLFcuw32'};
    featuresToShow = uniqueInOrder(featuresToShow);

    forSlides = 0;
    
    if showAllFigures
        [~,~,~,curFet] = curSortingFeatures;

        if length(curFet) == 1
            featuresToShow(end+1) = curFet;
        else
            featuresToShow(end+1) = {curFet};
        end    
        featuresToShow = uniqueC(featuresToShow); % in case curFet was a duplicate.   
    end
    
    nFeatureSets = length(featuresToShow);
    nMaxPairsPerBasis = 20;    
    minClusterSize = [];

    basis_id0 = 1;    
    
    show2D   = 1;
    show2Dr  = 0;
    show3D   = 0;
    
    showISIs = 0;
    showEllipsoids0   = 0;
    showPairDists0    = 0;
    showISIs0         = 0;
    showClusterWaveforms0 = 0;
        whitenWvfms = 1;
    
    showTicks0 = false;        
    showLabels0 = true;
       
    nDensityBins = 64;

    subM_2D = 3;
    subN_2D = 2;
    
    [isi_allRanges_ms, isi_allNbins, refrPeriod_ms_range] = ...
        getGlobals('isi_allRanges_ms', 'isi_allNbins', 'refrPeriod_ms_range');

    refrPeriod_ms0 = 1; 

    nEllipsePoints2D = 50;
    nEllipsePoints3D = 20;    

        
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin < 1
        Gid = 1877;
    end
    if ~showAllFigures
        show2D = 0;
    end


    
    %%% A) Load features (spike amplitudes, PCA components, etc)    
    % 1. raw spike amplitudes
%     varName = ['group_' num2str(Gid)];
    siteInfo = siteDataFor(Gid);
    samplingRateHz = siteInfo.dataFileInfo.samplingRateHz;
    ticksPerMs = samplingRateHz/1000;

    matchDB = curMatchDB;
%     
    spikePropertiesFile = getFileName('properties', Gid);
    if matchDB && ~exist(spikePropertiesFile, 'file');

        S_db = dbGetSpikes(Gid, [], 'tick', 1);
        spikeTimes_tk = S_db(:,1);        
        negAmps = S_db(:,3:6);
        spikeCellIds = S_db(:,2);        
        
    else                
        S_spkProp = load(spikePropertiesFile);
        negAmps = S_spkProp.negAmps;
        spikeTimes_tk = S_spkProp.position;
        spikeMahalaHeights = S_spkProp.mahalanobisPeak;                       
        
        if doAutoPruning || show2Dr
            clusterGroupingToUse = clustGrouping; % iff(useCombinedICclusters, 'cells_onlyIC', 'clusters');
        else
            clusterGroupingToUse = 'clustersPruned';
        end
%         clusterGroupingToUse = 'clusters';
                
        spikeCellIds = getCellSorting(Gid, clusterGroupingToUse); 
                
        MDlims = lims(spikeMahalaHeights, 0);
        mdBinEdges = MDlims(1) : .5 : MDlims(2);
        mdBinCent = binEdge2cent(mdBinEdges);
        allMDs_bin = histcnt(spikeMahalaHeights, mdBinEdges);                        
    end            
    spikeTimes_ms = double(spikeTimes_tk) / ticksPerMs;
    nChannels = size(negAmps, 2);   

    [channelMeans, channelCov, channelWhitenMtx] = getChannelMeansAndCovariance(Gid);            
    invCovMtx = inv(channelCov);
    [Evect,Eval] = eig(channelCov);
    EvectMax = Evect(:, indmax(diag(Eval)));
    channelScaling = 1./sqrt(diag(channelCov));    
    
    
    % load features.
    features = cell(1,length(featuresToShow));
    featureLabels = cell(1,length(featuresToShow));
    isSingleFeature = cellfun(@ischar, featuresToShow);    
    allSingleFeatureNames = [featuresToShow(isSingleFeature), featuresToShow{~isSingleFeature}];
    
    featureSetOptions = cell(1,length(featuresToShow));    
    featureSetOptions(isSingleFeature) = featuresToShow(isSingleFeature);
    featureSetOptions(~isSingleFeature) = cellfun(@(cs) cellstr2csslist(cs, '_'), featuresToShow(~isSingleFeature), 'un', 0);
    
    [features_S, featureLabels_S] = getGroupFeatures(Gid, allSingleFeatureNames, [], []);        

    for fet_i = find(isSingleFeature)
        features{fet_i} = features_S.(featuresToShow{fet_i});        
        featureLabels{fet_i} = featureLabels_S.(featuresToShow{fet_i});        
    end
    for fet_i = find(~isSingleFeature)
        featureList_i = featuresToShow{fet_i};
        for j = 1:length(featureList_i)
            indivFets = cellfun(@(fld) features_S.(fld), featureList_i, 'un', 0);
            features{fet_i} = [indivFets{:}]; 

            indivLabels = cellfun(@(fld) featureLabels_S.(fld), featureList_i, 'un', 0);
            featureLabels{fet_i} = [indivLabels{:}];       
        end
    end                
        
    loadWaveForms    = showClusterWaveforms0 || showPairDists0 || show2Dr || doAutoPruning;    
    loadExtWaveforms = show2Dr && showAllFigures;
    
    if loadWaveForms
        [spkWaveforms, t_ms] = getSpikeWaveforms(Gid, [], [], 'raw');
%         t_ms = t_ms/ticksPerMs;
        spkWaveforms = single(spkWaveforms);
        [nT, nChannels, nSpk] = size(spkWaveforms); %#ok<NASGU> % note: the spike waveforms should already be mean-subtracted
    end
    if loadExtWaveforms        
        S_ext = load( getFileName('waveforms_ext', Gid));
        t_ms_ext = S_ext.t_ms_ext;
        spkWaveforms_ext = S_ext.spikeWaveforms_ext;
        clear('S_ext');
        t_tk_ext = round(t_ms_ext*ticksPerMs);        
    end

    
    spikeViewOptions = {'individual spikes', 'density'};
    spikeViewOption0 = spikeViewOptions{1};
    spikeViewDepVars = { {}, {'densityType', 'densityScale'} };

    densityTypeOptions = {'spike density', 'cluster density'};
    densityTypeOption0 = densityTypeOptions{1};
    
    densityScaleOptions = {'linear', 'log_10'};
    densityScaleOption0 = densityScaleOptions{1};

    % Remove cells/clusters with small spike-counts        
    spikeCellIds = single(spikeCellIds);
    [uSpkCellIds, cellSpikeCounts_tmp] = uniqueCount( spikeCellIds );
    
    if ~isempty(minClusterSize)
        cellsToUse = cellSpikeCounts_tmp >= minClusterSize | (uSpkCellIds == 0); %always include MU (will receive all discarded spikes)
    else
        cellsToUse = true(1, length(uSpkCellIds));
    end
    if showAllFigures
        disp([uSpkCellIds(:), cellSpikeCounts_tmp(:), cellsToUse(:)])
        fprintf('Total # spikes : %d\n', length(spikeCellIds));
    end

    % Compute list of spikes for each cell.    
    cellSpikeIdx = arrayfun(@(id) find(spikeCellIds == id), uSpkCellIds, 'un', 0);    
    for cell_Idx = find(~cellsToUse(:)')                
        spikeCellIds(cellSpikeIdx{cell_Idx}) = 0; % assign tiny units to multiunit.
    end
    cellIds = uSpkCellIds(cellsToUse);
    
%     cellIds(cellIds == 0) = [];
    
    cellSpikeIdx = arrayfun(@(id) find(spikeCellIds == id), cellIds, 'un', 0);
    nCells = length(cellIds);
    cellSpikeIdx_orig = cellSpikeIdx;   
    
    cellSpikeCounts = cellfun(@length, cellSpikeIdx);
    
%     disp('Current sorting');
%     disp([cellIds(:), cellfun(@length, cellSpikeIdx(:))]);
    
    cellIds_groups = num2cell(cellIds);
    
    siteHasMU = any(cellIds == 0);        
    
    % calculate  Means & Covariances of each cluster
    [clustMs, clustCovs, fetInEigBasis] = getClusterMeansCovs(cellSpikeIdx, features);
    
    if showClusterWaveforms0 || showPairDists0 || show2Dr || doAutoPruning
       % calculate mean waveforms for each cluster (in CCW space)                
        meanWaveforms_raw = getClusterMeanWaveforms(cellSpikeIdx, spkWaveforms);                        
    end
        
    % edits
%     edits = struct('refrSpikesDeleted', [], 'spikesCutOut', []);        
    
    
    cellColors =  jet(nCells-siteHasMU);
    styleStrs2D = arrayfun(@(cell_i) {'o', 'markersize', 1, 'color', cellColors(cell_i,:)}, 1:nCells-siteHasMU, 'un', 0 );
    styleStrs3D = arrayfun(@(cell_i) {'o', 'color', cellColors(cell_i,:)}, 1:nCells-siteHasMU, 'un', 0 );
    
    if siteHasMU
        cellColors = [1 1 1; cellColors];
        styleStrs2D = [{{'o', 'color', 'k', 'markersize',1}}, styleStrs2D];
        styleStrs3D = [{{'o', 'color', 'k', 'markersize',1}}, styleStrs3D];                        
    end
                        

    fig_2D = 102;
    fig_2Dr = 103;
    fig_3D = 106;    
    fig_pdist = 107;
    fig_isi = 108;
    fig_showCellsChkboxes = 109;
    
    slideFigs = [201:204];
    
%     title_str = sprintf('(Gid = %d. %d clusters. %d spikes)', Gid, length();
    Gid_str = sprintf('(Gid = %d)', Gid);
%     fig_wvfm = 107;

    isi_range_max = max(isi_allRanges_ms);
%     isi_nbins_max = max(isi_allNbins);

    rangesIgnoreMU = true;
    
    if rangesIgnoreMU
        idxForFetRange = spikeCellIds > 0;
    else
        idxForFetRange = true(size(spikeCellIds));
    end

    fet_mins = cellfun(@(x) min(x(idxForFetRange,:), [], 1), features, 'un', 0);
    fet_maxs = cellfun(@(x) max(x(idxForFetRange,:), [], 1), features, 'un', 0);            

    nBases = length(features);    
%     nMaxFeatures = max( nFetsPerBasis );

    if show2D || show3D
        nMult = 3;
        nFetsPerBasis = cellfun(@(x) size(x,2), features);
        
        combos2D = arrayfun(@(n) getNMultPairs(n, nMult), nFetsPerBasis, 'un', 0);        
        nPairsPerBasis = cellfun(@(c) size(c, 1), combos2D);
        
        nSubplotsPerPage = subM_2D*subN_2D;        
        nPagesPerBasis = ceil(nPairsPerBasis/nSubplotsPerPage);    
        nPagesMax = max( nPagesPerBasis );            
    end
    
    if show2D
    
        do2DLegend = false;
        do2DColorBar = true;        
        
        for bi = 1:nBases
            nPages = size(combos2D{bi},1)/nSubplotsPerPage;
            if abs(nPages - round(nPages)) > 1e-5
                combos2D{bi}(ceil(nPages)*nSubplotsPerPage,2) = 0;  % pad with zeros.
            end                
        end
        
        combos2D_0 = combos2D{basis_id0};
        h_2D_ax = zeros(1,nSubplotsPerPage);
        h_2D_xlab = zeros(1,nSubplotsPerPage);
        h_2D_ylab = zeros(1,nSubplotsPerPage);
        h_2D_spks = zeros(nSubplotsPerPage,nCells);
        h_2D_refr = zeros(nSubplotsPerPage,nCells);
        h_2D_im   = zeros(1,nSubplotsPerPage);
%         axesTicksShowing = false(1,nSubplotsPerPage);
%         axesLabelsShowing = false(1,nSubplotsPerPage);

         
        [spkDensities, axesLims2D, axesTickLims2D] = deal(cell(1,nBases));
        for bi = 1:nBases            
            [spkDensities{bi}, axesLims2D{bi}, axesTickLims2D{bi}] = ...
                getDensityHists(features{bi}, combos2D{bi}, cellSpikeIdx, nDensityBins, siteHasMU );
        end            
                        
        % Compute colormaps for the density plots
        cmap_spikes = [0 0 0; jet(64)];
        cellColors_dens = cellColors;
        if siteHasMU
            cellColors_dens(1,:) = ones(1,3)/1.3;  % make multiunit grey instead of black.
        end
        [cmap_clusters, cmapSclFactor_cell, cmapSclFactor_grp] = getStackedCmap(cellColors_dens, 1);
                
        
        figure(fig_2D); set(fig_2D, 'Name', ['Features in 2D ' Gid_str], 'NumberTitle', 'off'); clf;         
        for plot_i = 1:nSubplotsPerPage
            h_2D_ax(plot_i) = subplot(subM_2D,subN_2D,plot_i);  hold on; box on;
            [t1, t2] = deal(combos2D_0(plot_i,1), combos2D_0(plot_i,2));
            for cell_i = 1:nCells
                h_2D_spks(plot_i,cell_i) = plot(0,0, styleStrs2D{cell_i}{:});
                h_2D_refr(plot_i,cell_i) = plot(0,0, 'ko-');
            end        
            axis( axesLims2D{basis_id0}(plot_i,:) );
            h_2D_xlab(plot_i) = xlabel(featureLabels{basis_id0}{t1}); 
            h_2D_ylab(plot_i) = ylabel(featureLabels{basis_id0}{t2});
            h_2D_im(plot_i) = imagesc(zeros(nDensityBins, nDensityBins));                        
        end
                
        colormap(cmap_spikes);
        
        
        if do2DLegend 
            p = get(h_2D_ax(subN_2D), 'position');
            dx = .02;
            dummy_pos = [p(1:2) + [p(3), p(4)], dx, dx];
            h_dummy_ax = axes('position', dummy_pos); hold on;
            for cell_i = 1:nCells
                h_dummy_plot(cell_i) = plot(h_dummy_ax, 0, 0, styleStrs2D{cell_i}{:}); %#ok<AGROW>
            end                        
            fsize = iff(nCells<15, 10, 6);
            hLegend = legend(arrayfun(@(x) num2str(x), cellIds, 'un', 0), 'location', 'NW', 'fontsize', fsize); %#ok<NASGU>
            set(h_dummy_ax, 'xtick', [], 'ytick', []);
            set([h_dummy_ax, h_dummy_plot], 'visible', 'off')
        end                        
        if do2DColorBar
            p6 = get(h_2D_ax(end), 'position');
            dx = .02;
            dummy_pos = [p6(1:2) + [p6(3), 0], dx, p6(4)];
            h_dummy_ax = axes('position', dummy_pos); 
            h_dummy_im = imagesc(zeros(1,2)); 
            hColorBar = colorbar;
            p_colbr = get(hColorBar, 'position'); 
            p_colbr(1) = p6(1)+p6(3)+dx;
            p_colbr(3) = p_colbr(3)*2;
            set(hColorBar, 'position', p_colbr);
            
            set(h_dummy_ax, 'xtick', [], 'ytick', [], 'position', p_colbr, 'visible', 'off');
            3;
%             set([h_dummy_ax, , hColorBar], 'visible', 'off')
%             set(h_dummy_ax, );            
        end

        switch spikeViewOption0, 
            case 'individual spikes', set(h_2D_im, 'visible', 'off')
            case 'density',           set([h_2D_spks h_2D_refr], 'visible', 'off')
            otherwise, error('bad option');
        end
            
                
        if ~showTicks0
            set(h_2D_ax, 'xtick', [], 'ytick', []);
        end
        if ~showLabels0
            set([h_2D_xlab, h_2D_ylab], 'string', '');
        end

        
    end

    % Figure with checkboxes for which cells to show:    
    showCells = true(1,nCells);
    if showAllFigures && show2D
        h_2DCellChkbox =  zeros(2,nCells);
        figure(fig_showCellsChkboxes); clf; 
        set(fig_showCellsChkboxes, 'color', [.95 .95 .95], 'Name', ['Show Cells ' Gid_str], 'NumberTitle', 'off');
        positionShowCellCheckboxes;                   
    end
        
    function positionShowCellCheckboxes
        function P = chkBox_pos(i, M, N)
            [n, m] = ind2sub([N, M], i);
            L = (n-1)/N;
            B = (M-m+1)/(M+1);
            W = 1/N;
            H = 1/(2*(M+1));
            P = [L B W H];            
        end
        nSpksPerCell = cellfun(@length, cellSpikeIdx);
        nCells = length(cellIds);
        nMaxPerRow = 7;
        nShowCellsRows = ceil(nCells/nMaxPerRow);                           
        for c_j = 1:nCells                                    
            cellId = cellIds(c_j);
            pos = chkBox_pos(c_j, nShowCellsRows, nMaxPerRow);
            pos2 = pos; pos2(2) = pos2(2)- [1/(2*(nShowCellsRows+1))];
%             pos2 = pos; pos(2) = pos(2)- [1/(2*(nShowCellsRows+1))]
            if h_2DCellChkbox(1,c_j) == 0        
                h_2DCellChkbox(2,c_j) = uicontrol('style', 'text',  'parent', fig_showCellsChkboxes, ... 
                    'units', 'normalized', 'position', pos2, 'fontsize', 8, 'string', num2str(nSpksPerCell(c_j)));                

                h_2DCellChkbox(1,c_j) = uicontrol('style', 'checkbox',  'parent', fig_showCellsChkboxes, 'value', 1, ... 
                    'units', 'normalized', 'position', pos, 'fontsize', 10, 'fontweight', 'bold', ...
                    'string', num2str(cellId), 'ForegroundColor', cellColors_dens(c_j,:), 'userData', cellId, ...
                    'callback', @showCellChkbox_callBack);
                3;
            else                                
                set( h_2DCellChkbox(1,c_j), 'position', pos, 'string', cellList_str, ...
                    'ForegroundColor', cellColors_dens(cellId,:), 'userData', cellIds_groups{c_j}, 'visible', 'on' );                
            end            
        end
        set(h_2DCellChkbox(:,nCells+1:end), 'visible', 'off');
        
    end

    
    function showCellChkbox_callBack(~, ~)

        updateCellVisibilities;
    end

    function updateCellVisibilities
        showCells = logical( cell2mat( get(h_2DCellChkbox(1,:), 'value') )');
        vHandles.callUpdateFunction();
        
    end


    if show2Dr || doAutoPruning % showRefrPerSpikesFigure
%         clusterPruningFeatures = curPruningFeatures('');
        refrBasisId = find(strcmp(pruningFeatures, featureSetOptions));
%         nFeatureSets;
        switch clustGrouping
            case 'clusters',        clustPruning_fileType = 'clusterPrunings';
            case 'cells_onlyIC',    clustPruning_fileType = 'ICcellPrunings';
            case 'cells_onlyPrunedIC', clustPruning_fileType = 'ICprunedCellPrunings';
        end
        
        clusterPruningsFile = getFileName(clustPruning_fileType, Gid, 1);
%         refrBasisId = 1;
%         assert(show2D==true);        

        DELETE_1 = 1;
        DELETE_2 = 2;
        DELETE_BOTH = 3;
        CUT_SPIKES = -1;
        
        refr_cellIdx = 1;        
        refr_cellIdx_prev = nan;
        
        if showAllFigures
            
            figure(fig_2Dr), set(fig_2Dr, 'Name', sprintf('Current rotated cluster (Gid = %d)', Gid), 'NumberTitle', 'off'); clf;                  
            %%%% 1. Axes with scatter plot of refractory spikes.  %%%%        

    %         figure( fig_refrXY); set(fig_refrXY, 'Name', 'Refractory spikes : Amplitude vs Dist from center', 'NumberTitle', 'off' ); clf;
            if ~forSlides    
                h_refr_xy_ax = mySubplot(16,2, [2 9], 1);
            else
                figure(201); clf; h_refr_xy_ax = axes('units', 'pixels', 'position', [60 60 300 200]);
            end
                
    %         h_refr_ax = subplot(10,1,1:7);
            h_refr_xy = scatter(0,0, 'o'); hold on;
            h_refr_xy_del = scatter(0,0, 'x', 'm'); 
            h_refr_xy_cut = scatter(0,0, 's', 'r'); 
            set(gca, 'xscale', 'log'); 
            xlabel('Percent of spikes beyond larger of refractory spikes'); ylabel('Min(spike amplitudes)');
            h_refr_xy_cur = plot(0,0, 'ko', 'markersize', 11); 
            h_refr_xy_tit = title(' '); 

            figure(fig_2Dr);
            h_refr_ax_tmp = mySubplot(16,2, [9 10],1);
            p_tmp = get(h_refr_ax_tmp, 'outerPosition')-[0 .05 0 0];
            delete(h_refr_ax_tmp);
            h_refr_panel = uipanel('units', 'normalized', 'position', p_tmp);        

            % choose which cell
            uicontrol('style', 'text',  'units', 'normalized', 'position', [.02, .86, .07, .04], 'string', 'Cell #', 'horiz', 'cent', 'fontsize', 10);
            h_refrCellSelect = uicontrol('style', 'slider', 'units', 'normalized', 'position', [.1 .86, .3, .04], ...
                'min', 1, 'max', nCells, 'sliderstep', (1/(nCells-1))* [1, 2], 'value', 1, 'callback', @updateRefrCell);
            h_refrCellDisplay = uicontrol('style', 'text',  'units', 'normalized', 'position', [.41, .86, .05, .04], 'string', num2str(cellIds(1)), 'horiz', 'cent', 'fontsize', 10);

            doAutoDebug = 1;
            cb_upperLim = iff(doAutoDebug, 3, 1);
            cb_upperLim_ylim = iff(doAutoDebug, [.5, 3.5], [1.5, 2.5]);
            cb_hgt = iff(doAutoDebug, .08, .04);
            h_cellBlock_ax = axes('units', 'normalized', 'position', [.02 .91, .48, cb_hgt], 'xtick', [], 'ytick', []);
            cellBlock_cols = zeros(3, nCells);
            h_cellBlock_im = imagesc(cellBlock_cols);
            set(h_cellBlock_ax, 'xtick', [], 'ytick', [], 'clim', [-4 4]);
            pos_i = 1;
            h_cellBlock_select = rectangle('position', [.5+(pos_i-1), cb_upperLim_ylim(1), 1, cb_upperLim], 'linewidth', 3, 'edgecolor', 'w');
            set(h_cellBlock_ax, 'clim', [-4 4], 'ylim', cb_upperLim_ylim);
            set(h_cellBlock_im, 'ButtonDownFcn', @clickSelectCell);
            3;


            % refractory spike action buttons (delete, cut, load, save);
            nButtonsAcross = 6;
            nButtonsDown   = 2;
            but_pos = @(i,j) getNormPosition (nButtonsDown, nButtonsAcross, i,j, .02, .02);
    %         but_pos = @(i,j) [(4*(j+1)-3)/(4*nButtonsAcross),  [i/nButtonsDown  .12*(3-i)-0.1], 1/(1.2*nButtonsAcross), .09];        
            h_delSpk_button(1) = uicontrol('style', 'togglebutton',  'string', 'Del #1', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,1), 'userData', 1, 'callback', @deleteRefractorySpike);
            h_delSpk_button(2) = uicontrol('style', 'togglebutton',  'string', 'Del #2', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,2), 'userData', 2, 'callback', @deleteRefractorySpike);
            h_applyCut_button = uicontrol('style', 'togglebutton',  'string', 'Cut off', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,3), 'callback', @cutRefractorySpikes);
            h_discardClust_button = uicontrol('style', 'togglebutton',  'string', 'M-Unit', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,4), 'callback', @discardCluster);
            uicontrol('style', 'pushbutton',  'string', 'Auto Cut', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,5), 'callback', @autoCut_callback);
            uicontrol('style', 'pushbutton',  'string', 'Classify', 'parent', h_refr_panel, ...
                'units', 'normalized', 'position', but_pos(1,6), 'callback', @autoClust_callback);

            uicontrol('style', 'pushbutton',  'string', 'LOAD',  'units', 'normalized', 'position', but_pos(2,1), 'parent', h_refr_panel, 'callback', @loadEdits);
            uicontrol('style', 'pushbutton',  'string', 'SAVE',  'units', 'normalized', 'position', but_pos(2,2), 'parent', h_refr_panel, 'callback', @saveCurrentEdits);
            uicontrol('style', 'pushbutton',  'string', 'CLEAR', 'units', 'normalized', 'position', but_pos(2,3), 'parent', h_refr_panel, 'callback', @clearCurCellEdits);
            uicontrol('style', 'pushbutton',  'string', 'RESET(all)', 'units', 'normalized', 'position', but_pos(2,4), 'parent', h_refr_panel, 'callback', @clearAllCellEdits);
            uicontrol('style', 'pushbutton',  'string', 'AUTO', 'units', 'normalized', 'position', but_pos(2,5), 'parent', h_refr_panel, 'callback', @autoClassifyAllClusters);
            hnds_refr_allButtons = [h_delSpk_button, h_applyCut_button, h_discardClust_button];
            hnds_refr_cutButtons = [h_delSpk_button, h_applyCut_button];

            activeColor = [.6 .6 .6];      inactiveColor = [1 1 1];
            setButtonBkgColor = @(v_hnd) set(v_hnd, 'backgroundcolor', iff(get(v_hnd, 'value'), activeColor, inactiveColor) );        

    %         set(fig_refrXY,  'WindowKeyPressFcn',    @selectRefractorySpikeHelper); 
             %  'windowButtonDownFcn', {@selectPointsInFigure, @selectRefracSpike}, ...                       

            set(fig_2Dr,  'WindowKeyPressFcn',  @selectRefractorySpikeHelper); ...
    %             'WindowScrollWheelFcn', @scrollThroughRefracSpikes);


            %%%% 2. Axes with (current) cluster ISI. %%%%        
            if ~forSlides            
                h_refrCell_isi_ax =  mySubplot(16,2, [13 16],1);        
            else
                figure(202); clf; h_refrCell_isi_ax = axes('units', 'pixels', 'position', [60 60 320 180]); 
            end
    %         h_refrCell_isi_bar = bar([0 0], zeros(2,2), 1, 'edgecolor', 'none', 'facecolor', 'b');
            h_refrCell_isi_bar = bar([1 2], zeros(2,3), 1, 'stacked');
            set(h_refrCell_isi_bar, {'faceColor'}, {'b';'r';'g'});
    %         h_refrCell_isi_txt = text(0,0, ' ', 'fontsize', 8, 'hor', 'center', 'vert', 'top');
            h_refrCell_isi_vline = line(0,0, 'linestyle', ':');
            h_refrCell_isi_tit = title('');
            h_refrCell_isi_xlab = xlabel('');
            h_refrCell_isi_ylab = ylabel('');
            set(h_refrCell_isi_ax, 'xtick', [], 'ytick', []);

            isi_rangeOptions = arrayfun(@(rng_ms) sprintf('%d ms', rng_ms), isi_allRanges_ms, 'un', 0);
            isi_nbinOptions = arrayfun(@(nbin) sprintf('%d bins', nbin), isi_allNbins, 'un', 0);

            h_refrISI_range = uicontrol('style', 'popupmenu', 'string', isi_rangeOptions, 'parent', fig_2Dr, 'value', 3, ...
                                  'tag', 'range_ms', 'units', 'normalized', 'position', [.02 .25, .065, .04], 'callback', @updateRefrCellISI);
            h_refrISI_nbins = uicontrol('style', 'popupmenu', 'string', isi_nbinOptions, 'parent', fig_2Dr, 'value', 2, ...
                                  'tag', 'nbins', 'units', 'normalized', 'position', [.1 .25, .085, .04], 'callback', @updateRefrCellISI);

            n_refrPeriod_step = (1/length(refrPeriod_ms_range))*[1, 5];        
            h_refrISI_period_slider = uicontrol('style', 'slider', 'parent', fig_2Dr, 'value', refrPeriod_ms0, ...
                                     'min', refrPeriod_ms_range(1), 'max', refrPeriod_ms_range(end), 'sliderstep', n_refrPeriod_step, ...
                                     'tag', 'refrPeriod_slider', 'units', 'normalized', 'position', [.2 .25, .22, .04], 'callback', @updateRefrCell_RefrPeriod);
            h_refrISI_period_box = uicontrol('style', 'edit', 'parent', fig_2Dr, 'string', sprintf('%.2f ms', refrPeriod_ms0), ...
                                     'tag', 'refrPeriod_text', 'units', 'normalized', 'position', [.42 .25, .08, .04], 'callback', @updateRefrCell_RefrPeriod);
            h_refrISI_showAll_chk = uicontrol('style', 'checkbox', 'parent', fig_2Dr, 'string', 'All', 'value', 1, ...
                                     'tag', 'refrShowAllISIs', 'units', 'normalized', 'position', [.45 .2, .05, .05], 'callback', @updateRefrCellISI);

            h_refrISIcontrols = [h_refrISI_range, h_refrISI_nbins, h_refrISI_period_slider, h_refrISI_period_box, h_refrISI_showAll_chk];

            %%%% 3. Axes with waveforms of current (refractory) spikes %%%%        
    %         figure(fig_refrWvfms); set(fig_refrWvfms, 'Name', 'Waveforms of refractory-period spikes', 'NumberTitle', 'off'); clf; 
            if ~forSlides
                h_refrWvfm_ax = mySubplot(16,2, [1 8],2);        
            else
                figure(203); clf; h_refrWvfm_ax = axes('units', 'pixels', 'position', [60 60 300 200]); hold on;
            end
            h_refrSpkLines = line([0 0; 0 0], [1 1; 1 1]); hold on;
            h_refrMahala = plot(0,0, '.-', 0,0, 'k.-'); 
            h_wvfm_th = plot(0,0, 'k-');        
            h_refrWvfm = plot(0,0, '.-', 0,0, '.-', 0,0, '.-', 0,0, '.-'); hold on;            
            
            set(h_refrSpkLines(1), 'marker', 'o'); set(h_refrSpkLines(2), 'marker', 's')
            h_refrWvfm_tit = title('');
            h_refrWvfm_ylab = ylabel('Std deviations from baseline');
            h_refrWvfm_xlab = xlabel('ms');
            [refrCell_spkIdxs, curOffset_add, posSelectedSpikes, idxSelectedSpikes] = deal(0);        
            axis(h_refrWvfm_ax, 'tight');

            waveformCategoryOptions = {'Refractory spikes', 'All spikes'};
            waveformScaleOptions = {'CCW', 'Scaled', 'Raw'};        
            refrWaveformCategory = waveformCategoryOptions{1};
            refrWaveformCategory_prev = '';
            refrWaveformScale = waveformScaleOptions{2};
            h_mahala_txt(1) = text(0,0, ' ', 'vertical', 'top', 'horiz', 'center');
            h_mahala_txt(2) = text(0,0, ' ', 'vertical', 'top', 'horiz', 'center');

            r_h = .4; r_w = .4;
            wvfm_ax_pos = get(h_refrWvfm_ax, 'position');        
            meanWvfm_ax_pos = [wvfm_ax_pos(1)+ wvfm_ax_pos(3)*(1-r_w), wvfm_ax_pos(2), wvfm_ax_pos(3)*r_w, wvfm_ax_pos(4)*r_h];
            h_refrMeanWvfm_ax = axes('position', meanWvfm_ax_pos);
            h_refrMeanWvfm = plot(1:2, zeros(4,2));
            axis(h_refrMeanWvfm_ax, 'tight');
            set(h_refrMeanWvfm_ax, 'xtick', [], 'yAxisLocation', 'right', 'ylimmode', 'auto');

            meanWaveformCurScale = '';
            meanWaveformCurCellIdx = nan;
            curMeanWvfm_raw = meanWaveforms_raw(:,:,1)';
            
            % controls for 'all spikes'/'refractory spikes'  and 'CCW/Scaled'Raw' 
            h_refrWvfm_category = uicontrol('style', 'popupmenu', 'string', waveformCategoryOptions, 'parent', fig_2Dr, 'value', 1, ...
                      'tag', 'wvfm_select', 'units', 'normalized', 'position', [.56 .95, .15, .04], 'callback', @updateListOfCellSpikes);
            uicontrol('style', 'popupmenu', 'string', waveformScaleOptions, 'parent', fig_2Dr, 'value', 2, ...
                      'tag', 'wvfm_scale', 'units', 'normalized', 'position', [.73 .95, .15, .04], 'callback', @updateWaveformScale);        
            h_refrWvfm_ylock = uicontrol('style', 'checkbox', 'string', 'Lock Y', 'parent', fig_2Dr, 'value', 0, ...
                                'tag', 'wvfm_lockY', 'units', 'normalized', 'position', [.90 .95, .08, .04]);        


            %%%% 4. Figure with (current) rotated/projected cluster. %%%%        
            if ~forSlides
                h_2Dr_clust_ax = mySubplot(16,2, [9 15], 2);
            else
                figure(204); clf; h_2Dr_clust_ax = axes('units', 'pixels', 'position', [60 60 220 220]); 
                axis equal;
            end
            
            h_2Dr_spks = plot(0,0, 'bo', 'markersize', 3, 'userData', 'skipSelect'); hold on;
            h_2Dr_refr = plot(0,0, 'ko-');
            h_2Dr_refrSel(1) = plot(0,0, 'ko-', 'lineWidth', 2, 'markerfacecolor', 'k', 'userData', 'skipSelect', 'visible', 'off');
            h_2Dr_refrSel(2) = plot(0,0, 'ks-', 'lineWidth', 2, 'markerfacecolor', 'k', 'userData', 'skipSelect', 'visible', 'off');
            h_2Dr_spksToDel = plot(0,0, 'o', 'color', [.6 .6 .6], 'markersize', 3, 'userData', 'skipSelect'); hold on;
            h_2Dr_spksDel = plot(0,0, 'o', 'color', [.1 .1 .1], 'markersize', 3, 'userData', 'skipSelect'); hold on;
            h_2Dr_tit = title('');
            h_2Dr_xlab = xlabel('');
    %         h_2Dr_ylab = ylabel('');
            h_autoRot_chkbox   = uicontrol('style', 'checkbox', 'units', 'normalized', 'position', [.55 .00, .08, .03], 'parent', fig_2Dr, 'string', 'Rotate');                
            h_showMD_chkbox    = uicontrol('style', 'checkbox', 'units', 'normalized', 'position', [.70 .00, .10, .03], 'parent', fig_2Dr, 'string', 'Show MD', 'callback', @updateMD);                
            h_showAllMD_chkbox = uicontrol('style', 'checkbox', 'units', 'normalized', 'position', [.80 .00, .10, .03], 'parent', fig_2Dr, 'string', 'Show All', 'callback', @updateMD);
            h_logMD_chkbox     = uicontrol('style', 'checkbox', 'units', 'normalized', 'position', [.90 .00, .10, .03], 'parent', fig_2Dr, 'string', 'Log count', 'callback', @updateMD);


            showMD_prev = nan;
            showMD = 0;
            figure(fig_2Dr);
            h_md_ax = mySubplot(16,2, [9 16],2);
            h_md_bar = bar([1 2], zeros(2,3), 1, 'stacked', 'edgecolor', 'none'); hold on;
            h_md_sm = plot(0,0, 'k.-');
            h_md_tit = title(' ');
            h_md_xlab = xlabel(' ');
            h_md_ylab = ylabel(' ');
            set(h_md_bar, {'faceColor'}, {'g';'r';'b'});
            set(h_md_ax, 'xlim', MDlims);                
            set([h_md_ax h_md_bar h_md_sm], 'visible', 'off');        


            % Button to show/hide all CCGs
            h_showCCGs_button = uicontrol('style', 'togglebutton',  'string', 'CCG', 'parent', fig_2Dr, ...
                'units', 'normalized', 'position', [.47, .86, .05, .04], 'callback', @showHideRefrCellCCGs);

            ccg_show = false;
            set([h_refr_panel h_refrISIcontrols], 'visible', 'off');

            h_ccg_panel = uipanel('units', 'normalized', 'position', [.01 .01, .75 .87], 'parent', fig_2Dr);
            ccgTypeOptions = {'CCG', 'ISI'};
            h_ccg_type_popup = uicontrol('style', 'popupmenu', 'string', ccgTypeOptions, 'parent', h_ccg_panel, 'value', 1, ...
                                  'tag', 'range_ms', 'units', 'normalized', 'position', [.02 .95, .12, .03], 'callback', @updateRefrCellCCG);
            h_refrCCG_range = uicontrol('style', 'popupmenu', 'string', isi_rangeOptions, 'parent', h_ccg_panel, 'value', 3, ...
                                  'tag', 'range_ms', 'units', 'normalized', 'position', [.2 .95, .15, .04], 'callback', @updateRefrCellCCG);
            h_refrCCG_nbins = uicontrol('style', 'popupmenu', 'string', isi_nbinOptions, 'parent', h_ccg_panel, 'value', 2, ...
                                  'tag', 'nbins', 'units', 'normalized', 'position', [.4 .95, .15, .04], 'callback', @updateRefrCellCCG);

            ccgFigMargin = [.01 .01 .01 .1];            
            nCCGs = nCells-1;
            ccg_N = floor(sqrt(nCCGs)); 
            ccg_M = ceil((nCCGs)/ccg_N); 
            [h_ccg_hist_ax, h_ccg_wvfm_ax, h_ccg_hist_bar, h_ccg_txt] = deal( zeros(nCCGs, 1) );
            h_ccg_wvfms = zeros(nCCGs, 4);        
            for cj = 1:nCells-1
                h_ccg_hist_ax(cj) = mySubplot(ccg_M, ccg_N, cj, [], 0, ccgFigMargin, 'parent', h_ccg_panel);
                h_ccg_hist_bar(cj) = bar(0,0, 1);            
                h_ccg_wvfm_ax(cj) = mySubplot(ccg_M, ccg_N, cj, [], 0, ccgFigMargin, 'parent', h_ccg_panel, 'color', 'none');
                linkprop([h_ccg_hist_ax(cj), h_ccg_wvfm_ax(cj)], 'position');
                h_ccg_wvfms(cj,:) = plot(h_ccg_wvfm_ax(cj), [1:2], zeros(2,4));            
                set(h_ccg_wvfm_ax(cj), 'color', 'none', 'xtick', [], 'ytick', [], 'xlim', [1 nT]);
                set(h_ccg_hist_ax(cj), 'xtick', []);
                h_ccg_txt(cj) = text(0,1, 'A', 'parent', h_ccg_hist_ax(cj), 'horiz', 'center', 'vert', 'top');            
            end
            3;

            idx_cellMarked = [];


            refr_cellIdx = 1;
            rot2D_allFeatures = fetInEigBasis{refrBasisId, refr_cellIdx};
            rot2D_projFeatures = rot2D_allFeatures;
            refrSpkInRangeSelect_idx = zeros(1,0);
            refrSpkSelect_idx = zeros(1,0);        

            if ~showTicks0
    %             set(h_2Dr_ax, 'xtick', [], 'ytick', []);
            end
            if ~showLabels0
    %             set([h_2Dr_xlab, h_2Dr_ylab], 'string', '');
            end        
            3;
            set([h_refr_panel h_refrISIcontrols], 'visible', 'on');
            set(h_ccg_panel, 'visible', 'off');

            if forSlides
                set([201:204],  'WindowKeyPressFcn',  @selectRefractorySpikeHelper); ...
            end

            
        end
    end
    
    
    function updateRefrCellCCG(~, ~)
        if ~get(h_showCCGs_button, 'value')  % only update if button is pressed
            return;
        end

        ccg_or_isi = ccgTypeOptions{ get(h_ccg_type_popup, 'value') };
    
        range_idx = get(h_refrCCG_range, 'value');          
        nbins_idx = get(h_refrCCG_nbins, 'value');        
        range_ms = isi_allRanges_ms(range_idx);
        nbins = isi_allNbins(nbins_idx);
%         mtxIdxs = @(j) sub2indV([nCells nCells], [ones(1, j-1)*j,   j+1:nCells; 1:j-1,     ones(1, nCells-j)*j]);        
%         mtxIdxs = @(j) sub2indV([nCells nCells], [ones(1, nCells-1)*j; 1:j-1, j+1:nCells]);

        switch ccg_or_isi
            case 'CCG',
                binEdges = linspace(-range_ms, range_ms, nbins+1);
                binCents = binEdge2cent(binEdges);                    
                masterBins = allCCGs_masterBin;
                masterBinIdxs = ccg_masterBinIdxs{range_idx, nbins_idx};                
            case 'ISI',                    
                binEdges = linspace(0, range_ms, nbins+1);
                binCents = binEdge2cent(binEdges);
                masterBins = allISIs_masterBin;
                masterBinIdxs = isi_masterBinIdxs{range_idx, nbins_idx};                
        end                        

        switch refrWaveformScale
            case 'CCW',    wvfmScaling = @(wvfm) channelWhitenMtx * wvfm;
            case 'Scaled', wvfmScaling = @(wvfm) bsxfun(@times, wvfm, channelScaling);
            case 'Raw',    wvfmScaling = @(wvfm) wvfm;                
        end
        
        cellIdxs = setdiff(1:nCells, [refr_cellIdx, find(refrClustAction == 1)]);
        nCCGs = length(cellIdxs);
        
        ccg_N = floor(sqrt(nCCGs)); 
        ccg_M = ceil(nCCGs/ccg_N);                 
        
        set([h_ccg_hist_ax(1:nCCGs),     h_ccg_wvfm_ax(1:nCCGs),     h_ccg_wvfms(1:nCCGs,:)],     'visible', 'on');
        set([h_ccg_hist_ax(nCCGs+1:end), h_ccg_wvfm_ax(nCCGs+1:end), h_ccg_wvfms(nCCGs+1:end,:)], 'visible', 'off');                               
        set(h_ccg_hist_ax(1:nCCGs), 'xlim', [binEdges(1), binEdges(end)], 'xtickmode', 'auto');
        
        3;
%         idxs = mtxIdxs(refr_cellIdx);
        allBinVals = masterBin2Bin(masterBins(refr_cellIdx,[1:refr_cellIdx-1, refr_cellIdx+1:nCells]), masterBinIdxs);
        for cidx = 1:length(cellIdxs)            
            h_ccg_hist_ax(cj) = mySubplot(ccg_M, ccg_N, cidx, [], h_ccg_hist_ax(cj), ccgFigMargin);
            ymax = max(allBinVals{cidx})*1.1+.1;
            set(h_ccg_hist_bar(cidx), 'xdata', binCents, 'ydata', allBinVals{cidx});
            set(h_ccg_hist_ax(cidx), 'ylim', [0 ymax]);

            wvfm_i = wvfmScaling(meanWaveforms_raw(:,:,cellIdxs(cidx))');
            for ch_i = 1:4
                set(h_ccg_wvfms(cidx, ch_i), 'xdata', 1:nT, 'ydata', wvfm_i(ch_i,:), 'linewidth', 2)
            end       
            set(h_ccg_txt(cidx), 'string', sprintf(' C(%d,%d)', cellIds([refr_cellIdx, cellIdxs(cidx)])), 'position', [-range_ms ymax], 'horiz', 'left')
            set(h_ccg_wvfm_ax(cidx), 'xlim', [1 nT]);
        end        
        
    end
    
    
    function showHideRefrCellCCGs(~,~)
    
        ccg_newShow = get(h_showCCGs_button, 'value');
        valChanged = ccg_show ~= ccg_newShow;
        ccg_show = ccg_newShow;
        if valChanged
            if ccg_show
                set(h_ccg_panel, 'visible', 'on');
                set([h_refr_panel h_refrISIcontrols h_autoRot_chkbox], 'visible', 'off');
                
                updateRefrCellCCG;
            else
                set(h_ccg_panel, 'visible', 'off');
                set([h_refr_panel h_refrISIcontrols h_autoRot_chkbox], 'visible', 'on');
            end            
        end        
        
    end


    if show2Dr || showISIs || doAutoPruning
        ClustData = getClusterData(Gid, cluster_data_idx);
    end
        
    if showISIs0 || show2Dr || doAutoPruning
                
        isi_range0 = 4;
        isi_nbin0 = 10;                                               
        nISIcolors = 64;
        if ~doAutoPruning
            isiColors = jet(nISIcolors);
        end
        refrPeriod_ms_max = max(refrPeriod_ms_range);
                
%     S.isiData = struct('isis_ms', allISIs_ms, 'isis_mb', allISIs_masterBin, 'mBinEdges', isi_masterBinEdges, 'mBinIdxs', isi_masterBinIdxs);
%     S.ccgData = struct('ccgs_ms', allCCGs_ms, 'ccgs_mb', allCCGs_masterBin, 'mBinEdges', ccg_masterBinEdges, 'mBinIdxs', ccg_masterBinIdxs);
%     S.autoRefrData = struct('allCell_refrSpkIdxs', allCell_refrSpkIdxs);
                        
        S_isi = ClustData.isiData;
        [allISIs_masterBin, isi_masterBinIdxs, isi_masterBinEdges ] = deal(...
            S_isi.isis_mb, S_isi.mbIdxs, S_isi.mbEdges);

        haveCCGData = isfield(ClustData, 'ccgData');
        if haveCCGData
            S_ccg = ClustData.ccgData;
            [allCCGs_masterBin, ccg_masterBinIdxs] = deal(...
                S_ccg.ccgs_mb, S_ccg.mbIdxs);
        end
        
        isi_range0_idx = find(isi_range0 == isi_allRanges_ms, 1); assert(~isempty(isi_range0_idx))
        isi_nbin0      = find(isi_nbin0  == isi_allNbins, 1);     assert(~isempty(isi_nbin0))
                
        if showAllFigures
            allISIs_bin = masterBin2Bin(allISIs_masterBin, isi_masterBinIdxs{isi_range0_idx, isi_nbin0});
            if haveCCGData
                allCCGs_bin = masterBin2Bin(allCCGs_masterBin, ccg_masterBinIdxs{isi_range0_idx, isi_nbin0});
            end
        end
        
        multiunitISI_ms = getISIs(spikeTimes_ms, true, isi_range_max);
        [refracViol_n, refracViol_pct, isiColorIdxs] = countRefractoryPeriodViolations(allISIs_masterBin, isi_masterBinEdges, refrPeriod_ms0, multiunitISI_ms, nISIcolors, cellSpikeCounts);
        
        S_refr = ClustData.autoRefrData;        
        if ~isequal(S_refr.clusterPruningFeatures, pruningFeatures)
            error('Refractory calculations were made with different features. Please recalculate');
        end
        [allCell_refrSpkIdxs, refrSpikesInfo, refrSpikeMarkedIdxs] = deal(...
            S_refr.allCell_autoRefrSpkIdxs, S_refr.refrSpikesInfo, S_refr.refrSpikeMarkedIdxs);
        allCellsFracsRemoved = cell(1, nCells);
        
        nRefrSpikes = cellfun(@length, allCell_refrSpkIdxs);
        
        % Variables to keep track of all refractory spikes & which ones we've deleted  
        curCellEdits = struct;               
        clustSpikesRemoved_idx   = cell(1, nCells);        
        refrSpkAction       = arrayfun(@(n) zeros(1, n), nRefrSpikes, 'un', 0);
        refrClustAction     = zeros(1,nCells);        
%         nRefrSpikesRemaining = nRefrSpikes;                                
               
        cellsRefrPeriod_ms  = ones(1,nCells)*refrPeriod_ms0;                      
        refrSpksInRange = cell(1,nCells);
        refrSpksInRange_list = cell(1,nCells);                        
        for cell_i = 1:nCells  % (xx) skip multi-unit, so start at #2        
            refrSpksInRange{cell_i} = refrSpikesInfo(cell_i).isi_ms <= cellsRefrPeriod_ms(cell_i);      
            refrSpksInRange_list{cell_i} = find(refrSpksInRange{cell_i}); 
        end
        nRefrSpikesInitially = cellfun(@length, refrSpksInRange_list)';
        nRefrSpikesRemaining = nRefrSpikesInitially;
        
        cellBlock_cols(1,:) = log(1 + nRefrSpikesInitially(:)');
        cellBlock_cols(1, uSpkCellIds == 0) = nan;
        cellBlock_cols(2,:) = cellBlock_cols(1,:);
        if showAllFigures
            autoClassifyAllClusters; % for cellBlock_cols(3)
        end
        
        if showAllFigures
            set(h_cellBlock_im, 'CData', cellBlock_cols);        
        
            if showISIs0 
                figure(fig_isi); set(fig_isi, 'Name', 'Interspike intervals', 'NumberTitle', 'off'); clf;

                [h_isi_ax, h_isi_bar, h_isi_txt, h_isi_tit] = deal( zeros(nCells) );
                h_isi_refrLine = zeros(nCells, 1);

                for cell_i = 1:nCells
                    for cell_j = 1:cell_i
                        h_isi_ax(cell_i,cell_j) = mySubplot(nCells,nCells, nCells-cell_i+1, cell_j); 
                        h_isi_bar(cell_i,cell_j) = bar(0, 0, 1, 'edgecolor', 'none', 'facecolor', isiColors(isiColorIdxs(cell_i, cell_j), :));
                        h_isi_txt(cell_i, cell_j) = text(0,0, ' ', 'fontsize', 8, 'hor', 'center', 'vert', 'top');
                        h_isi_refrLine(cell_i) = line(0,0, 'linestyle', ':');
                        h_isi_tit(cell_i, cell_j) = title('');
                        set(h_isi_ax(cell_i,cell_j), 'xtick', [], 'ytick', []);
                    end
                end           


                posLL = get(h_isi_ax(1,1), 'position');
                posUR = get(h_isi_ax(nCells,nCells), 'position');
                tickaxpos = [posLL(1:2), posUR(1:2)+posUR(3:4)-posLL(1:2)];                
                h_isi_tick_ax = axes('position', tickaxpos, 'color', 'none', 'xAxisLocation', 'top', 'tickDir', 'out', ...
                    'xlim', [.5 nCells+.5], 'ylim', [.5 nCells+.5], 'xtick', 1:nCells, 'ytick', 1:nCells, ...
                    'xticklabel', arrayfun(@num2str, cellIds, 'un', 0), 'yticklabel', arrayfun(@num2str, cellIds, 'un', 0));
                3;
            end
        
            if show2Dr                
                refresh(fig_2Dr);
            end
        end
        updateRefractorySpikes_curClust(1);
        
    end
    
    
    function scrollThroughRefracSpikes_helper(~,evnt)         

        scrollThroughRefracSpikes( evnt.VerticalScrollCount );

    end

    function scrollThroughCells(offset)
        cur_refr_cellIdx = refr_cellIdx;
        refr_cellIdx = bound(refr_cellIdx+offset, 1, nCells);
        if cur_refr_cellIdx ~= refr_cellIdx
            updateRefrCell;
        end
        
        
    end


    function scrollThroughRefracSpikes(scrollCount)         
        % refrSpkInRngSelect_idx == the index of the selected spike out of *all* refractory spikes of the cell
        % refrSpkSelect_idx      == the index of the selected spike out of the refractory spikes of the cell that are in range.  
                
        if ~isempty(refrCell_spkIdxs)
            if isempty(refrSpkInRangeSelect_idx)
                refrSpkInRangeSelect_idx = 1;
            end
%             switch refrWaveformCategory
%                 case 'Refractory spikes',
%                     all_idxs = refrSpksInRange_list{refr_cellIdx};
% %                     refrSpkSelect_idx = indmin( abs(refrSpkInRangeSelect_idx - all_idxs) );
%                 case 'All spikes';
%                     all_idxs = 1:cellSpikeCounts(refr_cellIdx);
% %                     refrSpkSelect_idx = refrSpkInRangeSelect_idx;
%             end
            
            refrSpkInRangeSelect_idx = refrSpkInRangeSelect_idx + scrollCount;
            updateRefractorySpikeIdx;
%             refrSpkInRangeSelect_idx  = bound( refrSpkInRangeSelect_idx + scrollCount, 1, length(all_idxs));
%             refrSpkSelect_idx  = all_idxs(refrSpkInRangeSelect_idx);
            
            updateSelectedRefractorySpike;
%             set(h_refr_xy_tit, 'string', '');
            set(h_refr_xy_tit, 'string', sprintf('Idx = %d. Id: %d.', refrSpkSelect_idx, refrSpkInRangeSelect_idx));
            
        end
    end    


    function selectRefractorySpikeHelper(src, evnt)
        3;
        if isfield(evnt, 'Key') 
            3;
            switch evnt.Key
                case {'1', '2'}, 
                    src = h_delSpk_button( str2double(evnt.Key) );
                    toggleValue(src);          % button up/down
                    deleteRefractorySpike(src);
                case '3',         
                    toggleValue(h_applyCut_button);
                    cutRefractorySpikes;
                case '9',                             
                    toggleValue(h_discardClust_button); 
                    discardCluster;
                case 'uparrow',    scrollThroughRefracSpikes(1);
                case 'downarrow',  scrollThroughRefracSpikes(-1);
                case 'leftarrow',  scrollThroughCells(-1);
                case 'rightarrow', scrollThroughCells(1);
                case 'pageup',     scrollThroughRefracSpikes(10);
                case 'pagedown',   scrollThroughRefracSpikes(-10);
                case 'home',       scrollThroughRefracSpikes(-inf);
                case 'end',        scrollThroughRefracSpikes(inf);
                case 'x', 
                    if gca == h_2Dr_clust_ax
                        selectPointsInFigure(src, evnt, @selectRefracSpike_clust);
                    elseif gca == h_refr_xy_ax
                        selectPointsInFigure(src, evnt, @selectRefracSpike_xy);
                    end                    
                    
%                 case 'backspace',  clearCurCellEdits
            end
        end
        3;
    end

    function selectRefracSpike_xy(glob_id, group_id, local_id, pt_x, pt_y)  %#ok<*INUSL>
%         refrSpkInRangeSelect_idx = ceil(glob_id/3);
            3;
%         inRange = refrSpksInRange{refr_cellIdx};
        
        switch group_id
            case 1, refrSpkInRangeSelect_idx = curCellEdits.idx_0(local_id);
            case 2, refrSpkInRangeSelect_idx = curCellEdits.idx_del(local_id);
            case 3, refrSpkInRangeSelect_idx = curCellEdits.idx_cut(local_id);
        end
%         refrSpkSelect_id_inRange = x(local_id);        
%         spikesInRange = find(refrSpksInRange{refr_cellIdx});
%         refrSpkInRangeSelect_idx = spikesInRange(refrSpkSelect_id_inRange);
        
%         idx_del1 = find(refrSpkAction{refr_cellIdx} == 1 | refrSpkAction{refr_cellIdx} == 3);
%         idx_del2 = find(refrSpkAction{refr_cellIdx} == 2 | refrSpkAction{refr_cellIdx} == 3);
%         idx_cut  = find(refrSpkAction{refr_cellIdx} == -1);
% 
%         h_refr_xy = plot(0,0, 'o'); hold on;
%         h_refr_xy_del = plot(0,0, 'mx'); 
%         h_refr_xy_cut = plot(0,0, 'rs'); 
        updateSelectedRefractorySpike;        

    end

    function selectRefracSpike_clust(glob_id, group_id, local_id, x, y) 
        
        switch refrWaveformCategory
            case 'Refractory spikes',                
                refrSpkSelect_id_inRange = ceil(glob_id/3);
                refrSpkInRangeSelect_idx = refrSpkSelect_id_inRange;                
%                 refrSpkInRangeSelect_idx = refrSpksInRange_list{refr_cellIdx}(refrSpkSelect_id_inRange);                
            case 'All spikes'
                refrSpkInRangeSelect_idx = glob_id;
        end
        
        updateSelectedRefractorySpike;
    end


    function updateRefractorySpikes_allClusters
        cur_refr_cellIdx = refr_cellIdx;
        for cl_i = 1:nCells
            refr_cellIdx = cl_i;
            updateRefractorySpikes_curClust(1);
        end
        
        refr_cellIdx = cur_refr_cellIdx;
        updateRefractorySpikes_curClust;
    end

%     function updateRefractorySpikes_allClusters2
%         for ci = 1:nCells
%             updateRefractorySpikesForClust(ci);
%         end        
%     end


    function clustEdits = updateRefractorySpikesForClust(clustIdx)
    
        % re-apply cuts & deletes. then update view.
        inRefrRange = refrSpikesInfo(clustIdx).isi_ms <= cellsRefrPeriod_ms(clustIdx);
        if ~isempty(inRefrRange)
            inRefrRange = inRefrRange';
        end
        refrSpksInRange{clustIdx} = inRefrRange;                
%         inRefrRange_tmp = true(size(inRefrRange_tmp));
        refrSpksInRange_list{clustIdx} = find(inRefrRange);
        clustEdits.idx_0    = find( inRefrRange & (refrSpkAction{clustIdx} == 0 ));
        clustEdits.idx_del  = find( inRefrRange & (refrSpkAction{clustIdx} >= 1 ));
        clustEdits.idx_del1 = find( inRefrRange & (refrSpkAction{clustIdx} == 1 | refrSpkAction{clustIdx} == 3) );
        clustEdits.idx_del2 = find( inRefrRange & (refrSpkAction{clustIdx} == 2 | refrSpkAction{clustIdx} == 3) );
        clustEdits.idx_cut  = find( inRefrRange & (refrSpkAction{clustIdx} == -1 ));                
        if ~isempty(refrSpikesInfo(clustIdx).dep)
            clustEdits.idx_cutOut = find( any( refrSpikesInfo(clustIdx).dep(:,clustEdits.idx_cut), 2) );                
        else
            clustEdits.idx_cutOut = [];                
        end
            
        
        %%%% Store the total number of spikes removed - both from indiv. deletes, and from cuts :
        % individual spike deletes
        indivSpksDelete = unique([allCell_refrSpkIdxs{clustIdx}(clustEdits.idx_del1); ...
                                  allCell_refrSpkIdxs{clustIdx}(clustEdits.idx_del2)+1]); %#ok<*FNDSB>
        % entire cuts
        spikeCuts = unique(cat(1, refrSpikeMarkedIdxs{clustIdx}{clustEdits.idx_cut}));        
        
        % store in the 'clustSpikesRemoved_idx' cell array.
        clustSpikesRemoved_idx{clustIdx} = unique([indivSpksDelete(:); spikeCuts(:)]);
        
        curRemainingRefrSpikes = refrSpksInRange{clustIdx};
        curRemainingRefrSpikes(clustEdits.idx_del) = 0;
        curRemainingRefrSpikes(clustEdits.idx_cut) = 0;
        curRemainingRefrSpikes(clustEdits.idx_cutOut) = 0;
        nRefrSpikesRemaining(clustIdx) = nnz(curRemainingRefrSpikes);                
        
    end

    function updateRefractorySpikes_curClust(suppressShowFlag)
                
        suppressShow = nargin > 0  && ~isempty(suppressShowFlag);
        curCellEdits = updateRefractorySpikesForClust(refr_cellIdx);
                
        updateListOfCellSpikes;    
        %%%% Update the refractory spike index (make sure is pointing to a spike that is in range  
        
        if ~suppressShow && showAllFigures
            updateRefractorySpikeIdx;
        
            updateRefrCellISI;
        
        
            %%%% Update figures & plots
            updateMD;
            updateRefractorySpikesView;
        end
    end


    function updateRefractorySpikesView
        
        cell_x = rot2D_projFeatures(1,:); 
        cell_y = rot2D_projFeatures(2,:); 
%         cellColor = cellColors_dens(refr_cellIdx,:);
        cellColor = 'b';
        set(h_2Dr_spks, 'xdata', cell_x, 'ydata', cell_y, 'color', cellColor, 'markersize', 3);        
        
        switch refrWaveformCategory
            case 'Refractory spikes',                
                [refr_x,refr_y] = linesFromAtoB(cell_x(refrCell_spkIdxs), cell_y(refrCell_spkIdxs), cell_x(refrCell_spkIdxs+1), cell_y(refrCell_spkIdxs+1));
                set(h_2Dr_refr, 'xdata', refr_x(:), 'ydata', refr_y(:), 'visible', 'on', 'markersize', 6 );        
                set(h_2Dr_spks, 'userData', 'skipSelect');
                set(h_refr_xy_cur, 'visible', 'on')
            case 'All spikes'
                set(h_2Dr_refr, 'visible', 'off')          
                set(h_2Dr_spks, 'userData', ''); 
                set(h_refr_xy_cur, 'visible', 'off')
        end
        
        if isempty(refrSpkInRangeSelect_idx) || isempty(refrSpkSelect_idx) %%% double check (that need both)
            set(h_2Dr_refrSel, 'visible', 'off');  
            set(hnds_refr_cutButtons, 'enable', 'off');
            
        else

            switch refrWaveformCategory
                case 'Refractory spikes',

                    %%%% Refractory Spike XY plot
                    spkIds = refrCell_spkIdxs(refrSpkInRangeSelect_idx)+ curOffset_add;
                    idx_cellMarked = refrSpikeMarkedIdxs{refr_cellIdx}{refrSpkSelect_idx};
                    if isempty(idx_cellMarked)
                        idx_cellMarked = getIdxSpikesBeyondSelSpike([cell_x(:), cell_y(:)], spkIds);
                    end
                    
                    %%%% Rotated Cluster plot
                    set(h_2Dr_spksToDel, 'xdata', cell_x(idx_cellMarked), 'ydata', cell_y(idx_cellMarked), 'visible', 'on');
                    set(h_2Dr_spksDel, 'xdata', cell_x(clustSpikesRemoved_idx{refr_cellIdx}), 'ydata', cell_y(clustSpikesRemoved_idx{refr_cellIdx}), 'visible', 'on',....
                        'markerfacecolor', 'none');

                    %             set(h_2Dr_tmp, 'xdata', [0 A(1), B(1)], 'ydata', [0 A(2), B(2)], 'visible', 'off');
                    nCutHere = length(idx_cellMarked);
                    nCutTotal = length(clustSpikesRemoved_idx{refr_cellIdx});
                    curCutStr = sprintf('This cut:  %d / %d (%.2f%%)', nCutHere, length(cell_x) , nCutHere/length(cell_x)*100);
                    totalCutStr = sprintf('All cuts:  %d / %d (%.2f%%)', nCutTotal, length(cell_x) , nCutTotal/length(cell_x)*100);
                    
                    set([h_2Dr_xlab h_md_xlab], 'string', {curCutStr, totalCutStr} );
                    set([h_2Dr_tit h_md_tit], 'string', sprintf('Cell # %d.  Refractory spike pair # %d (%d in range, %d total). [%d spikes]', ...
                        cellIds(refr_cellIdx), refrSpkSelect_idx, length(refrSpksInRange_list{refr_cellIdx}), nRefrSpikes(refr_cellIdx), cellSpikeCounts(refr_cellIdx) ));

%                     set([h_2Dr_xlab h_md_xlab], 'string', {curCutStr} );
%                     set([h_2Dr_tit h_md_tit], 'string', sprintf('Cluster # %d.  Refractory spike pair # %d / %d', ...
%                         cellIds(refr_cellIdx), refrSpkSelect_idx, nRefrSpikes(refr_cellIdx)));

                    set(hnds_refr_cutButtons, 'enable', 'on');
                    
                case 'All spikes',
                    nCutTotal = length(clustSpikesRemoved_idx{refr_cellIdx});
                    totalCutStr = sprintf('All cuts:  %d / %d (%.2f%%)', nCutTotal, length(cell_x) , nCutTotal/length(cell_x)*100);
                    
                    set([h_2Dr_xlab h_md_xlab], 'string', totalCutStr );
                    set([h_2Dr_tit h_md_tit], 'string', sprintf('Cell # %d. Spike # %d / %d (%d Refr spikes in range)', ...
                        cellIds(refr_cellIdx), refrSpkSelect_idx, cellSpikeCounts(refr_cellIdx), length(refrSpksInRange_list{refr_cellIdx}) ));
                    
                    set(hnds_refr_cutButtons, 'enable', 'off');
            end
            
            
        end
                          
        
        allX = refrSpikesInfo(refr_cellIdx).NFurtherFromCent /cellSpikeCounts(refr_cellIdx)*100;
        allY = refrSpikesInfo(refr_cellIdx).minSpikePeaks;            
        
        sizeData = ones(1,nRefrSpikes(refr_cellIdx))*6^2;
        sizeData(curCellEdits.idx_cutOut) = 2^2;
                        
        set(h_refr_xy,     'xdata', allX(curCellEdits.idx_0),   'ydata', allY(curCellEdits.idx_0), 'sizedata', sizeData(curCellEdits.idx_0) );        
        set(h_refr_xy_del, 'xdata', allX(curCellEdits.idx_del), 'ydata', allY(curCellEdits.idx_del), 'sizedata', sizeData(curCellEdits.idx_del) );
        set(h_refr_xy_cut, 'xdata', allX(curCellEdits.idx_cut), 'ydata', allY(curCellEdits.idx_cut), 'sizedata', sizeData(curCellEdits.idx_cut) );

        
        arrayfun(@(v) setButtonBkgColor(v), [h_delSpk_button, h_applyCut_button, h_discardClust_button]);        

        % update colors of cellblocks.
%         cellBlock_cols
        
        cellBlock_cols(2,:) = log(1 + nRefrSpikesRemaining(:)');
        cellBlock_cols(2,refrClustAction == 1 | uSpkCellIds == 0) = nan;
        set(h_cellBlock_im, 'CData', cellBlock_cols);
                
        % update background color of xy plot, if necessary
        bkCol = iff(refrClustAction(refr_cellIdx) == 1, [.9 .9 .9], [1 1 1]);
        set(h_refr_xy_ax, 'color', bkCol);

        
    end

    function updateSelectedRefractorySpike
    
        if ~isempty(refrSpkInRangeSelect_idx) 
            updateRefractorySpikeIdx;            
%             refrSpkInRangeSelect_idx = bound(refrSpkInRangeSelect_idx, 1, length(refrSpksInRange_list{refr_cellIdx}));
%             refrSpkSelect_idx = refrSpksInRange_list{refr_cellIdx}(refrSpkInRangeSelect_idx);
            
%             refrCell_spkIdxs = allCell_refrSpkIdxs{refr_cellIdx};            
            
%             set(h_delSpk_button(1), 'value', edits.refrSpikesDeleted{refr_cellIdx}(refrSpkInRangeSelect_idx,1));
%             set(h_delSpk_button(2), 'value', edits.refrSpikesDeleted{refr_cellIdx}(refrSpkInRangeSelect_idx,2));
    %         h_applyCut_button
            
            curSpkIdxs = refrCell_spkIdxs(refrSpkInRangeSelect_idx)+curOffset_add;            
            idxSelectedSpikes = cellSpikeIdx{refr_cellIdx}(curSpkIdxs);
            posSelectedSpikes = spikeTimes_tk(idxSelectedSpikes);            
                        
            r2D_autoRotate = get(h_autoRot_chkbox, 'value') && strcmp(refrWaveformCategory, 'Refractory spikes');    
            if r2D_autoRotate
                % project all features onto plane spanned by these two spikes.
                rot2D_projFeatures = projFeaturesOntoSpikes(rot2D_allFeatures, curSpkIdxs, 2);
            end            
            cell_x = rot2D_projFeatures(1,:);
            cell_y = rot2D_projFeatures(2,:);         
            
            if strcmp(refrWaveformCategory, 'Refractory spikes')
                set(h_2Dr_refrSel(1), 'xdata', cell_x(curSpkIdxs), 'ydata', cell_y(curSpkIdxs), 'visible', 'on');              
                set(h_2Dr_refrSel(2), 'xdata', cell_x(curSpkIdxs(2)), 'ydata', cell_y(curSpkIdxs(2)), 'visible', 'on');  
                curInfo = refrSpikesInfo(refr_cellIdx);
                set(h_refr_xy_cur, 'xdata', curInfo.NFurtherFromCent(refrSpkSelect_idx) /cellSpikeCounts(refr_cellIdx)*100, ...
                                   'ydata', curInfo.minSpikePeaks(refrSpkSelect_idx), ...
                                   'visible', 'on', 'markersize', 11 );    
%                  set(h_refr_xy_tit, 'string', sprintf('isis = %.2f. t1 = %.2f', curInfo.isi_ms(refrSpkSelect_idx), curInfo.spkTs_ms(refrSpkSelect_idx)) );
%                 assert(diff(double(posSelectedSpikes))/ticksPerMs == refrSpikesInfo(refr_cellIdx).isi_ms(refrSpkSelect_idx));
                 set(h_refr_xy_tit, 'string', '');
                               
                curAction = refrSpkAction{refr_cellIdx}(refrSpkSelect_idx);

                set(h_delSpk_button(1), 'value', curAction == 1 || curAction == 3 );        
                set(h_delSpk_button(2), 'value', curAction == 2 || curAction == 3 );
                set(h_applyCut_button,  'value',  curAction == -1 );
                               
            else
                set(h_2Dr_refrSel(1), 'xdata', cell_x(curSpkIdxs(1)), 'ydata', cell_y(curSpkIdxs(1)), 'visible', 'on');              
                set(h_2Dr_refrSel(2), 'visible', 'off');                                  
            end
    

%             if isempty(refrSpikeMarkedIdxs{refr_cellIdx}{refrSpkInRangeSelect_idx})  % most likely because was large, so didn't bother doing it before.
%                     3;
%             end           
                           
            % plot waveforms.
            
            updateISIwaveforms;

            updateRefractorySpikesView;

        else
            [posSelectedSpikes, idxSelectedSpikes] = deal([]);
            
            set([h_2Dr_tit h_md_tit], 'string', sprintf('Cell # %d.  No refractory spikes [%d spikes]', cellIds(refr_cellIdx), cellSpikeCounts(refr_cellIdx)) );
            set([h_2Dr_xlab h_md_xlab], 'string', '');            
            set(h_refr_xy_cur, 'visible', 'off');
            
            if strcmp(refrWaveformCategory, 'All spikes')
                updateISIwaveforms;
            end
            set([h_refrWvfm; h_refrMahala], 'visible', 'off');
        end
                    
                               
    end
    
    % Change the current RefrCell 
    function updateRefrCell(~, ~)
                
        refr_cellIdx = round(get(h_refrCellSelect, 'value'));
        refrCellChanged = refr_cellIdx_prev ~= refr_cellIdx;
                
        if refrCellChanged
            % update Cell # display
            set(h_refrCellDisplay, 'string', num2str(cellIds(refr_cellIdx)));
            
            % update ISI popump-menus
            set(h_refrISI_period_slider, 'value', cellsRefrPeriod_ms(refr_cellIdx));
            set(h_refrISI_period_box, 'string', sprintf('%.2f ms', cellsRefrPeriod_ms(refr_cellIdx))  );

%             if isempty(allCell_refrSpkIdxs{refr_cellIdx}) || isempty(refrSpksInRange_list{refr_cellIdx})
%                 refrSpkInRangeSelect_idx = zeros(1,0);                    
%             else
%                 refrSpkInRangeSelect_idx = refrSpksInRange_list{refr_cellIdx}(1);
%             end           
                    
            posSelectedSpikes = [];
            set(h_2Dr_spksToDel, 'visible', 'off')
            set(h_2Dr_spksDel, 'visible', 'off')

            rot2D_allFeatures = fetInEigBasis{refrBasisId, refr_cellIdx};
            rot2D_projFeatures = rot2D_allFeatures;
            
            set(h_cellBlock_select, 'position', [.5+(refr_cellIdx-1), cb_upperLim_ylim(1), 1, cb_upperLim] );                        
                        
            bkCol = iff(refrClustAction(refr_cellIdx) == 1, [.9 .9 .9], [1 1 1]);
            set(h_refr_xy_ax, 'color', bkCol);
            set(h_discardClust_button, 'value', refrClustAction(refr_cellIdx) ); 
            setButtonBkgColor(h_discardClust_button);
            if ~any(size(refrSpkInRangeSelect_idx))
                3;
            end
            
            curMeanWvfm_raw = meanWaveforms_raw(:,:,refr_cellIdx)';
            
            updateListOfCellSpikes;
            updateRefractorySpikeIdx;
            
            updateRefractorySpikes_curClust;                        
            updateSelectedRefractorySpike;    
            updateISIwaveforms;

            updateRefrCellCCG;
            refr_cellIdx_prev = refr_cellIdx;
        end
                
        updateRefrCellISI;
        updateISIwaveforms;
        
    end


    function clickSelectCell(~,~)
        cp = get(h_cellBlock_ax, 'currentPoint');
        newCellIdx = ceil(cp(1,1)-.5);
        set(h_refrCellSelect, 'value', newCellIdx);
        updateRefrCell;        
        
    end

    % Action Button #1, #2
    function deleteRefractorySpike(src, ~)
        
        if isempty(refrSpkSelect_idx)
            return;
        end
        assert( refrSpksInRange{refr_cellIdx}(refrSpkSelect_idx) );        
        
        del_id = get(src, 'userData');
        onOff = get(src, 'value');
        turningOn = onOff == 1;
        prevAction = refrSpkAction{refr_cellIdx}(refrSpkSelect_idx);
              
%         refrSpikesDeleted{refr_cellIdx}(refrSpkInRangeSelect_idx, id)        
%         curSpkIdxs = refrCell_spkIdxs(refrSpkInRangeSelect_idx)+[0, 1];        
        
        if turningOn && ((prevAction == -1) ||  (get(h_applyCut_button, 'value') == 1)) % 'cut off' was selected --> unselect it
            set(h_applyCut_button, 'value', 0)
        end

        if turningOn
            if ibetween(prevAction, 1, 2)  % 1/2 --> 3.
                refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = prevAction+del_id;
            else
                refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = del_id;
            end                        
        else % turning off 
            if (prevAction == 3)  % 3-->1/2.
                refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = prevAction-del_id;
            else
                refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = 0;
            end                        
        end        
        setButtonBkgColor(h_delSpk_button(del_id));        
                            
        updateRefractorySpikes_curClust;
    end
    
    % Action Button #3
    function cutRefractorySpikes(~, ~)
        
        if isempty(refrSpkInRangeSelect_idx)
            return;
        end
        assert( refrSpksInRange{refr_cellIdx}(refrSpkSelect_idx) );        
        
        onOff = get(h_applyCut_button, 'value');
                       
        prevAction = refrSpkAction{refr_cellIdx}(refrSpkSelect_idx);
        if onOff % turned from off-->on
            if prevAction > 0  % one or both of the spikes were deleted --> undelete them.
                set(h_delSpk_button, 'value', 0); 
                
            end
            if isempty(refrSpikeMarkedIdxs{refr_cellIdx}{refrSpkSelect_idx})  % most likely because was large, so didn't bother doing it before.
                refrSpikeIds = allCell_refrSpkIdxs{refr_cellIdx}(refrSpkSelect_idx)+[0,1];
                projFeatures = projFeaturesOntoSpikes(fetInEigBasis{refrBasisId, refr_cellIdx}, refrSpikeIds, 1);
                idx_marked_spks = getIdxSpikesBeyondSelSpike(projFeatures', refrSpikeIds);                
                refrSpikeMarkedIdxs{refr_cellIdx}{refrSpkSelect_idx} = idx_marked_spks;                                
            end

            
%             idxSpksDeleted{refr_cellIdx} = unique([idxSpksDeleted{refr_cellIdx}(:); idx_cellMarked(:)]);
            refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = -1;
        else  % turned from on-->off
%             idxSpksDeleted{refr_cellIdx} = setdiff(idxSpksDeleted{refr_cellIdx}, idx_cellMarked(:));
            refrSpkAction{refr_cellIdx}(refrSpkSelect_idx) = 0;
        end
        setButtonBkgColor(h_applyCut_button);
        updateRefractorySpikes_curClust;
        
%         idx_0 = refrSpkAction{refr_cellIdx} == 0;
%         idx_del = refrSpkAction{refr_cellIdx}(refrSpkInRangeSelect_idx) >= 1;
%         idx_cut = refrSpkAction{refr_cellIdx}(refrSpkInRangeSelect_idx) == -1;
        
    end


    % Action Button #4
    function discardCluster(~, ~)
        onOff = get(h_discardClust_button, 'value');                
        setButtonBkgColor(h_discardClust_button);
        refrClustAction(refr_cellIdx) = onOff;
        
        updateRefractorySpikesView;                     
    end

    % Action Button #5
    function autoCut_callback(~,~)
        autoCutRefractorySpikes(refr_cellIdx);
    end

    function autoClust_callback(~, ~)
        [assignClustToMU, refrPeriod_ms] = autoClassifyCluster(refr_cellIdx);        

        % adjust refrPeriod slider / box
        set(h_refrISI_period_slider, 'value', refrPeriod_ms);
        set(h_refrISI_period_box, 'string', sprintf('%.2f ms', refrPeriod_ms));

        % adjust 'discard cluster' button
        set(h_discardClust_button, 'value', assignClustToMU);
        discardCluster;        
    end
        
    function [assignClustToMU, refrPeriod_ms, data] = autoClassifyCluster(cellIdx)
        
        [assignClustToMU, refrPeriod_ms, data] = getClusterAutoClassification(cellIdx);        

        cellsRefrPeriod_ms(cellIdx) = refrPeriod_ms;            
        autoCutRefractorySpikes(cellIdx, refrPeriod_ms);
        
        refrClustAction(refr_cellIdx) = assignClustToMU;        
        updateRefractorySpikes_curClust;
    end



    function spkAction = getAutoCutForCluster(cellIdx, refrPeriod_ms)

        spkAction = zeros(1, nRefrSpikes(cellIdx));
        
        curDep = refrSpikesInfo(cellIdx).dep;        
        if (nRefrSpikes(cellIdx) == 0) || (isscalar(curDep) && isnan(curDep));
            return;
        end
        
        if (nargin < 2) || isempty(refrPeriod_ms)
            refrPeriod_ms = cellsRefrPeriod_ms(cellIdx);
        end
        idx_inRange = find(refrSpikesInfo(cellIdx).isi_ms <= refrPeriod_ms);            
        
        N = length(idx_inRange);
        if N > 0
            curDep_inRange = curDep(idx_inRange, idx_inRange);

            spk_i = N;
            spks_cut = zeros(1,N);
            while ~isempty(spk_i);
                spks_cut(spk_i) = 1;            
                spks_cut(curDep_inRange(:,spk_i)==1) = -1;
                spk_i = find(spks_cut==0, 1, 'last');
            end
            spkAction(idx_inRange(spks_cut == 1)) = -1;
            spkAction(idx_inRange(spks_cut ~= 1)) = 0;                
        end
        
    end

    function autoCutRefractorySpikes(cellIdx, refrPeriod_ms)
        if nargin < 2
            refrPeriod_ms = cellsRefrPeriod_ms(cellIdx);
        end
        
        depMtx = refrSpikesInfo(cellIdx).dep;
        isis_ms = refrSpikesInfo(cellIdx).isi_ms;
        
        spkAction = double( getAutoCutFromDepMtx(depMtx, isis_ms, refrPeriod_ms) );
        spkAction(spkAction == 1) = -1; % -1 == cut in this context.
        
%         if isempty(spkAction)
%             spkAction = zeros(1,0);
%         end        
        refrSpkAction{cellIdx} = spkAction;
        
        updateRefractorySpikes_curClust;        
    end




    function [assignToMU, bestRefrPeriod_ms, data] = getClusterAutoClassification(cellIdx)        
                
%         cutRemainTh = 0.66;
%         meanWvfmSecondPeakRatio = 0.5;
%         refrRange_ms = [.8:.05:1.2];
%         refrCut_pow = 2;
        
        mu_cluster = cellIds(cellIdx) == 0;
                

        % 	1a. [marked as MU from # refr violations]
        % 	1. autocut: remaining < 66%        
        
        true_refrPeriod_ms = [];
%         if useICrefrPeriod
            IC_cand_idx = find(IC_stats.IC_candidates == cellIds(cellIdx));
%             EC_clust_idx = find(IC_stats.EC_clustIds == cellIds(cellIdx), 1);
%             fracOfClusterIC = IC_stats.fracOfCluster_IC(EC_clust_idx);
            
            if useICrefrPeriod_useIndivClustRefrPeriod && ~isempty(IC_cand_idx) %% && (fracOfClusterIC > 0.4)                                            
                true_refrPeriod_ms = IC_stats.IC_refrPeriod_EC_clust_ms(IC_cand_idx);
            else
                true_refrPeriod_ms = IC_stats.IC_refrPeriod_all_EC_ms;
            end
            if isnan(true_refrPeriod_ms)
                true_refrPeriod_ms = [];
            end
            
            autoClassOpt.true_refrPeriod_ms = true_refrPeriod_ms;
            if useICrefrPeriod
                true_refrPeriod_ms_arg = true_refrPeriod_ms;
            else
                true_refrPeriod_ms_arg = [];
            end
%         end

        clust_refrMaxCalc = ClustData.isiData.refrMaxCalc(cellIdx);
        
        
        
        [mu_refrac, refrFracRemoved, bestRefrPeriod_ms, allFracsRemoved] = isMUfromRefracSpikes(...
            cellIdx, autoClassOpt.cutRemainTh, autoClassOpt.refrRange_ms, ...
            autoClassOpt.refr_leeway_ms, clust_refrMaxCalc, true_refrPeriod_ms_arg);

        
        
        % 	1. MD: (a) mode of MDs is close to threshold
        %            (b) std dev is > std(MD_all)
        [mu_mahala, idxMode, stdMD]  = isMUfromMahala(cellIdx, autoClassOpt.nFirstBins_mode);

        %  2. mean spike waveform: multiple negative peaks
        %    -> second peak (after one at 0) is > 50% of 0 peak	
        [mu_wvfm, relHgts]  = isMUfromMeanWvfm(cellIdx, autoClassOpt.meanWvfmSecondPeakRatio);

%         assignToMU = mu_cluster || mu_refrac || mu_mahala || mu_wvfm;
        assignToMU = mu_cluster || mu_refrac || ...
            (mu_mahala && autoClassOpt.assignToMU_mahalaSmallOrWide) || ...
            (mu_wvfm && autoClassOpt.assignToMU_waveformMultPeaks);        

        data = struct('Gid', Gid, 'clustId', cellIds(cellIdx), 'nspikes', cellSpikeCounts(cellIdx), 'refrFracRemoved', refrFracRemoved, 'stdMD', stdMD, 'relHgts', relHgts, 'bestRefrPeriod_ms', bestRefrPeriod_ms, 'allFracsRemoved', single(allFracsRemoved));
            
    end



    function [data, assignToMU] = autoClassifyAllClusters(~,~)
        
        initialRun = nargin == 0;
        
        data(nCells) = struct('Gid', [], 'clustId', [], 'nspikes', [], 'refrFracRemoved', [], 'stdMD', [], 'relHgts', [], 'bestRefrPeriod_ms', [], 'allFracsRemoved', []);
        
        assignToMU = false(1,nCells);
        bestRefrPer_ms = zeros(1,nCells);
        assignToMU(1) = true;        
        for ci = 1:nCells            
            [assignToMU(ci), bestRefrPer_ms(ci), data(ci)] = getClusterAutoClassification(ci);
            if ~initialRun  % on initial run, just get whether MU or not.
                            % if called, actually apply the cuts & MU assignments
                cellsRefrPeriod_ms(ci) = bestRefrPer_ms(ci);
                autoCutRefractorySpikes(ci);                
                
                refrClustAction(ci) = assignToMU(ci);
                
                allCellsFracsRemoved{ci} = data(ci).allFracsRemoved;
            end
        end        
        
        if showAllFigures
            autoClust_callback; % this is to update the sliders, etc. for the current cluster.
        end
        
        cellBlock_cols(3,assignToMU) = nan;
        cellBlock_cols(3,~assignToMU) = 0;    
        updateRefractorySpikes_allClusters;
        
    end



    
    function [isMU, fracRemoved, bestRefrPeriod_ms, allFracsRemoved] = isMUfromRefracSpikes(cellIdx, th, refrRange_ms, ...
            refr_leeway_ms, clust_refrMaxCalc, true_refrPeriod_ms)
        
        depMtx = refrSpikesInfo(cellIdx).dep;
        isis_ms = refrSpikesInfo(cellIdx).isi_ms;
               
        params.nSpksInClust = cellSpikeCounts(cellIdx);
        params.refrMaxCalc = clust_refrMaxCalc;
        params.true_refrPeriod_ms = true_refrPeriod_ms;               

        autoClassOpt.clustId = cellIds(cellIdx);
        autoClassOpt.Gid = Gid;
        [bestRefrPeriod_ms, fracRemoved, isMUfromDep, allFracsRemoved] = getBestRefrPeriodFromDepMtx(...
            depMtx, isis_ms, refrSpikeMarkedIdxs{cellIdx}, params, autoClassOpt);

        
        isMU = isMUfromDep || (1-fracRemoved < th);
        
        %{
            
            h2 = subplot(3,1,2);
            binE = linspace(0, refrRange_ms(end), 40)
            binV = histcnt(isis_ms, binE);
            binC = binEdge2cent(binE);
            bar(binC, binV, 1);
            xlims = get(h2, 'xlims');
            xlim(xlims)

        %}        
               
        return;
%         
%         bestRefrPeriod_ms = refrRange_ms(end);
%         N = length(refrRange_ms);
%         if (nRefrSpikes(cellIdx) == 0) % no refractory spikes.
%             isMU = false;
%             fracRemoved = 0;
%             return;
%         end
%             
%         isMultiUnit = isscalar(dep) && isnan(dep);
% %         isMultiUnit = ischar(dep) && strcmp(dep, 'Too many refractory spikes') 
% %         isMultiUnit = isempty(dep);
% 
%         if isMultiUnit % too many refractory spikes.
%             isMU = true;
%             fracRemoved = nan;
%             return;
%         end
%                   
%         fracsRemoved = zeros(1,N);
%         for ri = 1:N
%             spkAction = getAutoCutForCluster(cellIdx, refrRange_ms(ri) );        
%             idx_cut  = find( spkAction == -1 );                
%             spikeCuts = unique(cat(1, refrSpikeMarkedIdxs{cellIdx}{idx_cut}));        
% 
%             nSpikesRemovedInAutoCut = length(spikeCuts);
%             nSpikesTotal = cellSpikeCounts(cellIdx);
%             fracsRemoved(ri) = (nSpikesRemovedInAutoCut / nSpikesTotal);            
%         end
%         
%         refr_pow = 2;
%         gap_ms = 0.2;
%         m_expected = fracsRemoved(end) / (refrRange_ms(end) - gap_ms);
%         c_expected = fracsRemoved(end) - m_expected * refrRange_ms(end);
%         fracsRemoved_expected = m_expected * refrRange_ms + c_expected;
%         
%         idx_below = fracsRemoved <= fracsRemoved_expected;
%         eps = 1e-5;
%         drem_rel = fracsRemoved ./ fracsRemoved(end);
%         drem_abs = fracsRemoved - fracsRemoved(1);
%         drem = (fracsRemoved_expected - fracsRemoved) + eps;
%         
% %         drem = drem_abs .* drem_rel;
% %         drem = (fracsRemoved_expected - fracsRemoved);
%         dt   = (refrRange_ms(end) - refrRange_ms) + eps;
%         
% %         refrDecreases = idx_below .* (abs(drem).^(refr_pow) ./ dt);
%         refrDecreases = idx_below .* (drem);
%         idx_best_tmp = indmax(refrDecreases);
%         
%         % this often results in a refractory period that is slightly too restrictive - allow for a
%         % little leeway:
%         idx_1ms = find(refrRange_ms == 1, 1);
%         
%         nLeewayPts = round( refr_leeway_ms/diff(refrRange_ms(1:2)) );
%         nLeewayPts = min(nLeewayPts, rectified(idx_best_tmp-idx_1ms) );
%         
%         fracsRemoved_nearBest = fracsRemoved(idx_best_tmp + [-nLeewayPts:0]);
%         bestFracRemoved = min ( fracsRemoved_nearBest );
%         idx_best = find(bestFracRemoved == fracsRemoved_nearBest, 1, 'last') + idx_best_tmp - (nLeewayPts+1);
%                 
%         bestRefrPeriod_ms = refrRange_ms(idx_best);
%         fracRemoved = fracsRemoved(idx_best);        
        

    end

    function [isMU, idxMode, stdMD] = isMUfromMahala(cellIdx, nFirstBins_mode)
                    
        curCellMDs = spikeMahalaHeights(cellSpikeIdx{cellIdx});
        curCellMDs_bin = histcnt(curCellMDs, mdBinEdges);
        curMDs_sm = gaussSmooth(curCellMDs_bin, 1.5);

        stdMD = std(curCellMDs);
        
        idxMode = indmax(curMDs_sm);
        isMU = (idxMode <= nFirstBins_mode) || (stdMD > std(spikeMahalaHeights));
        
    end

    function [mu_wvfm, relHgts] = isMUfromMeanWvfm(cellIdx, th_ratio)
       3;
        mnWvfm = bsxfun(@times, meanWaveforms_raw(:,:,cellIdx), channelScaling');
        idx_t0 = indmin(abs(t_ms-0));
        mu_channel = false(1,nChannels);
        relHgts = zeros(1,nChannels);
        
        for ch_i = 1:nChannels
            idx_mins = findLocalMinima(mnWvfm(:,ch_i), 3);
            ind_0min = indmin(abs(idx_mins - idx_t0));                 
            val_at0 = mnWvfm(idx_mins(ind_0min),ch_i);
            if (val_at0 > -3) % should be at least 3 std devs below baseline
                continue;
            end
            
            nMins = length(idx_mins);                        
            
            if nMins > 1                
                for min_i = setdiff(1:nMins, ind_0min)
                    val = mnWvfm(idx_mins(min_i),ch_i);
                    if (val < 0)
                        curRatio = val / val_at0;
                        relHgts(ch_i) = max(relHgts(ch_i), curRatio);
%                         if (curRatio >= th_ratio)
%                             mu_channel(ch_i) = true;
%                         end
                    end                    
                end                                
            end            
            
        end        
%         mu_wvfm = any(mu_channel);              
        mu_wvfm = mean(relHgts) >= th_ratio;              
        
    end




    %%%% Other Action Buttons : Loading, Saving, etc %%%
    
    %%% Action button 2,1: LOAD
    function loadEdits(~, ~) %#ok<*INUSD>                
        if exist(clusterPruningsFile, 'file')                        
            S = load(clusterPruningsFile);
            savedBasisFeatureSet = S.refrBasisFeatureSet;
            if ~strcmp(featureSetOptions{refrBasisId}, savedBasisFeatureSet) 
                error('Saved file was generated using a different basis feature set');
            end
            if ~isequal(S.allCell_refrSpkIdxs, allCell_refrSpkIdxs)
                error('Different refractory spikes data detected');
            end
            
%             [  refrSpkAction,    refrClustAction,   clustSpikesRemoved_idx,   cellsRefrPeriod_ms,   nRefrSpikesRemaining] = deal(...
%              S.refrSpkAction, S. refrClustAction, S.refrSpikesRemoved , S.cellsRefrPeriod_ms, S.nRefrSpikesRemaining);
            [  refrSpkAction,    refrClustAction,   clustSpikesRemoved_idx,   cellsRefrPeriod_ms,   nRefrSpikesRemaining] = deal(...
             S.refrSpkAction, S. refrClustAction, S.clustSpikesRemoved_idx , S.cellsRefrPeriod_ms, S.nRefrSpikesRemaining);
         
        end
        
        updateRefractorySpikes_curClust;
    end

    % Action Button 2,2: SAVE
    function saveCurrentEdits(~, ~) %#ok<*INUSD>
        saveStruct = struct('allCell_refrSpkIdxs', {allCell_refrSpkIdxs}, ...
                             'refrSpkAction', {refrSpkAction}, ...
                             'refrClustAction', refrClustAction, ...
                             'clustSpikesRemoved_idx', {clustSpikesRemoved_idx}, ...
                             'cellsRefrPeriod_ms', cellsRefrPeriod_ms, ... 
                             'cellsRefrPeriod_range_ms', single(autoClassOpt.refrRange_ms), ...
                             'cellsFracsRemoved', {allCellsFracsRemoved}, ...
                             'refrBasisFeatureSet', featureSetOptions{refrBasisId}, ...
                             'nRefrSpikesRemaining', nRefrSpikesRemaining, ...
                             'clustIds', cellIds); %#ok<NASGU>
        
         save(clusterPruningsFile, '-struct', 'saveStruct');        
        
             
        if strcmp(clustGrouping, 'clusters'); % don't actually create the pruned cluster file for MergePruning or PruneMergePruning 
            generatePrunedClusterFile(Gid, clustGrouping);
        end
        
        verifyPruningClearedRefrPeriod(Gid, clustGrouping);
        3;
    end

    % Action Button #2,3: CLEAR
    function clearCurCellEdits(~, ~)
        refrSpkAction{refr_cellIdx}(:) = 0;
        refrClustAction(refr_cellIdx) = 0;
        
        set(hnds_refr_allButtons, 'value', 0);
        updateRefractorySpikes_curClust;
        
    end


    % Action Button #2,4: CLEAR ALL
    function clearAllCellEdits(~, ~)
        for ci = 1:nCells
            refrSpkAction{ci}(:) = 0;            
        end
        refrClustAction(:) = 0;  
        
        set(hnds_refr_allButtons, 'value', 0);
        updateRefractorySpikes_allClusters;
    end
    %%%%%% 


    % Current RefrCell ISI update #2
    function updateRefrCellISI(~, evnt)
        
        rng_idx = get(h_refrISI_range, 'value');
        nbin_idx = get(h_refrISI_nbins , 'value');

        isi_range_ms = isi_allRanges_ms(rng_idx);
        isi_nbins    = isi_allNbins(nbin_idx);
        isiBinEdges = linspace(0, isi_range_ms, isi_nbins+1);
        isiBinCent = binEdge2cent(isiBinEdges);

        refrISI_binOrig = masterBin2Bin(allISIs_masterBin(refr_cellIdx, refr_cellIdx), isi_masterBinIdxs{rng_idx, nbin_idx});                
        refrISI_binOrig = refrISI_binOrig{1};
        %         refrISI_bin = allISIs_bin{refr_cellIdx, refr_cellIdx};
                                                                                                                    
        xlims = [isiBinEdges(1), isiBinEdges(end)];
        
        remainCellSpkIdxs = cellSpikeIdx{refr_cellIdx};
        remainCellSpkIdxs(clustSpikesRemoved_idx{refr_cellIdx}) = [];
        curCell_ISI_val = getISIs( spikeTimes_ms(remainCellSpkIdxs), true, isi_range_max );
        discardLastBin = (isi_range_ms ~= isi_range_max);
        curCell_ISI_bin = histcnt( curCell_ISI_val, isiBinEdges, [], discardLastBin);                       
                
        showAllISIs = get(h_refrISI_showAll_chk, 'value');
        if showAllISIs        
            refrISI_binDiff = curCell_ISI_bin(:)' - refrISI_binOrig;
            refrISI_binRemoved = refrISI_binDiff; refrISI_binRemoved(refrISI_binRemoved > 0) = 0;
            refrISI_binAdded   = refrISI_binDiff; refrISI_binAdded( refrISI_binAdded < 0) = 0;
            refrISI_binRemaining = refrISI_binOrig + refrISI_binRemoved;

            set(h_refrCell_isi_bar(1), 'xdata', isiBinCent, 'ydata', refrISI_binRemaining);
            set(h_refrCell_isi_bar(2), 'xdata', isiBinCent, 'ydata', -refrISI_binRemoved, 'visible', 'on');
            set(h_refrCell_isi_bar(3), 'xdata', isiBinCent, 'ydata', refrISI_binAdded, 'visible', 'on');
        else
            set(h_refrCell_isi_bar(2:3), 'xdata', isiBinCent, 'ydata', zeros(size(isiBinCent)), 'visible', 'off');  % for some reason need to set x & y data, otherwise weird things happen... (?)
            set(h_refrCell_isi_bar(1), 'xdata', isiBinCent, 'ydata', curCell_ISI_bin);            
        end
%         txt_col = iff(refracViol_n(cell_I_idx, cell_J_idx) > 10 || refracViol_pct(cell_I_idx, cell_J_idx) >.5, 'r', 'k');
%         txt_prefix = iff(cell_I == cell_J, sprintf('[%d]', cell_I), sprintf('[%d,%d]', cell_J, cell_I));

%                     s1 = sprintf('%s<%.1fms;', txt_prefix, refrPeriod_ms);
%                     s2 = sprintf('%d spks, (%.1f%%)', refracViol_n(cell_I_idx, cell_J_idx), refracViol_pct(cell_I_idx, cell_J_idx));
        set(h_refrCell_isi_vline, 'visible', 'off');
        set(h_refrCell_isi_ax, 'xlim', xlims, 'xtickmode', 'auto', 'ytickmode', 'auto');
        ylims = get(h_refrCell_isi_ax, 'ylim');
                    
        set(h_refrCell_isi_vline, 'xdata', [1,1]*cellsRefrPeriod_ms(refr_cellIdx), 'ydata', ylims, 'visible', 'on', ...
            'color', 'k', 'linestyle', '--', 'linewidth', 2);
        
        refrISI_binOrig_tot = sum(refrISI_binOrig);
        curCell_ISI_bin_tot = sum(curCell_ISI_bin);        
%         set(h_refrCell_isi_tit, 'string', sprintf('Orig: %d. Cur: %d', refrISI_binOrig_tot, curCell_ISI_bin_tot));

        set(h_refrCell_isi_tit, 'string', 'ISI distribution');
        set(h_refrCell_isi_xlab, 'string', 'ms')
        set(h_refrCell_isi_ylab, 'string', '# spikes');
        3;
%         set(h_isi_tit(cell_I_idx,cell_J_idx), 'string', sprintf('ISI for cell %d', cellIds(cell_I_idx) ));
                           
%         set(h_isi_txt(cell_I_idx, cell_J_idx), 'position', [isi_range_ms/2, ylims(2)], 'string', {s2}, 'color', txt_col, 'visible', 'on');

        
        
    end
    
    % Current RefrCell ISI update #3
    function updateRefrCell_RefrPeriod(src, evnt)        
        
        switch get(src, 'tag')
            case 'refrPeriod_slider', 
                newVal = get(src, 'value');                
                                
            case 'refrPeriod_text',
                newVal = sscanf( get(src, 'string'), '%f' );
                if isempty(newVal) || isnan(newVal)
                    set(src, 'string', num2str( cellsRefrPeriod_ms(refr_cellIdx), '%.2f ms') );
                    return;
                end
        end
        if ~any(newVal == refrPeriod_ms_range)
            newVal = refrPeriod_ms_range( indmin(abs(newVal - refrPeriod_ms_range)) );
        end            
                
        set(h_refrISI_period_slider, 'value', newVal);
        set(h_refrISI_period_box, 'string', sprintf('%.2f ms', newVal) );
        set(h_refrCell_isi_vline, 'xdata', [1,1]*newVal)
        
        cellsRefrPeriod_ms(refr_cellIdx) = newVal;
        updateRefractorySpikes_curClust;
    end

    

    function updateISIwaveforms

        doCurSpikeWvfm = 1;
        doMeanWvfm = 1;
%         set(h_refrWvfm, 'marker', 'o', 'markerfacecolor', 'auto', 'markersize', 4)
%         for jj = 1:4
%             set(h_refrWvfm(jj), 'markerfacecolor', get(h_refrWvfm(jj), 'color'), 'linewidth', 1);
%         end
%         set(h_refrWvfm_ax, 'position', [60 60 350 220]);
%         set(h_refrWvfm_ax, 'box', 'on');
        
        
        if isempty(posSelectedSpikes)
            doCurSpikeWvfm = 0;
        end
        
        if doCurSpikeWvfm
            spikeDetectTh = 6;
            posSelectedSpikes = double(posSelectedSpikes);
            set([h_refrWvfm; h_refrMahala; h_refrSpkLines(:)], 'visible', 'on');

            switch refrWaveformCategory
                case 'Refractory spikes', 
                    % retrieve spike waveforms of the two adjacent spikes
                    [pos1, pos2] = dealV(posSelectedSpikes);
                    [idx1, idx2] = dealV(idxSelectedSpikes);
                    if pos1 > pos2
                        [pos1, pos2] = deal(pos2, pos1);
                        [idx1, idx2] = deal(idx2, idx1);
                    end
                    assert(pos2-pos1 <= length(t_tk_ext));
                    tkStrt = pos1-20;
                    tkEnd  = pos2+35;
                    nTk = tkEnd-tkStrt+1;
                    wvfm_tks = tkStrt:tkEnd;

            %         idxTExt0 = find(t_tk_ext == 0);
                    t_1 = pos1+t_tk_ext;
                    t_2 = pos2+t_tk_ext;        
                    idx_ok1 = find( ibetween(t_1, tkStrt, tkEnd) );
                    idx_ok2 = find( ibetween(t_2, t_1(idx_ok1(end))+1, tkEnd) );                
                    % make sure overlap correctly:
%                     assert( isequal( spkWaveforms_ext(idx_ok1(end), :, idx1), spkWaveforms_ext(idx_ok2(1)-1, :, idx2) ) );

                    wvfm_amps_raw = double( vertcat(spkWaveforms_ext(idx_ok1, :, idx1), spkWaveforms_ext(idx_ok2, :, idx2) )');                
                    nWvfms = 2;

                case 'All spikes';
                    pos1 = posSelectedSpikes(1);
                    idx1 = idxSelectedSpikes(1);
                    tkStrt = pos1-30;
                    tkEnd  = pos1+40;                
                    wvfm_tks = tkStrt:tkEnd;

                    t_1 = pos1+t_tk_ext;
                    idx_ok1 = find( ibetween(t_1, tkStrt, tkEnd) );
                    wvfm_amps_raw = double( spkWaveforms_ext(idx_ok1, :, idx1)' );
                    nWvfms = 1;                
            end
            wvfm_x_tks = wvfm_tks/ticksPerMs;
        else
            wvfm_amps_raw = zeros(4,0);
        end
                
        switch refrWaveformScale
            case 'CCW'
                mahala_scale = 1;
                wvfm_amps = channelWhitenMtx * wvfm_amps_raw;   % multiply by channel-whitening matrix
                curMeanWvfm = channelWhitenMtx * curMeanWvfm_raw;
                
            case 'Scaled';
                mahala_scale = 1;
                wvfm_amps = bsxfun(@times, wvfm_amps_raw, channelScaling);   % divide by channel std deviation
                curMeanWvfm = bsxfun(@times, curMeanWvfm_raw, channelScaling(:));
                
            case 'Raw';
                mahala_scale = 10;
                wvfm_amps = wvfm_amps_raw;
                curMeanWvfm = curMeanWvfm_raw;
                
        end

        if doCurSpikeWvfm
            % get Mahalanobis distances for each spike
            [aboveThreshold, mahalaDists] = aboveEllipsoidalThreshold(wvfm_amps_raw, spikeDetectTh, EvectMax, invCovMtx);
            mn_above = mahalaDists; mn_above(~aboveThreshold) = nan;
            mn_below = mahalaDists; 


            pos1_ms = pos1/ticksPerMs;
            offset_tk = pos1;
            offset_ms = offset_tk/ticksPerMs;% wvfm_x_tks(1); % vs. 0;
            wvfm_xlims = [wvfm_x_tks(1), wvfm_x_tks(end)]-offset_ms;

            % update plots with : waveforms, mndist, and spike line indicators.
            for ch_i = 1:4
                set(h_refrWvfm(ch_i), 'xdata', wvfm_x_tks-offset_ms, 'ydata', wvfm_amps(ch_i,:));            
            end        
    %         set(h_refrWvfm_tit, 'string', sprintf('Spikes # %d and %d. %.2f ms apart', double(pos2-pos1)/ticksPerMs ));
    %         set(h_refrWvfm_tit, 'string', sprintf('%.2f ms apart. First spike at tick %d (%10.2f ms)', double(pos2-pos1)/ticksPerMs, pos1, pos1/ticksPerMs ));
%             pos1_str = sec2hms(pos1 / samplingRateHz);
            pos1_str = sprintf('%.2f ms', pos1/ticksPerMs);
            if nWvfms == 1
                set(h_refrWvfm_tit, 'string', sprintf('Spike at %s', pos1_str) );
            elseif nWvfms == 2                
%                 set(h_refrWvfm_tit, 'string', sprintf('%.2f ms apart [First spk at %s]', double(pos2-pos1)/ticksPerMs, pos1_str ));
                set(h_refrWvfm_tit, 'string', sprintf('%.2f ms apart', double(pos2-pos1)/ticksPerMs));
            end            
            
            set(h_refrMahala(1), 'xdata', wvfm_x_tks-offset_ms, 'ydata', mn_below*mahala_scale, 'color', [.6 .6 .6], 'linewidth', 1);        
            set(h_refrMahala(2), 'xdata', wvfm_x_tks-offset_ms, 'ydata', mn_above*mahala_scale, 'color', 'k', 'linewidth', 1);

            set(h_wvfm_th, 'xdata', wvfm_xlims, 'ydata', spikeDetectTh*[1 1]*mahala_scale);        
            adjustY = ~get(h_refrWvfm_ylock, 'value');
            if adjustY            
                allY = [wvfm_amps(:); mahalaDists(:)*mahala_scale];
                y_range = lims(allY, .05);                
            else
                y_range = get(h_refrWvfm_ax, 'ylim');
            end

            line_y = [y_range(1)-10, spikeDetectTh, y_range(2)+10];
            set(h_refrSpkLines(1), 'xdata', pos1_ms.*[1,1,1]-offset_ms, 'ydata', line_y, 'color', 'k');
            if nWvfms == 2    
                pos2_ms = pos2/ticksPerMs;
                set(h_refrSpkLines(2), 'xdata', pos2_ms.*[1,1,1]-offset_ms, 'ydata', line_y, 'color', 'k', 'visible', 'on');
            else
                set(h_refrSpkLines(2), 'visible', 'off');
            end
            set(h_refrWvfm_ax, 'ylim', y_range, 'xlim', wvfm_xlims);

            hgts = spikeMahalaHeights( cellSpikeIdx{refr_cellIdx}( refrCell_spkIdxs(refrSpkInRangeSelect_idx)+curOffset_add ) );
            set(h_mahala_txt(1), 'position', [pos1_ms-offset_ms,y_range(2)], 'string', num2str( hgts(1), '%.1f' ), 'horiz', 'right' );
            if nWvfms == 2
                set(h_mahala_txt(2), 'position', [pos2_ms-offset_ms,y_range(2)], 'string', num2str( hgts(2), '%.1f' ), 'horiz', 'left', 'visible', 'on' );
            else
                set(h_mahala_txt(2), 'visible', 'off')
            end
            set(h_mahala_txt, 'visible', 'on');
        end        
        % update Mean Waveform
%         h_refrMeanWvfm_ax 
        if doMeanWvfm && ~strcmp(meanWaveformCurScale, refrWaveformScale) || (meanWaveformCurCellIdx ~= refr_cellIdx)
            for ch_i = 1:4
                set(h_refrMeanWvfm(ch_i), 'xdata', [1:nT], 'ydata', curMeanWvfm(ch_i,:))
            end
            axis(h_refrMeanWvfm_ax, 'tight')
            meanWaveformCurScale = refrWaveformScale;
            meanWaveformCurCellIdx = refr_cellIdx;
        end
        3;
        
    end


    function updateListOfCellSpikes(src, evnt)             
        if showAllFigures %% ugly
            refrWaveformCategory = waveformCategoryOptions{get(h_refrWvfm_category, 'value')};            
        else
            refrWaveformCategory = 'Refractory spikes';
        end
            
        
        switch refrWaveformCategory
            case 'Refractory spikes', 
                refrCell_spkIdxs = allCell_refrSpkIdxs{refr_cellIdx}(  refrSpksInRange_list{refr_cellIdx}  );
                curOffset_add = [0, 1];
            case 'All spikes';
                refrCell_spkIdxs = 1:cellSpikeCounts(refr_cellIdx);
                curOffset_add = 0;
        end
        
        if (nargin > 0) % was called by the 'category' popup box - update waveforms immediately 
                        % (otherwise, was called when change the cell # this will be called later);            
            updateRefractorySpikeIdx;
            updateRefractorySpikes_curClust;
            updateSelectedRefractorySpike;
            updateISIwaveforms;
        end
%         refrWaveformCategory_prev = refrWaveformCategory;
        
    end

    function updateWaveformScale(src, evnt)
        refrWaveformScale = waveformScaleOptions{get(src, 'value')};
        updateISIwaveforms;        
    end
    
    function updateMD(~,~)
        showMD = get(h_showMD_chkbox, 'value');        
        showAll = get(h_showAllMD_chkbox, 'value');
        showLog = get(h_logMD_chkbox, 'value');
        if showMD
            set([h_md_ax h_md_bar h_md_sm], 'visible', 'on');        
            set(h_md_ax, 'xtickmode', 'auto')
            set(h_2Dr_clust_ax , 'visible', 'off')
            
            curCellSpkIdxs = cellSpikeIdx{refr_cellIdx};
            
            curCellMDs = spikeMahalaHeights(curCellSpkIdxs);
            curCellMDs_bin = histcnt(curCellMDs, mdBinEdges);
            allMDs_remain_bin = allMDs_bin(:)-curCellMDs_bin(:);
            
            curCellMDs_removed = spikeMahalaHeights(curCellSpkIdxs(clustSpikesRemoved_idx{refr_cellIdx}));
            curCellMDs_removed_bin = histcnt(curCellMDs_removed, mdBinEdges);
            curCellMDS_remain_bin = curCellMDs_bin(:)-curCellMDs_removed_bin(:);
            
            if showLog
                curCellMDS_remain_bin = log10(1+curCellMDS_remain_bin);
                curCellMDs_removed_bin = log10(1+curCellMDs_removed_bin);
%                 allMDs_remain_bin = log10(1+allMDs_remain_bin);
                
                allMDs_remain_bin = log10(1+allMDs_bin)-log10(1+curCellMDs_bin);
                
            end
            
            curMDs_sm = gaussSmooth(curCellMDS_remain_bin, 1.5);
            set(h_md_sm, 'xdata', mdBinCent, 'ydata', curMDs_sm);
            
            set(h_md_bar(1), 'xdata', mdBinCent, 'ydata', curCellMDS_remain_bin);
            set(h_md_bar(2), 'xdata', mdBinCent, 'ydata', curCellMDs_removed_bin);
            if showAll
                set(h_md_bar(3), 'xdata', mdBinCent, 'ydata', allMDs_remain_bin);
            else
                set(h_md_bar(3), 'xdata', mdBinCent, 'ydata', zeros(size(allMDs_remain_bin)));
            end
            set(h_md_ylab, 'string', sprintf('stddev = %.2f (%.2f)', std(curCellMDs), std(spikeMahalaHeights) ));
            
        else            
            set([h_md_ax h_md_bar h_md_sm], 'visible', 'off');        
            set(h_2Dr_clust_ax , 'visible', 'on')
%             set([h_2Dr_clust_ax h_2Dr_spks, h_2Dr_refr, h_2Dr_refrSel, h_2Dr_spksToDel, h_2Dr_spksDel], 'visible', 'on')
            
        end
        
        
    end

    function updateRefractorySpikeIdx
        
        if isempty(refrCell_spkIdxs)
            [refrSpkInRangeSelect_idx, refrSpkSelect_idx] = deal([]);
        else
            if isempty(refrSpkInRangeSelect_idx)
                refrSpkInRangeSelect_idx = 1;
            end
            refrSpkInRangeSelect_idx  = bound(refrSpkInRangeSelect_idx,  1, length(refrCell_spkIdxs) );
            switch refrWaveformCategory
                case 'Refractory spikes', 
                    refrSpkSelect_idx = refrSpksInRange_list{refr_cellIdx}(refrSpkInRangeSelect_idx);
                case 'All spikes'
                    refrSpkSelect_idx = refrSpkInRangeSelect_idx;
            end
            
        end
        
    end

    combos3D = nchoosek([1:nChannels],3);
    if show3D
        basis_id0 = 1;
        x_tet = 1;  xs = features{basis_id0}(:,x_tet);
        y_tet = 2;  ys = features{basis_id0}(:,y_tet);
        z_tet = 3;  zs = features{basis_id0}(:,z_tet);
        
        figure(fig_3D); set(fig_3D, 'Name', 'Features in 3D', 'NumberTitle', 'off'); clf; hold on;
        h_cell3D = zeros(1,nCells);        

        for comb_i = 1:size(combos3D,1);
            [t1, t2, t3] = deal(combos3D(comb_i,1), combos3D(comb_i,2), combos3D(comb_i,3));
            for bi = 1:nBases
                axesLims3D{bi}(comb_i,:) = [fet_mins{bi}(t1), fet_maxs{bi}(t1), fet_mins{bi}(t2), fet_maxs{bi}(t2), fet_mins{bi}(t3), fet_maxs{bi}(t3)];                  %#ok<NASGU,AGROW>
            end
        end        
        
        for cell_i = 1:nCells
            s = styleStrs3D{cell_i};
            h_cell3D(cell_i) = plot3(0,0,0, s{:});
        end
        h_3D_ax = gca;
        box on; grid on;
        h_3D_xlab = xlabel(featureLabels{basis_id0}{1});
        h_3D_ylab = ylabel(featureLabels{basis_id0}{2}); 
        h_3D_zlab = zlabel(featureLabels{basis_id0}{3});
        
        
        view(3);
        legend(arrayfun(@(x) num2str(x), cellIds, 'un', 0));
        h_3D_tit = title(' ');        
    end
    
    
    if showEllipsoids0
        if show2D
            for p_i = 1:nSubplotsPerPage
                for c_i = 1:nCells
                    h_elps2D(p_i, c_i) = plot(h_2D_ax(p_i), 0,0, '-', 'color', cellColors(c_i,:)); %#ok<AGROW>
                end
            end            
        end
        if show3D
            for c_i = 1:nCells
                h_elps3D(c_i) = surf(h_3D_ax, [0;0], [0;0], zeros(2)); %#ok<AGROW>
            end
        end        
        
    end
        


    
    if showPairDists0
        basis_id1 = 3; % channel-whitened basis
        clust_overlaps = zeros(nCells, nCells);
        clust_distances = zeros(nCells, nCells);                       
        for cell_i = 1:nCells
            for cell_j = 1:cell_i
                M1 = clustMs{basis_id1,cell_i};   M2 = clustMs{basis_id1,cell_j};
                C1 = clustCovs{basis_id1,cell_i}; C2 = clustCovs{basis_id1,cell_j};                
                clust_overlaps(cell_i,cell_j) = quadProdGaussians(M1, C1, M2, C2)/quadProdGaussians(M1, C1, M1, C2);
                clust_distances(cell_i,cell_j) = norm(M1-M2);
            end
        end
        clust_overlaps(clust_overlaps == 0) = 1e-45;
        clust_overlaps = -log10(clust_overlaps);
        
               

        dendrogramOptions = {'waveforms_CC', 'waveforms_ED', 'gaussian_overlap', 'gaussian_dist'};
%         dendrogramOptions2 = {'waveforms_CC', 'waveforms_ED', 'gaussian_overlap', 'gaussian_dist'};
        cellIds_str = arrayfun(@num2str, cellIds, 'un', 0);
        % pre-compute all dendrograms
        for di = 1:length(dendrogramOptions)
            
            switch dendrogramOptions{di}
                case 'waveforms_CC', 
                    pdists_i = pdist(meanWaveformsCCW', 'correlation');  % == 1-corr(x,y)
                case 'waveforms_ED', 
                    pdists_i = pdist(meanWaveformsCCW', 'euclidean');  
                case 'gaussian_overlap', 
                    pdists_i = squareform (clust_overlaps); 
                case 'gaussian_dist'
                    pdists_i = squareform (clust_distances);
            end
            tree_i = linkage(pdists_i,'weighted'); 
            dendrograms.(dendrogramOptions{di}) = calcDendrogram(tree_i, 0, 'colorthreshold','default', 'label', cellIds_str);                        
            
        end
        
        
        % Initialize the figure with dendrograms & waveforms.
        figure(fig_pdist); set(fig_pdist, 'Name', 'Pairwise clust_distances', 'NumberTitle', 'off'); clf;         


        % Dendograms with controls at the top showing which for measures to display the dendrograms  
        nDendrogramsToDisplay = 2;
        dend_ids = [1, 3];

        h_dend_ax = zeros(1,nDendrogramsToDisplay);
        h_dendrograms = cell(1,nDendrogramsToDisplay);
        for dend_i = 1:nDendrogramsToDisplay        
            % the popup menu to choose which measure
            dendSelect_pos = [(4*dend_i-3)/(4*nDendrogramsToDisplay), .93,  1/(2*nDendrogramsToDisplay), .05];
            uicontrol('style', 'popupmenu', 'string', dendrogramOptions, 'parent', fig_pdist, 'value', dend_ids(dend_i), ...
                      'userdata', dend_i, 'units', 'normalized', 'position', dendSelect_pos, 'callback', @selectDendrogramType_callback);                                     

            % the actual plot with the dendrogram
            h_dend_ax    (dend_i) = mySubplot(2,nDendrogramsToDisplay,1,dend_i);
            h_dendrograms{dend_i} = plotDendrogram([], dendrograms.(dendrogramOptions{dend_ids(dend_i)}) );
        end
        3;
        
        % Plot the axis for showing the cell-waveforms        
        h_cellWvfms_ax = mySubplot(2,10,2,[1 7], [], [0 .1 0 .1]); hold on; 
        wvfm_basis_id = 2;
        wvfm_dt = diff(t_ms(1:2));
        wvfm_len = wvfm_dt*length(t_ms);
        extWvfmTicks = bsxfun(@plus, t_ms(:), wvfm_len*[0:3]);
        extWvfmTicks = insertNansAt(extWvfmTicks(:), [1:3]*nT+1);        
        meanWaveforms_plot_indiv{1} = insertNansAt(meanWaveforms_raw,    [1:3]*nT+1);        
        meanWaveforms_plot_indiv{2} = insertNansAt(meanWaveformsCCW, [1:3]*nT+1);        
        meanWaveforms_plot_grouped = meanWaveforms_plot_indiv{wvfm_basis_id};
        nT_plot = length(extWvfmTicks);

        h_meanWvfms = zeros(1,nCells);        
        for cell_i = 1:length(cellIds)
            h_meanWvfms(cell_i) = plot(h_cellWvfms_ax, extWvfmTicks, meanWaveforms_plot_indiv{wvfm_basis_id}(:,cell_i), '.-', 'color', cellColors_dens(cell_i,:));
        end 
        xlim([t_ms(1), t_ms(1)+ wvfm_len*4]);
%         waveformLims = [min(meanWaveforms_raw, [], 1); max(meanWaveforms_raw, [], 1)];
        for i = 1:3
            h_wvfmSepLine(i) = drawVerticalLine(t_ms(1) - wvfm_dt/2 + wvfm_len*i); 
        end

        
        % Popup-menu for choosing regular waveforms or CCW waveforms
        wvfmOptions = {'raw', 'CCW'};
        uicontrol('style', 'popupmenu', 'string', wvfmOptions, 'parent', fig_pdist, 'value', 2, ...
                  'units', 'normalized', 'position', [.65, .3, .05, .05], 'callback', @selectWaveformType_callback);                                     
        
%         meanWaveformsCCW
        % Axes for the current cell pair cross correlogram.
        h_cellPairCCG_ax = mySubplot(2,10,2,[8 10], [], [0 .1 0 .1]); hold on; 
        h_cellPairCCG    = bar(0, 0, 1, 'edgecolor', 'none', 'facecolor', 'b');

        h_cellCCG_range = uicontrol('style', 'popupmenu', 'string', isi_rangeOptions, 'parent', fig_pdist, 'value', 1, ...
                                    'tag', 'range_ms', 'units', 'normalized', 'position', [.65 .20, .065, .05], 'callback', @updateWvfmCellsCCG);
        h_cellCCG_nbins = uicontrol('style', 'popupmenu', 'string', isi_nbinOptions, 'parent', fig_pdist, 'value', 1, ...
                                    'tag', 'nbins', 'units', 'normalized', 'position', [.65 .15, .065, .05], 'callback', @updateWvfmCellsCCG);


        % Plot the checkboxes for showing/hiding each cell's waveform        
        chkBox_pos = @(i) [ (i) / (nCells+3), .005, 1/(nCells+3), .04];
        h_cellWvfmChkbox = zeros(1,nCells);        
        positionCellWvfmCheckboxes;                
        
        % Display the buttons for merging & splitting cells
        uicontrol('style', 'pushbutton', 'string', 'MERGE', ...
                  'units', 'normalized', 'position', chkBox_pos(nCells+1), 'parent', fig_pdist, 'callback', @mergeButton_callBack);
        uicontrol('style', 'pushbutton', 'string', 'SPLIT', ...
                  'units', 'normalized', 'position', chkBox_pos(nCells+2), 'parent', fig_pdist, 'callback', @splitButton_callBack);
        
        updateCellWaveformsVisibility;
        
    end


    function selectWaveformType_callback(src, evnt)
        id = get(src, 'value');
        wvfm_basis_id = id;
        updateCellWaveforms;
    end
    
    function selectDendrogramType_callback(src, evnt)
        id = get(src, 'userdata');
        dend_id = get(src, 'value');        
        plotDendrogram(h_dendrograms{id}, dendrograms.(dendrogramOptions{dend_id}) );
        
    end
            
    function showCellWaveform_callBack(src, evnt)
        button_cellIds = get(src, 'userdata');
        grp_idx = find( cellfun(@(x) x(1) == button_cellIds(1), cellIds_groups) );
        onOff = get(src, 'value');
        set(h_meanWvfms(grp_idx), 'visible', iff(onOff, 'on', 'off'));
        
        updateCellWaveformsVisibility;
    end


    function mergeButton_callBack(src, evnt)        
        idx_selected = find( cell2mat( get(h_cellWvfmChkbox, 'value')) );                
        if length(idx_selected) > 1
            set(h_cellWvfmChkbox(idx_selected(2:end)), 'value',0)
            mergeCellGroups(idx_selected);
            % un-select all except #1;            
        end
            
    end

    function mergeCellGroups(cellGroupIds_toMerge)
        id1 = cellGroupIds_toMerge(1);
        allCellIdsInNewGroup = [cellIds_groups{cellGroupIds_toMerge}];
        cellIds_groups{id1} = allCellIdsInNewGroup;
        cellIds_groups(cellGroupIds_toMerge(2:end)) = [];
        
        positionCellWvfmCheckboxes;        
        updateCellWaveforms;
    end

    function splitButton_callBack(src, evnt)
        idx_selected = find( cell2mat( get(h_cellWvfmChkbox, 'value')) );                
        if ~isempty(idx_selected)
            idx_grpsWithMultipleCells = idx_selected(cellfun(@length, cellIds_groups(idx_selected)) > 1);
            if ~isempty(idx_grpsWithMultipleCells)
                cellIds_toBeSelectedAfter = [cellIds_groups{idx_selected}];
                
                splitCellGroups(idx_grpsWithMultipleCells);            
                
                idx_chkBoxToBeSelected = find(cellfun(@(cids) any(cellIds_toBeSelectedAfter == cids(1)),  cellIds_groups));
                set(h_cellWvfmChkbox, 'value',0)
                set(h_cellWvfmChkbox(idx_chkBoxToBeSelected), 'value',1)
                updateCellWaveformsVisibility;
            end
        end
         
    end


    function updateWvfmCellsCCG(src, evnt)
        wvfmsOn = logical ( cell2mat( get(h_cellWvfmChkbox, 'value')) );                
        
        if nnz(wvfmsOn) == 2
            [ccg_cellIdx1, ccg_cellIdx2] = dealV(find(wvfmsOn));
                                
            rng_idx = get(h_cellCCG_range, 'value');
            nbin_idx = get(h_cellCCG_nbins , 'value');        
            
            range_ms = isi_allRanges_ms(rng_idx);
            nbins    = isi_allNbins(nbin_idx);
            ccgBinEdges = linspace(-range_ms, range_ms, nbins+1);
            ccgBinCent = binEdge2cent(ccgBinEdges);

            wvfmCCG_bin = masterBin2Bin(allCCGs_masterBin(ccg_cellIdx2, ccg_cellIdx1), ccg_masterBinIdxs{rng_idx, nbin_idx});                            
            set(h_cellPairCCG, 'xdata', ccgBinCent, 'ydata', wvfmCCG_bin{1}, 'visible', 'on');            
            
        else
            set(h_cellPairCCG, 'visible', 'off')
        end
        
    end

    function updateCellWaveformsVisibility
        wvfmsOn = logical ( cell2mat( get(h_cellWvfmChkbox, 'value')) );                
        set(h_meanWvfms( wvfmsOn), 'visible', 'on')
        set(h_meanWvfms(~wvfmsOn), 'visible', 'off')        
        for cell_grp_i = find(wvfmsOn(:)')            
            cellId1 = find(cellIds == cellIds_groups{cell_grp_i}(1), 1);
            set(h_meanWvfms(cell_grp_i), 'ydata', meanWaveforms_plot_grouped(:,cell_grp_i), ...
                'color', cellColors_dens(cellId1,: ) );
        end
        
        updateWvfmCellsCCG;
        
        
        allY = get(h_meanWvfms( wvfmsOn), 'ydata');
        if ~isempty(allY)
            if iscell(allY)
                allY = [allY{:}];
            end
            % update y-limits and vertical lines        
            lo = min(allY);
            hi = max(allY);
            if ~isempty(lo)
                set(h_wvfmSepLine, 'ydata', [lo, hi]);       
                set(h_cellWvfms_ax, 'ylim', [lo, hi]);
            end
        end
    end

    function splitCellGroups(cellGroupIds_toSplit)
        recoveredCellIds = num2cell([cellIds_groups{cellGroupIds_toSplit}]);        
        cellIds_groups(cellGroupIds_toSplit) = [];
        cellIds_groups = [cellIds_groups, recoveredCellIds];
        ascend_idx = ord( cellfun(@(x) x(1), cellIds_groups), 'ascend');
        cellIds_groups = cellIds_groups(ascend_idx);        
        
        positionCellWvfmCheckboxes;                
        updateCellWaveforms;
%         updateDendrogram;
    end

    
    function updateCellWaveforms

        meanWaveforms_plot_grouped = zeros(nT_plot, length(cellIds_groups));
        nSpksPerCell = cellfun(@length, cellSpikeIdx_orig);
        for cell_grp_i = 1:length(cellIds_groups)
            grpCellIds = arrayfun(@(cid) find(cellIds == cid), cellIds_groups{cell_grp_i});
            curWvfms = meanWaveforms_plot_indiv{wvfm_basis_id}(:,grpCellIds);
            if length(grpCellIds) > 1
                wgts = nSpksPerCell(grpCellIds);
                wgts = wgts/sum(wgts);                
                meanWaveforms_plot_grouped(:,cell_grp_i) = sum(bsxfun(@times, curWvfms, wgts(:)'), 2) ;
            else
                meanWaveforms_plot_grouped(:,cell_grp_i) = curWvfms;                
            end
        end
        
%         waveformLims = [min(meanWaveforms_plot_grouped, [], 1); max(meanWaveforms_plot_grouped, [], 1)];
        updateCellWaveformsVisibility;
%         updateDendrogram;
        
    end


    function updateDendrogram
        
        waveformPDists_v = waveformPDists;
%         waveformPDists_v = squareform(waveformPDists);
        wvfmPDistLinkage = linkage(waveformPDists_v,'weighted'); 
        cellIds_str = arrayfun(@num2str, cellIds, 'un', 0);
        wvfm_dendrogram = calcDendrogram(wvfmPDistLinkage, 'colorthreshold','default', 'label', cellIds_str);
%         [h_dendrogram, cols] = plotDendrogram([], wvfm_dendrogram);
        h_dendrogram = plotDendrogram(h_dendrogram, wvfm_dendrogram);
        
    end

    function positionCellWvfmCheckboxes
        nCellGroups = length(cellIds_groups);
%         if nGrps < length(h_cellWvfmChkbox)
%             set(h_cellWvfmChkbox(nGrps+1:end) = 0;
%         end        
        if isempty(h_cellWvfmChkbox)
            h_cellWvfmChkbox = zeros(1,nCellGroups);
        end
        nPerGrp = cellfun(@length, cellIds_groups);
        curChkBox_pos = @(i) [ ((sum(nPerGrp(1:(i-1))))+1) / (sum(nPerGrp)+3), .005, nPerGrp(i)/(sum(nPerGrp)+3), .04];
%         curChkBox_pos = @(i) [ (i) / (nCellGroups+3), .005, 1/(nCellGroups+3), .04];
        for cell_grp_i = 1:nCellGroups                        
            cellList_str = cellstr2csslist(arrayfun(@num2str, cellIds_groups{cell_grp_i}, 'un', 0));    
            cellId1 = find(cellIds == cellIds_groups{cell_grp_i}(1));
            if h_cellWvfmChkbox(cell_grp_i) == 0        
                h_cellWvfmChkbox(cell_grp_i) = uicontrol('style', 'checkbox',  'parent', fig_pdist, ... 
                    'units', 'normalized', 'position', curChkBox_pos(cell_grp_i), 'fontsize', 10, 'fontweight', 'bold', ...
                    'string', cellList_str, 'ForegroundColor', cellColors_dens(cellId1,:), 'userData', cellIds_groups{cell_grp_i}, ...
                    'callback', @showCellWaveform_callBack);
            else                
                pos = curChkBox_pos(cell_grp_i);
                set( h_cellWvfmChkbox(cell_grp_i), 'position', pos, ...
                    'string', cellList_str, 'ForegroundColor', cellColors_dens(cellId1,:), 'userData', cellIds_groups{cell_grp_i}, 'visible', 'on' );                
            end            
        end
        set(h_cellWvfmChkbox(nCellGroups+1:end), 'visible', 'off');
        
    end



    
    
    
    
    function cArr = CC(cArr)
        cArr = cellfun(@(x) x(:), cArr, 'un', 0);
        cArr = cat(1, cArr{:});
    end        

    
    

    function updateAllFigures(showDims, showEllipseRefrac, ellipse_nStd, featureSetOption, page, x_tetName, y_tetName, z_tetName, ...
         showTicks, showLabels, viewOption, densityTypeOption, densityScaleOption, ...
         showISIs, isi_pairsOption, isi_range_str, isi_nbins_str, refrPeriod_ms, isi_crossType, ...
         r2D_cellId )
     
        [show2D, show2Dr, show3D] = dealV(showDims);
        [showEllipsoids, showRefracSpikes] = dealV(showEllipseRefrac);
        

        basis_id = find(strcmp(featureSetOption, featureSetOptions));
        
%         showOvlps = num2cell(showOvlps);
%         [showPairDists, showISIs, showEllipsoids] = showOvlps{:};        
                
        basisChanged = ~strcmp(featureSetOption_prev, featureSetOption);
        showCellsChanged = any(showCells_prev ~= showCells);
        pageChanged = page_prev ~= page;
        
        viewChanged = ~strcmp(viewOption_prev, viewOption);
        showTicksChanged = (showTicks_prev ~= showTicks); 
        showLabelsChanged = (showLabels_prev ~= showLabels); 
        showRefracSpikesChanged = (showRefracSpikes_prev ~= showRefracSpikes);        
                
        densityTypeChanged = ~strcmp(densityTypeOption_prev, densityTypeOption);
        densityScaleChanged = ~strcmp(densityScaleOption_prev, densityScaleOption);
%         pairDistTypeOptionChanged = ~strcmp(pairDistTypeOption_prev, pairDistTypeOption);
        activeCellsIdx = find(showCells);
        isi_pairsChanged = ~strcmp(isi_pairsOption_prev, isi_pairsOption);
        isi_range_ms = sscanf(isi_range_str, '%d ms');
        isi_nbins = str2double(isi_nbins_str);        
        isiParamsChanged = any(isiHistParams_prev ~= [isi_range_ms, isi_nbins]);
        isiCrossTypesChanged = ~strcmp(isi_crossType_prev, isi_crossType);
        
       
        combos2D_here =  combos2D{basis_id};
        
        if show2D            

            if nPairsPerBasis(basis_id) <= nSubplotsPerPage
                cOff = 0; %co-ordinate offset
            else
                cOff = (page-1)*nSubplotsPerPage;
                if cOff >= nPairsPerBasis(basis_id)
                    cOff = floor((nPairsPerBasis(basis_id)-1)/nSubplotsPerPage)*nSubplotsPerPage; % in case page # is too high for this basis.                
                end
%                 cOff = min(cOff,
%                 nPairsPerBasis(basis_id)-nSubplotsPerPage+1); 
            end
            plotsAvailable = combos2D_here(cOff + [1:nSubplotsPerPage],1);
            activePlotsIdx = find(plotsAvailable(:)');
            inactivePlotsIdx = find(~plotsAvailable(:)');
                            
            % 1. Adjust axes limits.
            if basisChanged || pageChanged
                for plot_j = activePlotsIdx
                    curXlim = get(h_2D_ax(plot_j), 'xlim'); newXlim = axesLims2D{basis_id}(plot_j+cOff,1:2);
                    curYlim = get(h_2D_ax(plot_j), 'ylim'); newYlim = axesLims2D{basis_id}(plot_j+cOff,3:4);
                    if any([curXlim, curYlim] ~= [newXlim, newYlim])
                        set(h_2D_ax(plot_j), 'xlim', newXlim, 'ylim', newYlim)
                    end
                end
            end
            
            
            if strcmp(viewOption, 'individual spikes');
                if viewChanged || pageChanged || basisChanged
                    set(h_2D_im, 'visible', 'off');
                    set(h_2D_spks(:, activeCellsIdx), 'visible', 'on');
                end
                set(h_2D_refr(:, activeCellsIdx), 'visible', iff(showRefracSpikes, 'on', 'off'));
                set(h_2D_refr(:, ~showCells), 'visible', 'off');
                
                if basisChanged || pageChanged 
                    for plot_j = activePlotsIdx
                        idx = combos2D_here(plot_j+cOff,:);
                        assert(idx(1) > 0);
                        for c_j = activeCellsIdx
                            cell_x = features{basis_id}(cellSpikeIdx{c_j},idx(1)); cell_y = features{basis_id}(cellSpikeIdx{c_j},idx(2));
                            set(h_2D_spks(plot_j, c_j), 'xdata', cell_x, 'ydata', cell_y);                            
                            if showRefracSpikes
%                                 refrIdx = allCell_refrSpkIdxs{c_j};                           
%                                 [refr_x, refr_y] = linesFromAtoB(cell_x(refrIdx), cell_y(refrIdx), cell_x(refrIdx+1), cell_y(refrIdx+1));
                                refr_x = [cell_x(allCell_refrSpkIdxs{c_j}), cell_x(allCell_refrSpkIdxs{c_j}+1), nan(size(allCell_refrSpkIdxs{c_j}))]';
                                refr_y = [cell_y(allCell_refrSpkIdxs{c_j}), cell_y(allCell_refrSpkIdxs{c_j}+1), nan(size(allCell_refrSpkIdxs{c_j}))]';
                                set(h_2D_refr(plot_j, c_j), 'xdata', refr_x(:), 'ydata', refr_y(:) );
                            end                                                            
                            
                        end 

                    end                    
                end
                if basisChanged || pageChanged || showCellsChanged            
                    set(h_2D_spks(activePlotsIdx, showCells), 'visible', 'on');
                    set(h_2D_spks(activePlotsIdx, ~showCells), 'visible', 'off');                    
                    set(h_2D_spks(inactivePlotsIdx, :), 'visible', 'off');
                end

            elseif any(strcmp(viewOption, {'density', 'cluster density'}))
                if viewChanged || pageChanged || basisChanged
                    set(h_2D_im(activePlotsIdx), 'visible', 'on');                
                    set(h_2D_spks, 'visible', 'off');
                end
                                        
                if densityTypeChanged
                    switch densityTypeOption
                        case 'spike density',  colormap(h_2D_ax(1), cmap_spikes);
                        case 'cluster density', colormap(h_2D_ax(1), cmap_clusters);
                    end
                end

                                
                curClims = [nan, nan];
                
                if viewChanged || basisChanged || pageChanged || densityTypeChanged || showCellsChanged || densityScaleChanged
                
                    for plot_j = activePlotsIdx
                        ax = axesLims2D{basis_id}(plot_j+cOff,:);     
                        ax_tck = axesTickLims2D{basis_id}(plot_j+cOff,:);     

                        [cdata, curClims] = getSpikeDensityCData(spkDensities{basis_id}, plot_j+cOff, ...
                            showCells, densityTypeOption, densityScaleOption, curClims, cmapSclFactor_cell, cmapSclFactor_grp);

                        set(h_2D_im(plot_j), 'xdata', ax_tck(1:2), 'ydata', ax_tck(3:4), 'cdata', cdata' );                    
                        set(h_2D_ax(plot_j), 'xlim', ax(1:2), 'ylim', ax(3:4));
                    end

                    set(h_2D_im(inactivePlotsIdx), 'visible', 'off');
                end
                
                if any(isnan(curClims)),
                    curClims = [0 1];
                end                
                
                switch densityTypeOption
                    case 'spike density', 
                        set(h_dummy_im, 'cdata', [curClims]);
                        set([h_2D_ax h_dummy_ax], 'climmode', 'auto');
                        set(hColorBar, 'yTickMode', 'auto', 'yTickLabelMode', 'auto');
                    case 'cluster density', 
                        set(h_dummy_im, 'cdata', [0, nCells]);
                        set(h_dummy_ax, 'climMode', 'auto');
                        set(h_2D_ax, 'clim', [0 nCells]);
                        drawnow;
%                         set(hColorBar, 'yTick', [0:nCells-1]+5, 'yTickLabel', cellfun(@num2str, num2cell(0:nCells-1), 'un', 0));
                        skp = ceil(nCells/10);%iff(nCells<10, nCells, round( nCells/10 );
                        cellIdIdx = 1:skp:nCells;
                        set(hColorBar, 'yTick', [cellIdIdx]-.5, 'yTickLabel', arrayfun(@num2str, cellIds(cellIdIdx), 'un', 0));
                        3;
                end          
%                 drawnow;
%                 colormap(h_2D_ax(plot_j), cmap );
                4;
            else
                error('Unknown view option');
            end
            

            if showEllipsoids0
                if showEllipsoids
                    for plot_j = activePlotsIdx
                        idx = combos2D_here(plot_j,:);
                        for c_j = activeCellsIdx
                            [x_elps,y_elps] = ellipsoidFromCov(clustMs{basis_id,c_j}(idx), clustCovs{basis_id,c_j}(idx, idx), ellipse_nStd, nEllipsePoints2D);
                            set(h_elps2D(plot_j, c_j), 'xdata', x_elps, 'ydata', y_elps, 'visible', 'on', 'color', 'k', 'linestyle', '-')
                        end 
                    end
                    set(h_elps2D(:, [1, find(~showCells)]), 'visible', 'off') 
                else
                    set(h_elps2D, 'visible', 'off')
                end
            end

            
            if basisChanged || pageChanged || showTicksChanged
                axesTicksToBeShown = (plotsAvailable(:) & showTicks);                
                xticks = get(h_2D_ax, 'xtick'); if ~iscell(xticks), xticks = {xticks}; end;
                axesTicksShowing = ~cellfun(@isempty, xticks);
                axesTicksToTurnOn  = axesTicksToBeShown & ~axesTicksShowing;
                axesTicksToTurnOff = ~axesTicksToBeShown & axesTicksShowing;
                if any(axesTicksToTurnOn)
                    set(h_2D_ax(axesTicksToTurnOn), 'xtickmode', 'auto', 'ytickmode', 'auto');
                end
                if any(axesTicksToTurnOff)
                    set(h_2D_ax(axesTicksToTurnOff), 'xtick', [], 'ytick', []);
                end
            end
            
            if basisChanged || pageChanged || showLabelsChanged
                axesLabelsToBeShown = (plotsAvailable(:) & showLabels);           
                for plot_j = find(axesLabelsToBeShown')
                    idx = combos2D_here(plot_j+cOff,:);
                    set(h_2D_xlab(plot_j), 'string', featureLabels{basis_id}{idx(1)});
                    set(h_2D_ylab(plot_j), 'string', featureLabels{basis_id}{idx(2)});                    
                end
                if any(~axesLabelsToBeShown)
                    set([h_2D_xlab(~axesLabelsToBeShown), h_2D_ylab(~axesLabelsToBeShown)], 'string', '');
                end                
            end
                                    
        end
                
        if show2Dr            

            
        end
                                    
            
        
        
        if show3D
            x_tet_id = str2double(x_tetName(9));    xs = features{basis_id}(:,x_tet_id);
            y_tet_id = str2double(y_tetName(9));    ys = features{basis_id}(:,y_tet_id);
            z_tet_id = str2double(z_tetName(9));    zs = features{basis_id}(:,z_tet_id);
            axID_3D = [x_tet_id, y_tet_id, z_tet_id];
            axIDsChanged = any(axID_3D ~= axID_3D_prev);
            
            cellsOn = sort(activeCellsIdx);
            ovlpAvailable = (length(cellsOn) == 2) && exist('bestfeatures', 'var') && ~isempty(bestfeatures{cellsOn(2), cellsOn(1)});            
            
            for ci = activeCellsIdx
                set(h_cell3D(ci), 'xdata', xs(cellSpikeIdx{ci}), 'ydata', ys(cellSpikeIdx{ci}), 'zdata', zs(cellSpikeIdx{ci}), 'visible', 'on');
            end
            set(h_cell3D(~showCells), 'visible', 'off');
            cellIds_str = sprintf(' %d,', cellIds(showCells));
            tit_str = sprintf('Gid = %d; cells %s', Gid, cellIds_str(1:end-1));
            if ovlpAvailable
                ovlp_str = sprintf('Gauss overlap = %.2g', clust_overlaps(cellsOn(2), cellsOn(1)));
                set(h_3D_tit, 'string', {tit_str, ovlp_str})
            else
                set(h_3D_tit, 'string', tit_str)
            end
            
            xlims = [fet_mins{basis_id}(x_tet_id), fet_maxs{basis_id}(x_tet_id)];
            ylims = [fet_mins{basis_id}(y_tet_id), fet_maxs{basis_id}(y_tet_id)];
            zlims = [fet_mins{basis_id}(z_tet_id), fet_maxs{basis_id}(z_tet_id)];
            set(h_3D_ax, 'xlim', xlims, 'ylim', ylims, 'zlim', zlims );

            if showEllipsoids
    
                idx = [x_tet_id, y_tet_id, z_tet_id];
                for c_j = activeCellsIdx
                    X = [xs(cellSpikeIdx{c_j}), ys(cellSpikeIdx{c_j}), zs(cellSpikeIdx{c_j})];                    
                    Ms_1 = mean(X,1)';
                    Cov_1 = cov(X);
                    
                    Ms_2 = clustMs{basis_id,c_j}(idx);
                    Cov_2 = clustCovs{basis_id,c_j}(idx, idx);
                    
                    assert(all( abs(Ms_1-Ms_2) < 1e-5 ) );
                    assert(all( abs(Cov_1(:)-Cov_2(:)) < 1e-5 ) )                    
                    
%                     [x,y,z] = ellipsoidFromCov(clustMs{c_j}(idx), clustCovs{c_j}(idx, idx), ellipse_nStd, nEllipsePoints3D);
                    [x,y,z] = ellipsoidFromCov(Ms_1, Cov_1, ellipse_nStd, nEllipsePoints3D);
                    col = max(cellColors(c_j,:)-.1,0);
                    set(h_elps3D(c_j), 'xdata', x, 'ydata', y, 'zdata', z, 'visible', 'on', ...
                         'facecolor', col, 'facealpha', .1, 'edgecolor', col, 'edgealpha', .2 );
                end
                set(h_elps3D([1, find(~showCells)]), 'visible', 'off') 
                
%                 set(h_3D_xlab, 'string', ' ');
            else
                set(h_elps3D, 'visible', 'off')
            end
            
            if basisChanged || showLabelsChanged || axIDsChanged
                if showLabels
                    set(h_3D_xlab, 'string', sprintf('%s %d', featureLabels{basis_id}{x_tet_id}));
                    set(h_3D_ylab, 'string', sprintf('%s %d', featureLabels{basis_id}{y_tet_id}));                    
                    set(h_3D_zlab, 'string', sprintf('%s %d', featureLabels{basis_id}{z_tet_id}));                    
                else
                    set([h_3D_xlab, h_3D_ylab, h_3D_zlab], 'string', '');
                end                
            end
            
        end
            3;
            

        
        if showISIs && (  (showCellsChanged || isiParamsChanged || isi_pairsChanged || basisChanged) ...
...                            || (strcmp(isi_pairsOption, '2Dr auto') && r2DCellChanged ) ...
                           || (~isempty(strfind(isi_pairsOption, 'cross')) && isiCrossTypesChanged ) )

           % isi_pairsChanged            
            isi_range_idx = find(isi_range_ms == isi_allRanges_ms, 1);
            isi_nbin_idx = find(isi_nbins == isi_allNbins, 1);

            onlyShowActiveCells = true;                    
                        
            if onlyShowActiveCells                
                if ~strcmp(isi_pairsOption, '2Dr auto')
                    isi_cellIdxs = activeCellsIdx;
                else
                    isi_cellIdxs = refr_cellIdx;
                end
%                 h_toHide = [h_isi_ax(~showCells, :);  h_isi_ax(:,~showCells)'; 
%                             h_isi_bar(~showCells, :); h_isi_bar(:,~showCells)'; ]; 
                h_toHide = [h_isi_ax(:); h_isi_bar(:); h_isi_txt(:); h_isi_refrLine(:)];
                set(h_toHide(h_toHide > 0), 'visible', 'off');
            else
                isi_cellIdxs = 1:nCells;
            end                
                
            isi_cellIds = cellIds(isi_cellIdxs);
            nAcCells = length(isi_cellIdxs);

            isiBinEdges = linspace(0, isi_range_ms, isi_nbins+1);
            isiBinCent = binEdge2cent(isiBinEdges);
            ccgBinEdges = linspace(-isi_range_ms, isi_range_ms, isi_nbins+1);
            ccgBinCent = binEdge2cent(ccgBinEdges);
            
            
            if isiParamsChanged || isiCrossTypesChanged     
                allISIs_bin = masterBin2Bin(allISIs_masterBin, isi_masterBinIdxs{isi_range_idx, isi_nbin_idx});                
                if strcmp(isi_crossType, 'cross correlogram')
                    allCCGs_bin = masterBin2Bin(allCCGs_masterBin, ccg_masterBinIdxs{isi_range_idx, isi_nbin_idx});                
                end
            end

            for cell_I_idx_idx = 1:length(isi_cellIdxs)
                cell_I_idx = isi_cellIdxs(cell_I_idx_idx);
                cell_I = cellIds(cell_I_idx);
                
                switch isi_pairsOption
                    case 'auto & cross',            isi_cell_Js_idxs = isi_cellIdxs(isi_cellIdxs <= cell_I_idx);                                                
                    case {'auto only', '2Dr auto'}, isi_cell_Js_idxs = isi_cellIdxs(isi_cellIdxs == cell_I_idx);
                    case 'cross only',              isi_cell_Js_idxs = isi_cellIdxs(isi_cellIdxs <  cell_I_idx);                        
                    otherwise, error('wrong option')
                end
                                
                for cell_J_idx_idx = 1:length(isi_cell_Js_idxs)
                    cell_J_idx = isi_cell_Js_idxs(cell_J_idx_idx);
                    cell_J = cellIds(cell_J_idx);              
                    
%                     J_idx = find(cell_J_idx == isi_cellIdxs, 1);
                    set(CC({h_isi_ax(cell_I_idx,cell_J_idx), h_isi_bar(cell_I_idx,cell_J_idx)}), 'visible', 'on');   
%                     set(h_isi_ax(cell_I_idx,cell_J_idx), 'color', 'w');
                    
                    switch isi_pairsOption
                        case 'auto & cross',               M = nAcCells; N = nAcCells; m = nAcCells-cell_I_idx_idx+1; n = cell_J_idx_idx;
                        case {'auto only', '2Dr auto'},    M = floor(sqrt(nAcCells)); N = ceil(nAcCells/M); m = cell_I_idx_idx; n = [];
                        case 'cross only',                 M = nAcCells-1; N = nAcCells-1; m = nAcCells-cell_I_idx_idx+1; n = cell_J_idx_idx;
                    end
                    mySubplot(M,N,m,n, h_isi_ax(cell_I_idx,cell_J_idx) );
                    
%                     h_isi_ax(cell_I_idx,cell_J_idx,:) = mySubplot(nAcCells,nAcCells, nAcCells-cell_I_idx+1, cell_J_idx, h_isi_ax(cell_I_idx,cell_J_idx,:) );
                    if (cell_I ~= cell_J) && strcmp(isi_crossType, 'cross correlogram')
                        xdata = ccgBinCent;
                        xlims = [ccgBinEdges(1), ccgBinEdges(end)];
                        ybardata = allCCGs_bin{cell_I_idx, cell_J_idx};
                    else
                        xdata = isiBinCent;
                        xlims = [isiBinEdges(1), isiBinEdges(end)];
                        ybardata = allISIs_bin{cell_I_idx, cell_J_idx};
                    end
                    set(h_isi_bar(cell_I_idx,cell_J_idx), 'xdata', xdata, 'ydata', ybardata);
                    
                    txt_col = iff(refracViol_n(cell_I_idx, cell_J_idx) > 10 || refracViol_pct(cell_I_idx, cell_J_idx) >.5, 'r', 'k');
                    txt_prefix = iff(cell_I == cell_J, sprintf('[%d]', cell_I), sprintf('[%d,%d]', cell_J, cell_I));

                    s1 = sprintf('%s<%.1fms;', txt_prefix, refrPeriod_ms);
                    s2 = sprintf('%d spks, (%.1f%%)', refracViol_n(cell_I_idx, cell_J_idx), refracViol_pct(cell_I_idx, cell_J_idx));
                    ylims = get(h_isi_ax(cell_I_idx,cell_J_idx), 'ylim');
                    if (cell_I ~= cell_J) && strcmp(isi_crossType, 'cross correlogram')
                        set(h_isi_txt(cell_I_idx, cell_J_idx), 'visible', 'off');
                    else
                        set(h_isi_txt(cell_I_idx, cell_J_idx), 'position', [isi_range_ms/2, ylims(2)], 'string', {s2}, 'color', txt_col, 'visible', 'on');
                    end
                    
                    if (cell_I == cell_J) && show2Dr
                        set(h_isi_refrLine(refr_cellIdx), 'xdata', [1,1]*cellsRefrPeriod_ms(refr_cellIdx), 'ydata', ylims, 'visible', 'on' );
                    end                    
                    set(h_isi_ax(cell_I_idx,cell_J_idx), 'xlim', xlims);
                    if (length(isi_cellIdxs) == 1) || (strcmp(isi_pairsOption, 'cross only') && (length(isi_cellIdxs) == 2))
                        set(h_isi_ax(cell_I_idx,cell_J_idx), 'xtickmode', 'auto', 'ytickmode', 'auto')
                        set(h_isi_tit(cell_I_idx,cell_J_idx), 'string', sprintf('ISI for cell %d', cellIds(cell_I_idx) ));
                    else
%                         set(h_isi_ax(cell_I_idx,cell_J_idx), 'xtickmode', 'auto', 'ytickmode', 'auto')
                        set(h_isi_ax(cell_I_idx,cell_J_idx), 'xtick', [], 'ytick', []);                        
                        set(h_isi_tit(cell_I_idx,cell_J_idx), 'string', '');                        
                    end
                    
                end
            end
            
            if ~isempty(strfind(isi_pairsOption, 'cross'))
                offset = isempty(strfind(isi_pairsOption, 'auto'));
                posLL = get(h_isi_ax(isi_cellIdxs(1+offset),isi_cellIdxs(1) ), 'position');
                posUR = get(h_isi_ax(isi_cellIdxs(end),isi_cellIdxs(end-offset) ), 'position');
                xticks = isi_cellIds(1:end-offset);
                yticks = isi_cellIds(1+offset:end);
                xTickLabels = arrayfun(@num2str, xticks, 'un', 0);
                yTickLabels = arrayfun(@num2str, yticks, 'un', 0);
                tickaxpos = [posLL(1:2), posUR(1:2)+posUR(3:4)-posLL(1:2)];                
                set(h_isi_tick_ax, 'position', tickaxpos, ...
                    'xlim', [.5 length(xticks)+.5], 'ylim', [.5 length(yticks)+.5], 'xtick', 1:length(xticks), 'ytick', 1:length(yticks), ...
                    'xticklabel', xTickLabels, 'yticklabel', yTickLabels, 'visible', 'on');
            else
                set(h_isi_tick_ax, 'visible', 'off')
            end
            
            
        end
                       
        [featureSetOption_prev, showTicks_prev, showLabels_prev, viewOption_prev, isiHistParams_prev,        densityTypeOption_prev, page_prev, axID_3D_prev, showRefracSpikes_prev, showCells_prev, densityScaleOption_prev, isi_pairsOption_prev, isi_crossType_prev] = deal(...
         featureSetOption,      showTicks,      showLabels,      viewOption,      [isi_range_ms, isi_nbins], densityTypeOption,      page,      axID_3D,      showRefracSpikes,      showCells,      densityScaleOption,      isi_pairsOption,      isi_crossType);        
    end


    if doAutoPruning
        [pruningAutoData, assignToMU] = autoClassifyAllClusters(1);
        saveCurrentEdits;
        if (Gid < 6000) 
            fn = getFileName('prunedClustersStats', Gid, 1);
            save(fn, 'pruningAutoData');
        end
        
        return;        
    end



    [featureSetOption_prev, viewOption_prev, densityTypeOption_prev, ...
        densityScaleOption_prev, isi_pairsOption_prev] = deal('');
    [showTicks_prev, showLabels_prev, showRefracSpikes_prev, isi_crossType_prev] = deal(nan);
    axID_3D_prev = [1 2 3];    
    axID_3D = axID_3D_prev;   
    showCells_prev = 0;
    isiHistParams_prev = [0 0];    
    page0 = 1;
    page_prev = 0;
    vHandles = [];
    
    cond_tf = @(tf0) iff(tf0, [false, true], false);
    
    tf = [false, true];
    % Arg #1 : showDim        
    tf2D  = cond_tf(show2D);
    tf2Dr = cond_tf(show2Dr);
    tf3D  = cond_tf(show3D);   
    showDimNames = {'2D', '2Dr', '3D'};
    showDimOptions = {tf2D, tf2Dr, tf3D};
    showDim0 = [show2D, show2Dr, show3D];    
    depVars_2D = {{}, {}};    
    
    if show2Dr
        depVars_2Dr = {{}, {'2Dr_cell'}};  % angles.
    else
        depVars_2Dr = {};
    end
    r2D_cellId0 = cellIds(2);
    r2D_auto0 = true;
    
    if show3D
        depVars_3D = {{}, {'x_coord', 'y_coord', 'z_coord'}};
        featurestr = {'tetrode 1', 'tetrode 2', 'tetrode 3', 'tetrode 4'};
        tet0s = [1 2 3];
    else
        depVars_3D = { };
        featurestr = {'tetrode 1'};
        tet0s = [1 1 1];
    end    
    showDimsDepVars = {depVars_2D, depVars_2Dr, depVars_3D};

        
    % Arg #2: show ellipsoids & refractory spikes
    tf_ellipsoids = cond_tf(showEllipsoids0);
    tf_refrac = [false, true];
    showEllipseRefrac = {tf_ellipsoids, tf_refrac};
    showEllipseRefrac0 = [false, false];
    showEllipseRefracLabels = {'ellipsoids', 'refractory spikes'};
    
    depVars_ellipsoids = iff(showEllipsoids0, {{}, {'nStdEllipse'}},  {{}});
    showEllipseRefracDepVars = { depVars_ellipsoids, {} };
        
       
    
    % Arg : showRefracSpikes
%     showRefracSpikes = [false, true];
    
    ellipse_nStdRange = [.5:.1:5];
    ellipse_nStd0 = 3;
    
    
    % Arg #2 : basis.
    featureSetOption0 = featureSetOptions{1};

        
    
    % Arg #7 : showPairDists

    depVars_isis = {{}, {'ISI_range', 'ISI_nbins', 'ISI_pairs', 'ISI_refrPer_ms', 'ISI_crossType'}};        
    
%     showOvlpDepVars = {depVars_pdist, depVars_isis};    
%     pairDistTypeOptions = {'distance', 'overlap', 'waveformPDists'};

    % Args #12+ ISI parameters
    isi_pairsOptions = {'auto only', 'auto & cross', 'cross only', '2Dr auto'};    
    if ~showISIs0
        isi_pairsOptions(1) = [];
    end
    isi_pairsOption0 = isi_pairsOptions{1};
        
    isi_allRanges_str = arrayfun(@(x) [num2str(x) ' ms'], isi_allRanges_ms, 'un', 0);
    isi_allNbins_str = arrayfun(@num2str, isi_allNbins, 'un', 0);
    isi_range0 = isi_allRanges_str{4};
    isi_nbins0 = isi_allNbins_str{3};
    isi_crossTypes = {'auto correlogram', 'cross correlogram'};
    isi_crossTypes0 = isi_crossTypes{1};    

%     varGroups = {   {'alsoShow', 'basis', 'page', 'x_coord', 'y_coord, z_coord', 'xy_ticks', 'xy_labels', 'spikeView', 'densityType', 'densityScale'}, ...
%                     {'showISIs', 'ISI_pairs', 'ISI_range', 'ISI_nbins', 'ISI_refrPer_ms', 'ISI_crossType'} ...                    
%                 };        
    if ~show2D
        return;
    end
    
    args = {{'showDims', showDimOptions, showDim0, showDimsDepVars, showDimNames}, ...
            {'alsoShow', showEllipseRefrac, showEllipseRefrac0, showEllipseRefracDepVars, showEllipseRefracLabels}, ...
            {'nStdEllipse', ellipse_nStdRange, ellipse_nStd0}, ...
            {'basis', featureSetOptions, featureSetOption0}, ...
            {'page', [1:nPagesMax], page0}, ...
            {'x_coord', featurestr, featurestr{tet0s(1)}}, ...
            {'y_coord', featurestr, featurestr{tet0s(2)}}, ...
            {'z_coord', featurestr, featurestr{tet0s(3)}}, ...
            {'xy_ticks', [false, true], showTicks0}, ...
            {'xy_labels', [false, true], showLabels0}, ...
            {'spikeView',  spikeViewOptions, spikeViewOption0, spikeViewDepVars}, ...
            {'densityType', densityTypeOptions, densityTypeOption0}, ...
            {'densityScale', densityScaleOptions, densityScaleOption0}, ...            
            {'showISIs', tf, showISIs0, depVars_isis}, ...
            {'ISI_pairs', isi_pairsOptions, isi_pairsOption0}, ...
            {'ISI_range', isi_allRanges_str, isi_range0}, ...
            {'ISI_nbins', isi_allNbins_str, isi_nbins0}, ...
            {'ISI_refrPer_ms', refrPeriod_ms_range, refrPeriod_ms0}, ...
            {'ISI_crossType', isi_crossTypes, isi_crossTypes0}, ...
            {'2Dr_cell', cellIds, r2D_cellId0},            
            }; 

        
        
        
        
        
    vHandles = manipulate(@updateAllFigures, args, 'Title', 'Control Panel', 'FigId', 101);

    
    
end



function val = iff(bool, valIfTrue, valIfFalse)
    if bool
        val = valIfTrue;
    else
        val = valIfFalse;
    end
end

function toggleValue(hnd)
    val = get(hnd, 'value');
    set(hnd, 'value', 1-val);
end
      


function meanWaveforms = getClusterMeanWaveforms(cellSpikeIdx, spkWaveforms)
    [nT, nChannels, nSpk] = size(spkWaveforms);     %#ok<NASGU>
    nCells = length(cellSpikeIdx);

    meanWaveforms = zeros(nT, nChannels, nCells);
    for cell_i = 1:nCells
        meanWaveforms(:, :, cell_i) = mean(spkWaveforms(:,:, cellSpikeIdx{cell_i}),3);
    end           
end


    

function [refracViol_n, refracViol_pct, isiColorIdxs] = countRefractoryPeriodViolations(...
    allISIs_masterbin, masterBinEdges, refrPeriod_ms, multiUnitISI_ms, nColrs, cellSpkCounts)
    ncells = length(allISIs_masterbin);
    refracViol_n = zeros(ncells);    
    refracViol_pct = zeros(ncells);    

    if nargout == 3
        allU_beforeCutoff_pct = nnz( multiUnitISI_ms < refrPeriod_ms ) / length(multiUnitISI_ms) * 100;
        isiColorIdxs = zeros(ncells);
    end                    

    relevantBinIdxs = find(masterBinEdges(2:end) < refrPeriod_ms);
    
    for cell_i = 1:ncells
        for cell_j = 1:cell_i
            refracViol_n    (cell_i, cell_j) = sum(allISIs_masterbin{cell_i, cell_j}(relevantBinIdxs));
            n = iff(cell_i == cell_j, cellSpkCounts(cell_i), cellSpkCounts(cell_i)+cellSpkCounts(cell_j));
            refracViol_pct(cell_i, cell_j) = refracViol_n(cell_i, cell_j) / n * 100 ;
            
            if nargout == 3
                isiColorIdx =  round ( refracViol_pct(cell_i, cell_j) / allU_beforeCutoff_pct *(nColrs-1) )+1;
                isiColorIdxs(cell_i, cell_j) = min(isiColorIdx, nColrs);
            end                          
        end
    end        
end



function uX = uniqueC(X)
    Y = X;
    idx_c = cellfun(@iscell, Y);
    for i = find(idx_c)
        Y{i} = cellstr2csslist(Y{i});
    end
    [uy, m] = unique(Y, 'first');    
    uX = X(sort(m));
    
end




%{            
criteria for multiunit:

	0. [marked as MU from # refr violations]
	1. MD: (a) mode of MDs is close to threshold
           (b) std dev is > std(MD_all)
    2. mean spike waveform: multiple negative peaks
		-> second peak (after one at 0) is > 50% of 0 peak	
	3. autocut: remaining < 66%

%}


%         if useICrefrPeriod
%             IC_cand_idx = find(IC_stats.IC_candidates == cellIds(cellIdx));
% %             EC_clust_idx = find(IC_stats.EC_clustIds == cellIds(cellIdx), 1);
% %             fracOfClusterIC = IC_stats.fracOfCluster_IC(EC_clust_idx);
%             
%             if useICrefrPeriod_useIndivClustRefrPeriod && ~isempty(IC_cand_idx) %% && (fracOfClusterIC > 0.4)                                            
%                 true_refrPeriod_ms = IC_stats.IC_refrPeriod_EC_clust_ms(IC_cand_idx);
%             else
%                 true_refrPeriod_ms = IC_stats.IC_refrPeriod_all_EC_ms;
%             end
%             if isnan(true_refrPeriod_ms)
%                 true_refrPeriod_ms = [];
%             end
%             
%         end


%         clustDataOpt = struct;
%         if strcmp(pruningMode, 'pruneMergePruning')        
%             clustDataOpt = struct('idealInitPruning', useICrefrPeriod);
%         end
%         ClustData = getClusterData(Gid, cluster_data_idx, [], [], clustDataOpt);