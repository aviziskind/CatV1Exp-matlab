 function tmp_cmap(Gid, skipIfDoneFlag)

    if nargin < 1
        return;
    end
    
    skipIfDone = exist('skipIfDoneFlag', 'var') && ~isempty(skipIfDoneFlag);
    
    clusterPruningFile = getFileName('clusterPrunings', Gid);    
    spikeMergesFile = getFileName('clusterMerges', Gid, 1);

    if skipIfDone && exist(spikeMergesFile, 'file')
        fprintf('Group %d clusters have already been merged. Exiting ...\n', Gid);           
%         mergeClustersIntoCells(Gid);
        return;        
    end

    
    
    wvfmWidth = 1;
%     useSavedData = 1;

%     kmean_dist = 'sqEuclidean';
%     nKmeansReps = 200;
%     kmean_dist = 'correlation';
%     nKmeansReps = 50;

    leaveOutMultiUnit = 1;       

%     distMeasureType = 'cluster';
%     distMeasureType = 'spike';

%     curDistMeasure = 'cc';        
%     spreadMeasure0 = 'GLF';
    
    showWvfm0 = 1;
    showOTC0 = 0;
    showSPF0 = 0;
    showSpikeDensities = 1;
    showRasters = 1;
        doColorBar = 1;

    curNcellMode = 'scale';
    % waveforms: 
    
    meanSpikeWvfmFile = getFileName('mwvfm_clustP', Gid);
%     [spikeWaveforms_raw, t_ms] = getSpikeWaveforms(Gid, [], [], 'raw');
%     spikeWaveforms_ccw = getSpikeWaveforms(Gid, [], [], 'ccw');
%     spikeWaveforms_raw = single(spikeWaveforms_raw);
    allMeanWvfms = load(meanSpikeWvfmFile);
    t_ms = allMeanWvfms.t_ms;
    nT = length(t_ms);
    
    if leaveOutMultiUnit
        idx_keep = [allMeanWvfms.meanWaveforms.cellId] ~= 0;
    else
        idx_keep = 1:length(allMeanWvfms.meanWaveforms);
    end
    
    clustMeanWaveforms_raw = [allMeanWvfms.meanWaveforms(idx_keep).wvfm_raw]';
    clustMeanWaveforms_scl = [allMeanWvfms.meanWaveforms(idx_keep).wvfm_scl]';
    clustMeanWaveforms_ccw = [allMeanWvfms.meanWaveforms(idx_keep).wvfm_ccw]';
    
    clustMeanWaveforms_raw_plot = insertNansAt(clustMeanWaveforms_raw', [1:3]*nT+1)';
    clustMeanWaveforms_scl_plot = insertNansAt(clustMeanWaveforms_scl', [1:3]*nT+1)';
    clustMeanWaveforms_ccw_plot = insertNansAt(clustMeanWaveforms_ccw', [1:3]*nT+1)';
    
    idxNegAmps = find(t_ms == 0,1) + [0:3]*(nT+1);
    negAmps = clustMeanWaveforms_scl_plot(:,idxNegAmps);    


    requirePrunedClusters = 1;
    controlFig = 201;
    wvfmsFig = 202;
    isisFig = 203;
    ccgsFig = 204;
    rastersFig = 205;
    spkDensFig = 211;
    mergeControlFig = 411;
    showClustsFig = 412;
    

    
    nSpkTh = [];

    
    % get spike times (ms)
    siteInfo = siteDataFor(Gid);
    samplingRateHz = siteInfo.dataFileInfo.samplingRateHz;
    expDuration_sec = siteInfo.dataFileInfo.duration_sec;
    ticksPerMs = samplingRateHz/1000;
    spikePropertiesFile = getFileName('properties', Gid);
    S_spkProp = load(spikePropertiesFile);
    spikeTimes_ms = double(S_spkProp.position)/ticksPerMs;    
    
    try 
        curCSort = 'clustersPruned';
        [clustIds, clustSpkIdxs] = getCellSorting(Gid, curCSort);        
    catch MErr
        if requirePrunedClusters
            rethrow(MErr);
        end
        curCSort = 'clusters';
        [clustIds, clustSpkIdxs] = getCellSorting(Gid, curCSort);        
        warning('Using un-pruned clusters');
    end
    
    if leaveOutMultiUnit
        idx_nonMU = clustIds ~= 0;
        clustIds = clustIds(idx_nonMU);
        clustSpkIdxs = clustSpkIdxs(idx_nonMU);
    end

%         clustNSpikes = clustNSpikes(idx_clustsUse);
%         clustFiringRates = clustFiringRates(idx_clustsUse);    
    if ~isempty(nSpkTh)        
        clustNSpikes_tmp = cellfun(@length, clustSpkIdxs);
        clustIds = clustIds(clustNSpikes_tmp > nSpkTh);
        clustSpkIdxs = clustSpkIdxs(clustNSpikes_tmp> nSpkTh);
    end
    
    clustNSpikes = cellfun(@length, clustSpkIdxs);
    clustFiringRates = clustNSpikes / expDuration_sec;    
        
    clustIds_orig = clustIds;
    clustSpkIdxs_orig = clustSpkIdxs;
    clustNSpikes_orig = clustNSpikes;    
    
    [channelMeans, channelCov] = getChannelMeansAndCovariance(Gid);            
    channelScaling = 1./sqrt(diag(channelCov));    
    meanChannelStd = mean(1./channelScaling);

    clustNSpikes = cellfun(@length, clustSpkIdxs);    
    
    nClust = length(clustSpkIdxs);
    if nClust <= 1
        switch nClust
            case 0, fprintf('Group %d has no usable clusters.\n', Gid);   
                saveStruct1 = struct('curNCells', 0);                 %#ok<NASGU>
            case 1, fprintf('Group %d has only 1 cluster -- no merging is required.\n', Gid);        
                saveStruct1 = struct('absClusterGrouping', {{1}}, ...
                    'absUseClusters', true, 'clustIds', clustIds_orig, ...
                    'useClusters', true, 'curCellGrouping', {{1}}, 'curNCells', 1 );  %#ok<NASGU>
        end
        
        save(spikeMergesFile, '-struct', 'saveStruct1');        
        mergeClustersIntoCells(Gid);
        return;        
        
    end    
    clustsAvailable = true(1,nClust);
    clustIdxs = find(clustsAvailable);
    inactiveClustIdxs = [];
    
    nClust_orig = nClust;
    SavedData = getClusterData(Gid, 2);

    nSpikesUsed = sum(clustNSpikes);    
    nSpikesTot = length(spikeTimes_ms);
    Gid_str = sprintf('Gid = %d. Nspk = %d (%d total)', Gid, nSpikesUsed, nSpikesTot);
        
    
    nPointsWvfmZoom = 10;
    wvfm_tk = linspace(0, 1, nT*4);
    wvfm_tk_zoom = linspace(0, 1, 4*(nPointsWvfmZoom*2+1) );
    
    wvfm_tk_plot = insertNansAt(wvfm_tk(:), [1:3]*nT+1);        
    wvfm_tk_zoom_plot = insertNansAt(wvfm_tk_zoom(:), [1:3]*(nPointsWvfmZoom*2+1)+1);                
    
            
    [clustMeanWaveforms_raw_plot_orig, clustMeanWaveforms_scl_plot_orig, clustMeanWaveforms_ccw_plot_orig, clustMeanWaveforms_ccw_orig] = deal(...
     clustMeanWaveforms_raw_plot,      clustMeanWaveforms_scl_plot,      clustMeanWaveforms_ccw_plot,      clustMeanWaveforms_ccw);
    
    
    [curMeasureName, curCoordsName, curNCells, curMS, curMSsettings, curCellGrouping, ...
        curScl_idx, curScl_norm] = deal([]);                
    [prevMeasureName, prevCoordsName] = deal('');
    prevNCells = nan;
    [dWBn] = deal(nan(1,nClust));    

    allCellGroupings = cell(1,nClust);

    
%     spreadOptions = 
    spreadOptions = {'GLF'}; %{'PCA' 'GLF'};    
    coordOptionNames = {'Waveform', 'NegAmps'};    
    coordOptionVals = {clustMeanWaveforms_ccw, negAmps};

    measureOptionNames = {'CD', 'ED'};    
    measuresScale0 =  [0.1,  0.3]; 
    measureOptions_dists = {'correlation', 'euclidean'};
    NScl = 201;
    
%     [KM, KMsettings] = getKMeansStruct(coordOptionNames, coordOptionVals, measureOptionNames, measureOptions_dists, measuresScale0, NScl);    
            
    circOptions = {'circles (N)', 'ellipses (std)'};
    waveformScaleOptions = {'CCW', 'Scaled', 'Raw'};            
    wvfmColorOptions = {'Lines', 'Ori tuning'};    
    clustColors = cell(1, nClust);            
    useClusters = true(nClust, nClust);
    spreadMeasure0 = spreadOptions{1};
%     coefSpreadData = getClustSpreadData(Gid, spreadOptions, clustSpkIdxs, clustIdxs);
    
    thetas = linspace(0, 2*pi, 50);
    sd = siteDataFor(Gid);
    oris_deg = sd.ori_deg;
    spfs = sd.spPeriod_pix;
    logSpfs = log(spfs);
    if flashedOrDrifting(Gid) == 1 % flashed
        oris_deg = oris_deg * 2;  % for making polar plots (0:180 --> 0:360)
    end
    oris_rad = deg2rad(oris_deg);

%     curCoefSpreadData = coefSpreadData.(spreadMeasure0);

    circMs0 = 'circles (N)';
    circMs_idx0 = find(strcmp(circMs0, circOptions));
    nOri = length(oris_deg);    
    nSpf = length(spfs);

    multipleOris = nOri > 1;
    multipleSpfs = nSpf > 1;
        
    
%     st = strtok2(siteInfo.stimType, ':');
%     fn_ext = switchh(st, {'Movie:Flashed_Gratings', 'Grating:Spatial Frequency Batch','Grating:Orientation Batch'}, ...
%             {'movie_fg', 'grating_dSf', 'grating_dOr'});
     
%     fn_ext = iff(flashedOrDrifting(Gid) == 1, 'movie_fg', 'grating_dSf');
%     indivCellRecords_filename = [CatV1Path 'indivClustersPruned_' curSortingFeatures('') '_' fn_ext '.mat'];    
    cellDataAvailable = 0;%exist(indivCellRecords_filename, 'file');
    
    if cellDataAvailable
        allNms = arrayfun(@(cid) getName('celldata', Gid, cid), clustIds, 'un', 0);
        varlist = who('-file', indivCellRecords_filename);
        for var_i = 1:length(allNms)
            if ~any(strcmp(allNms{var_i}, varlist))
                cellDataAvailable = false;
                break;
            end            
        end
        
        if cellDataAvailable            
            S_indiv = load(indivCellRecords_filename, allNms{:});
                        
            for clust_i = length(clustIds):-1:1            
                osp = S_indiv.(allNms{clust_i}).OSP;
                
                isOptPhysDataAvailable = isfield(osp, 'R2');
                if isOptPhysDataAvailable
                    physDataOptions = {'standard', 'optimized'};
                    physData_opt(clust_i) = getPhysData(osp.R, osp.R_full, oris_rad, '.');
                    physData_std(clust_i) = getPhysData(osp.R2, osp.R_full2, oris_rad, 'o');
                else
                    physDataOptions = {'standard'};
                    physData_opt = [];
                    physData_std(clust_i) = getPhysData(osp.R, osp.R_full, oris_rad, '.');                    
                end
            end
            if isOptPhysDataAvailable
                physData_opt_orig = physData_opt;
            end
            physData_std_orig = physData_std;            
            curPhysData = physData_std;
        end
        
    end
    
    
    S_pru = load(clusterPruningFile);
    clustRefrPeriod_ms_orig = S_pru.cellsRefrPeriod_ms(clustIds_orig+1); % add 1 so that MU is #1, cell #1 is #2, etc...
    clustRefrPeriod_ms = clustRefrPeriod_ms_orig;
    3;

%         cellColrs = jet(ceil(nClust/3));
%         cellColrs = repmat(cellColrs, [3 1]);


    

    absClusterGrouping = num2cell(1:nClust);          
    absUseClusters = true(1, nClust);           
    [absClusterGrouping_disp, absClusterGrouping_disp_saved] = deal(absClusterGrouping);     
    [absUseClusters_disp, absUseClusters_saved] = deal(absUseClusters);
    unselectedColor = [.941 .941 .941];
    selectedColor = [.9, .9, 0];
    fontSize = iff(nClust < 35, 10, 7);
    useChkboxTxt = nClust >= 35;
    
    
    h_clustMergeChkbox = zeros(1,nClust);        
    h_clustMergeChkboxTxt = zeros(1,nClust);
    h_ChkboxLink = cell(1,nClust);
%     mrgChkboxesVisible = zeros(1,nClust);

    % Plot the checkboxes for showing/hiding each cell's waveform       
    nExtraButtonColumns = 3.2;
    if useChkboxTxt
        chkBox_Bot = .4;
        chkBox_Hgt = .4;        
    else
        chkBox_Bot = .1;
        chkBox_Hgt = .8;                
    end
    
    chkTxt_diff = [0, -.3, 0, -.1];
    chkBox_pos = @(i) [ (i-1+.1) / (nClust+nExtraButtonColumns), chkBox_Bot, .90/(nClust+nExtraButtonColumns), chkBox_Hgt];
    chkBoxTxt_pos = @(i) chkBox_pos(i) + chkTxt_diff;
    
    button_pos = @(i,j) [ (nClust+j-1+.1) / (nClust+nExtraButtonColumns), .05 + .5*(2-i), .90/(nClust+nExtraButtonColumns), .40];    

    
    figure(mergeControlFig); clf; set(mergeControlFig, 'name', 'Merge Clusters', 'NumberTitle', 'off');    
%     screenSz = get(0, 'screenSize');
%     h_clustMergePanel = uipanel('parent', mergeControlFig, 'units', 'pixels', 'position', [5 0 screenSz(3)-10, 30], 'backgroundcolor', [.8 .8 .8]);
    h_clustMergePanel = uipanel('parent', mergeControlFig, 'units', 'normalized', 'position', [0 0 1 1], 'backgroundcolor', [.8 .8 .8]);
    
    
    % Display the buttons for merging & splitting cells
    uicontrol('style', 'pushbutton', 'string', 'Delete', ...
              'units', 'normalized', 'position', button_pos(1,1), 'parent', h_clustMergePanel, 'callback', @deleteButton_callBack);
    uicontrol('style', 'pushbutton', 'string', 'Merge', ...
              'units', 'normalized', 'position', button_pos(1,2), 'parent', h_clustMergePanel, 'callback', @mergeButton_callBack);
    uicontrol('style', 'pushbutton', 'string', 'Split', ...
              'units', 'normalized', 'position', button_pos(1,3), 'parent', h_clustMergePanel, 'callback', @splitButton_callBack);
    h_calcButton = uicontrol('style', 'pushbutton', 'string', 'CALC', ...
              'units', 'normalized', 'position', button_pos(2,1), 'parent', h_clustMergePanel, 'callback', @recalculateAfterClustMerges);
    uicontrol('style', 'pushbutton', 'string', 'Revert', ...
              'units', 'normalized', 'position', button_pos(2,2), 'parent', h_clustMergePanel, 'callback', @revertButton_callback);    
    uicontrol('style', 'pushbutton', 'string', 'Clear', ...
              'units', 'normalized', 'position', button_pos(2,3), 'parent', h_clustMergePanel, 'callback', @clearAllSelectedButton_callBack);

          
    for ci = 1:nClust
        h_clustMergeChkbox(ci) = uicontrol('style', 'checkbox',  'parent', h_clustMergePanel, ... 
                'units', 'normalized', 'position', chkBox_pos(ci), 'fontsize', fontSize, 'fontweight', 'bold', ...
                'string', num2str(clustIds(ci)), 'userdata', ci, 'tag', 'chk', 'backgroundcolor', unselectedColor, ... 
                ... 'ForegroundColor', cellColors_dens(cellId1,:), 
                 'callback', @clustCheckBoxSelect...
                );
            
        if useChkboxTxt
            txt_pos = chkBoxTxt_pos(ci);
            h_clustMergeChkboxTxt(ci) = annotation('textbox', txt_pos, 'string', num2str(clustIds(ci)), 'lineStyle', 'none', ...
                'horiz', 'center', 'vert', 'top', 'fontsize', 10, 'backgroundcolor', unselectedColor, 'fontweight', 'bold', ...
                'userdata', ci, 'tag', 'txt', 'buttondownfcn', @clustCheckBoxSelect);            
            h_ChkboxLink{ci} = linkprop([h_clustMergeChkbox(ci), h_clustMergeChkboxTxt(ci)], {'BackgroundColor', 'visible'});
        end
            
    end
        
    3;
    positionCellSelectCheckboxes;                
    
%     absClusterGrouping = {[1, 2], 3, 4, 5, 6, 7, 8};
    3;
    


%     figure(controlFig); clf; set(controlFig, 'name', Gid_str, 'NumberTitle', 'off');
% 
%     h_KM_coordSelect = uicontrol('style', 'popupmenu', 'string', coordOptionNames, 'parent', controlFig, 'value', 1, ...
%       'units', 'normalized', 'position', [.05 .95, .12, .04], 'callback', @updateMeasureType);
%     
%     h_KM_measureSelect = uicontrol('style', 'popupmenu', 'string', measureOptionNames, 'parent', controlFig, 'value', 1, ...
%                           'units', 'normalized', 'position', [.20 .95, .12, .04], 'callback', @updateMeasureType);
% 
%                       
% %                       h_spreadSelect = uicontrol('style', 'popupmenu', 'string', spreadOptions, 'parent', controlFig, 'value', curCoefSpreadData.coefIndex, ...
% %                           'units', 'normalized', 'position', [.20 .95, .12, .04], 'callback', @updateClusterSpreadMeasure);
% 
%                       
%     h_circSelect = uicontrol('style', 'popupmenu', 'string', circOptions, 'parent', controlFig, 'value', circMs_idx0, ...
%                           'units', 'normalized', 'position', [.35 .95, .20, .04], 'callback', @updateClusterCircles);
%     uicontrol('style', 'pushbutton',  'string', 'LOAD',  'units', 'normalized', 'position', [.6 .92 .15, .07], 'parent', controlFig, 'callback', @loadMerges);
%     uicontrol('style', 'pushbutton',  'string', 'SAVE',  'units', 'normalized', 'position', [.8 .92 .15, .07], 'parent', controlFig, 'callback', @saveCurrentMerges);
% 
%     % Show spread of cluster points.
%     figureMargin = [0 0 0 .1];
%     lineColors = lines(nClust);
% 
%     h_clust_spread_ax(1) = mySubplot(2,2, 1,1, [], figureMargin);
%     hold on; box on;
% 
%     [h_clust_spread_txt, h_clust_spread_elps, h_clust_spread_circ] = deal( zeros(1,nClust) );
%     normClustSizes = clustNSpikes/max(clustNSpikes);
%     
%     minCircSize = 12;
%     maxCircSize = 30;
%     
%     V0_mn = curCoefSpreadData.clustMeans;
%     V0_xelps = curCoefSpreadData.x_elps;
%     V0_yelps = curCoefSpreadData.y_elps;
%     for cl_i = 1:nClust        
%         h_clust_spread_txt(cl_i)  = text(V0_mn{cl_i}(1), V0_mn{cl_i}(2), num2str(clustIds(cl_i)), ...
%             'color', lineColors(cl_i,:), 'vert', 'mid', 'horiz', 'cent', 'parent', h_clust_spread_ax(1));
%         h_clust_spread_circ(cl_i) = plot(V0_mn{cl_i}(1), V0_mn{cl_i}(2), 'o', 'color', lineColors(cl_i,:), ...
%             'markersize', minCircSize + (maxCircSize-minCircSize)*normClustSizes(cl_i), 'linewidth', 2);
%         h_clust_spread_elps(cl_i) = plot(V0_xelps{cl_i}, V0_yelps{cl_i}, 'color', lineColors(cl_i,:), 'visible', 'off');
%         3;
%     end        
%     set(h_clust_spread_ax, 'xtick', [], 'ytick', [], 'xlim', curCoefSpreadData.xlims_ellipse, 'ylim', curCoefSpreadData.ylims_ellipse ) ;
%     h_clust_spread_tit = title(' ');                
%     
%     % plot for histograms of all CDs or EDs
%     h_ms_hist_ax = mySubplot(2,2, 1,2, [], figureMargin);
%     h_ms_hist = bar(zeros(10,2), 1, 'stacked');        
%     set(h_ms_hist(2), 'facecolor', 'g')
%     h_scl_line(1) = line([1 1], [0 1], 'color', 'r');
% 
% 
%     % plot for quality vs # clusters.
%     h_qual_ax = mySubplot(2,2, 2,1, [], figureMargin);
%     hold on; box on;
%     h_qual = plot(1:nClust, zeros(1,nClust), 'b.-');    
%     ylabel('relative quality');
%     xlabel('# clusters');    
%     xlim([.9, nClust+.1]);
%     
% %     title(['kmeans dist : ' kmean_dist]);
% 
%     % plot for best # clusters vs scale
%     h_ncell_ax = mySubplot(2,2, 2,2, [], figureMargin);
%     h_ncell_vs_scl = plot(0, 0, 'b.-');
% %     drawVerticalLine(CD_scale); 
%     h_ncell_tit = title(' ');         
%     ylabel('best number of clusters');        
%     h_scl_line(2) = drawVerticalLine(1, 'color', 'r');
%     
%     [h_wvfm_ax, h_wvfm_tit, h_wvfm, h_otc_ax, h_otc, h_otc_cent, h_spf_ax, h_spf, h_osp_ax, h_osp_im, h_osp_tit] = deal( zeros(1,nClust) );
%     h_wvfm_0line = zeros(1,nClust);
%     h_wvfm_1lines = zeros(2,nClust);        
% %     h_wvfm_std = zeros(2,nClust);        
%     
%     h_axlink = deal( cell(1,nClust) );
% 
%     %%% WAVEFORMS & ORI-TUNING CURVES 
%     figure(wvfmsFig); clf;   set(wvfmsFig, 'name', 'Waveforms', 'NumberTitle', 'off');
%     tfToOnOff = @(tf) iff(tf, 'On', 'Off');
%     tfCellData = tfToOnOff(cellDataAvailable);
%     h_wvfmScale = uicontrol('style', 'popupmenu', 'string', waveformScaleOptions, 'parent', wvfmsFig, 'value', 2, ...
%                           'units', 'normalized', 'position', [.02 .95, .12, .04], 'callback', @updateWaveformYvals);
%     h_showWvfm_chk = uicontrol('style', 'checkbox', 'string', 'Waveforms', 'parent', wvfmsFig, 'value', showWvfm0, ...
%                           'units', 'normalized', 'position', [.15 .93, .15, .06], 'callback', @updateWvfmsVisibility_callback, 'enable', tfCellData);
%     h_showOtc_chk = uicontrol('style', 'checkbox', 'string', 'Ori', 'parent', wvfmsFig, 'value', showOTC0, ...
%                           'units', 'normalized', 'position', [.30 .93, .10, .06], 'callback', @updateOtcVisibility_callback, 'enable', tfToOnOff(cellDataAvailable && multipleOris) );
%     h_showSpf_chk = uicontrol('style', 'checkbox', 'string', 'Spf', 'parent', wvfmsFig, 'value', showSPF0, ...
%                           'units', 'normalized', 'position', [.40 .93, .10, .06], 'callback', @updateSpfVisibility_callback, 'enable', tfToOnOff(cellDataAvailable && multipleSpfs) );
%     h_showOSP_chk = uicontrol('style', 'checkbox', 'string', 'OSP', 'parent', wvfmsFig, 'value', showSPF0, ...
%                           'units', 'normalized', 'position', [.50 .93, .10, .06], 'callback', @updateOSPVisibility_callback, 'enable', tfToOnOff(cellDataAvailable && multipleOris && multipleSpfs) );
%     h_wvfmColor = uicontrol('style', 'popupmenu', 'string', wvfmColorOptions, 'parent', wvfmsFig, 'value', 1, ...
%                           'units', 'normalized', 'position', [.60 .95, .15, .04], 'callback', @updateWvfms, 'enable', tfCellData);
%     h_physDataSelect = uicontrol('style', 'popupmenu', 'string', physDataOptions, 'parent', wvfmsFig, 'value', 1, ...
%                           'units', 'normalized', 'position', [.75 .95, .10, .04], 'callback', @updateTuningCurves, 'enable', tfToOnOff(cellDataAvailable && isOptPhysDataAvailable));
%     h_wvfmZoomNeg = uicontrol('style', 'checkbox', 'string', 'Zoom', 'parent', wvfmsFig, 'value', 0, ...
%                           'units', 'normalized', 'position', [.88 .95, .11, .04], 'callback', @updateWaveformYvals, 'enable', tfCellData);
% %     h_wvfmStdDev = uicontrol('style', 'checkbox', 'string', 'Std', 'parent', wvfmsFig, 'value', 0, ...
% %                           'units', 'normalized', 'position', [.88 .90, .11, .04], 'callback', @updateWvfms); %, 'enable', 0, 'visible', 'off');
% %     h_relative_chk = uicontrol('style', 'checkbox', 'string', 'Relative', 'parent', wvfmsFig, 'value', 0, ...
% %                           'units', 'normalized', 'position', [.88 .90, .11, .04], 'callback', @updateTuningCurves, 'enable', tfCellData);                                                                 
% %     relativeScale = get(h_relative_chk, 'value');                      
%     
% %       h_showISIs = uicontrol('style', 'togglebutton', 'string', 'ISI', 'parent', wvfmsFig, 'value', 1, ...
% %                            'units', 'normalized', 'position', [.85 .95, .1, .04], 'callback', @updateISIs);
% 
% 3;
%     wvfmFigMargin = [0 0 0 .2];    
%     curOSP_on = 0;
%     M1 = floor(sqrt(nClust)); N1 = ceil(nClust/M1);
%     wrp = @(x) [x(:); x(1)];
%     
% %     oriTc_tk = linspace(0, 1, nOri);
%     
%     oris_rad_ext = wrp(oris_rad);
%         
%     for cl_i = 1:nClust
%         h_wvfm_ax(cl_i) = mySubplot(M1, N1, cl_i, [], 0, wvfmFigMargin);
%         h_wvfm(cl_i) = plot(wvfm_tk_plot, clustMeanWaveforms_scl_plot(cl_i,:), 'linewidth', wvfmWidth); hold on;
% %         h_wvfm_std(1,cl_i) = plot(wvfm_tk_plot, clustMeanWaveforms_scl_plot(cl_i,:), ':' );
% %         h_wvfm_std(2,cl_i) = plot(wvfm_tk_plot, clustMeanWaveforms_scl_plot(cl_i,:), ':' );
%         h_wvfm_0line(cl_i) = drawHorizontalLine(0, 'linestyle', '-');
%         h_wvfm_1lines(:,cl_i) = drawHorizontalLine([-1, 1], 'linestyle', ':');
%         
% %         set(h_wvfm_ax(cl_i), 'xtick', [], 'ButtonDownFcn', @updateISI_cellIdx, 'userdata', cl_i); 
%         set(h_wvfm_ax(cl_i), 'xtick', []); 
%         h_wvfm_tit(cl_i) = title(' ');
%         
%         axesToLink = h_wvfm_ax(cl_i);
%         if cellDataAvailable
%             if multipleOris 
%                 h_otc_ax(cl_i) = mySubplot(M1, N1, cl_i, [], 0, wvfmFigMargin);            
%                 [otc_x,otc_y] = pol2cart(oris_rad_ext,wrp(curPhysData(cl_i).oriTC));
%                 h_otc(cl_i) = plot(h_otc_ax(cl_i), otc_x, otc_y, '.:' );
%                 set(h_otc_ax(cl_i), 'xtick', [], 'ytick', [], 'xlim', [-1, 1], 'ylim', [-1 1], 'color', 'none', 'box', 'off', 'nextplot', 'add');
%                 h_otc_cent(cl_i) = plot(h_otc_ax(cl_i), 0, 0, 'ko', 'markersize', 4, 'markerfacecolor', 'k');                
%                 axesToLink = [axesToLink, h_otc_ax(cl_i)]; %#ok<AGROW>
%             end
%                
%             if multipleSpfs
%                 h_spf_ax(cl_i) = mySubplot(M1, N1, cl_i, [], 0, wvfmFigMargin);            
%                 h_spf(cl_i) = plot(h_spf_ax(cl_i), logSpfs, curPhysData(cl_i).spfTC, 'o:' );
%                 set(h_spf_ax(cl_i), 'xtick', [], 'ytick', [], 'xlim', lims(logSpfs, .02), 'ylim', [0 1], 'color', 'none', 'box', 'off', 'nextplot', 'add');
%                 axesToLink = [axesToLink, h_spf_ax(cl_i)]; %#ok<AGROW>
%             end
%             
%             h_axlink{cl_i} = linkprop(axesToLink, {'position', 'visible'});            
%             if multipleOris && multipleSpfs 
%                 h_osp_ax(cl_i) = mySubplot(M1, N1, cl_i, [], 0, wvfmFigMargin);         
%                 h_osp_im(cl_i) = imagesc(logSpfs, oris_deg, curPhysData(cl_i).OSP );
%                 set(h_osp_ax(cl_i), 'xtick', [], 'ytick', []);            
%                 h_osp_tit(cl_i) = title(sprintf('# %d', clustIds_orig(cl_i)));
%             end
%         end        
%         
% 
%     end
%     set([h_osp_ax, h_osp_im], 'visible', 'off');
%     h_selectCellRect = annotation('rectangle', [0 0 .1 .1], 'visible', 'off');
%     
%     if ~showOTC0
%         set([h_otc h_otc_cent], 'visible', 'off')
%     end
%     if ~showSPF0
%         set(h_spf, 'visible', 'off');
%     end
%     
%     
%     if showSpikeDensities
%         siteHasMU = any(clustIds == 0);        
% 
%         cellColors =  jet(nClust-siteHasMU);    
%         if siteHasMU
%             cellColors = [1 1 1; cellColors];
%         end                
%         
%         featureSets = {'Neg', 'PCAcuw4', 'GLFcuw4'};
%         densityTypes = {'spike density', 'cluster density'};
%         scaleTypes = {'linear', 'log'};
%         featureSet_idx0 = find(strcmp(featureSets, 'GLFcuw4'), 1);
%         nDensityBins = 64;
%         
%         densityType_prev = '';
%         densityType = densityTypes{1};
%         nSets = length(featureSets); 
%         fset_id0 = 1;
% 
%         [features, featureLabels] = getGroupFeatures(Gid, featureSets, 1, 1);        
%         
%         nMult = 3;
%         nFetsPerBasis = cellfun(@(x) size(x,2), features);        
%         combos2D = arrayfun(@(n) getNMultPairs(n, nMult), nFetsPerBasis, 'un', 0);        
%         
%         [spkDensities, axesLims2D, axesTickLims2D] = deal(cell(1,nSets));
%         for bi = 1:nSets            
%             [spkDensities{bi}, axesLims2D{bi}, axesTickLims2D{bi}] = ...
%                 getDensityHists(features{bi}, combos2D{bi}, clustSpkIdxs, nDensityBins, siteHasMU);
%         end            
%         spkDensities_orig = spkDensities;
%              
%         cmap_spikes = [0 0 0; jet(64)];
%         [cmap_clusters, cmapSclFactor_cell, cmapSclFactor_grp] = getStackedCmap(cellColors, 1);
%         
%         
%         figure(spkDensFig); clf; set(spkDensFig, 'Name', 'Cluster Densities', 'NumberTitle', 'off');
%         
%         h_featureSetSelect = uicontrol('style', 'popupmenu', 'string', featureSets, 'parent', spkDensFig, 'value', featureSet_idx0, ...
%                               'units', 'normalized', 'position', [.05 .95, .20, .04], 'callback', @updateSpikeDensities);
%         h_densityTypeSelect = uicontrol('style', 'popupmenu', 'string', densityTypes, 'parent', spkDensFig, 'value', curCoefSpreadData.coefIndex, ...
%                               'units', 'normalized', 'position', [.30 .95, .20, .04], 'callback', @updateSpikeDensities);
%         h_scaleSelect = uicontrol('style', 'popupmenu', 'string', scaleTypes, 'parent', spkDensFig, 'value', 2, ...
%                               'units', 'normalized', 'position', [.55 .95, .20, .04], 'callback', @updateSpikeDensities);
%         h_showLabels_chk = uicontrol('style', 'checkbox', 'string', 'Labels', 'parent', spkDensFig, 'value', 1, ...
%                                      'units', 'normalized', 'position', [.80 .95, .1, .04], 'callback', @updateSpikeDensityLabels);
%         h_showTicks_chk = uicontrol('style', 'checkbox', 'string', 'Ticks', 'parent', spkDensFig, 'value', 0, ...
%                                      'units', 'normalized', 'position', [.90 .95, .1, .04], 'callback', @updateSpikeDensityTicks);
%                           
%         subM = 3; subN = 2;
%         nSubplots = subM*subN;
%         [h_2D_ax, h_2D_im, h_2D_xlab, h_2D_ylab] = deal( zeros(1,nSubplots) );
%         densFigMargin = [0 0 .05 .05];
%         for plot_i = 1:nSubplots
%             h_2D_ax(plot_i) = mySubplot(subM,subN,plot_i, [], 0, densFigMargin); 
%             tks = axesTickLims2D{fset_id0}(plot_i,:);
%             h_2D_im(plot_i) = imagesc(tks(1:2), tks(3:4), zeros(nDensityBins, nDensityBins));
%             
%             [t1, t2] = deal(combos2D{fset_id0}(plot_i,1), combos2D{fset_id0}(plot_i,2));
%             axis(h_2D_ax(plot_i), axesLims2D{fset_id0}(plot_i,:) );
%             axis(h_2D_ax(plot_i), 'xy');
%             h_2D_xlab(plot_i) = xlabel(featureLabels{fset_id0}{t1}); 
%             h_2D_ylab(plot_i) = ylabel(featureLabels{fset_id0}{t2});
%             set(h_2D_ax(plot_i), 'xtick', [], 'ytick', []);
%         end
%         
%         if doColorBar
%             p6 = get(h_2D_ax(end), 'position');
%             dx = .001;
%             dummy_pos = [p6(1:2) + [p6(3), 0], dx, p6(4)];
%             h_dummy_ax = axes('position', dummy_pos); 
%             h_dummy_im = imagesc(zeros(1,2)); 
%             hColorBar = colorbar;
%             p_colbr = get(hColorBar, 'position'); 
%             p_colbr(1) = p6(1)+p6(3)+dx;
%             p_colbr(3) = p_colbr(3)*2;
%             set(hColorBar, 'position', p_colbr);
% 
%             set(h_dummy_ax, 'xtick', [], 'ytick', [], 'position', p_colbr, 'visible', 'off');
%         end
%         
%     end
%     3;
    
    %%% ISIs
    figure(isisFig); clf(isisFig); set(isisFig, 'Name', 'ISI & Pairwise Data', 'NumberTitle', 'off');
    [allRanges_ms, allNbins] = getGlobals('isi_allRanges_ms', 'isi_allNbins');        

    maxISIrange = max(allRanges_ms);
    isiFigMargin = [0 0 0 0];    
    showRefrTot = 0;
    if showRefrTot
        h_refr_tot_ax = mySubplot(11,1, 1,1, 0, isiFigMargin);
        h_refr_tot_im = imagesc(zeros(2));
    end

    h_pair_ax = mySubplot(11,1, [2,6],1, 0, isiFigMargin);
    h_pair_im = imagesc(zeros(2)); axis xy;
    h_pair_txt = zeros(nClust,nClust);
    for ci = 1:nClust
        for cj = 1:nClust
            h_pair_txt(ci,cj) = text(ci, cj, '', 'horiz', 'center', 'vert', 'mid');
        end
    end
    pairCmap = [jet(100)];
    colormap(pairCmap);
    colorbar;


    rangeOptions = arrayfun(@(rng_ms) sprintf('%d ms', rng_ms), allRanges_ms, 'un', 0);
    nbinOptions = arrayfun(@(nbin) sprintf('%d bins', nbin), allNbins, 'un', 0);
    isiTypeOptions = {'ISI', 'auto'};
    [isi_masterBinEdges, isi_masterBinIdxs] = getHistMasterBins(allRanges_ms, allNbins, 0);  
    3;
    
%     h_refrISI_range = uicontrol('style', 'popupmenu', 'string', rangeOptions, 'parent', isisFig, 'value', 3, ...
%                           'tag', 'range_ms', 'units', 'normalized', 'position', [.05 .38, .25, .04], 'callback', @updateISI_hist_NbinRange);
%     h_refrISI_nbins = uicontrol('style', 'popupmenu', 'string', nbinOptions, 'parent', isisFig, 'value', 2, ...
%                           'tag', 'nbins', 'units', 'normalized', 'position', [.35 .38, .25, .04], 'callback', @updateISI_hist_NbinRange);
%     h_refrISI_type = uicontrol('style', 'popupmenu', 'string', isiTypeOptions, 'parent', isisFig, 'value', 1, ...
%                           'tag', 'nbins', 'units', 'normalized', 'position', [.65 .38, .25, .04], 'callback', @updateISI_hist_typeClustActive);
                          
    h_refr_isi_ax = mySubplot(11,1, [8,11], 1, 0, isiFigMargin);
    h_refr_isi_bar = bar(1:2, zeros(1,2), 1);
%     h_refr_isi_line = line([0 0], [1 1], 'linestyle', ':', 'color', 'r');
    h_refr_isi_tit = title(' ');

    pairwiseOptions = {'ISI refr spikes', 'CCG Pruning Mtx', 'Refractory Periods', 'Error Matrix', 'Cluster Overlap', 'Waveform CC', 'Neg Amps CC', 'Spk Raster CC'};
%     h_pairMs_select = uicontrol('style', 'popupmenu', 'string', pairwiseOptions, 'parent', isisFig, 'value', 1, ...
%                           'units', 'normalized', 'position', [.3 .95, .4, .04], 'callback', @updatePairwiseMeasure);
    curPairwiseOption = pairwiseOptions{1};
                  
    [pairNRefractorySpikes, curPairData] = deal(zeros(nClust));
       

    curISI_cell = 1;
    curCellZoom = [];

    allRefr_isi_mb = SavedData.isiData.isis_mb;
%     refr_isi_mbIdxs = SavedData.isiData.mbIdxs;
    refr_isi_mbEdges = SavedData.isiData.mbEdges;
    
    curCellISI_mb = [];
    [ccg_mbs, ccg_mbIdxs] = deal({}, []);

            
    pairErrorMtx = SavedData.pairData.errorMtx;
    pairGaussOverlaps = SavedData.pairData.overlaps;    
    pairPruningMatrix_orig = SavedData.crossRefrData.crossPruningData.nSpksRemoved;
    pairRefrPeriods_ms_orig = SavedData.crossRefrData.crossPruningData.refrPeriod_ms;
    pairRefrPeriods_ms_orig(logical(eye(nClust))) = clustRefrPeriod_ms;
    pairPruningMatrix = pairPruningMatrix_orig;
    pairRefrPeriods_ms = pairRefrPeriods_ms_orig;
%     pairWvfmCCs = 1-KM.Waveform.CD.pairDists;
%     pairNegAmpsCCs = 1-KM.NegAmps.CD.pairDists;
        
%     applySavedData(SavedData);
    

    isi_range0 = 4;
    isi_nbin0 = 10;

    isi_range0_idx = find(isi_range0 == allRanges_ms, 1); assert(~isempty(isi_range0_idx))
    isi_nbin0      = find(isi_nbin0  == allNbins, 1);     assert(~isempty(isi_nbin0))
    
    

    % CONTROL WHICH CLUSTERS ARE VISIBLE IN CURRENT "CELL"
    figure(showClustsFig); clf;  set(showClustsFig, 'name', 'Active Clusters', 'NumberTitle', 'off');
    
    h_clustSelectPanel = uipanel('units', 'normalized', 'position', [0 0, 1, .9], 'parent', showClustsFig);        
%     uicontrol('style', 'checkbox', 'value', 1,     
    
    h_clustSelectAll = uicontrol('style', 'checkbox', 'parent', showClustsFig, 'units', 'normalized', ...
        'position', [.05 .92 .3,  .06], 'string', '', 'enable', 'on', 'value', 1, 'callback', @clustSelectAll);

    uicontrol('style', 'pushbutton', 'parent', showClustsFig, 'units', 'normalized', ...
        'position', [.4  .92 .3, .06], 'string', 'Reset', 'enable', 'on', 'value', 1, 'callback', @resetUseClusters, 'fontsize',8);
    uicontrol('style', 'pushbutton', 'parent', showClustsFig, 'units', 'normalized', ...
        'position', [.75  .92 .22, .06], 'string', 'X', 'enable', 'on', 'value', 1, 'callback', @selectActiveClusters, 'fontsize',8);

    
    h_clustSelect = zeros(1,nClust);    
    clust_chk_pos = @(i) getNormPosition (nClust, 1, i,1, [.05, .02]);
    for cj = 1:nClust
        h_clustSelect(cj) = uicontrol('style', 'checkbox', 'parent', h_clustSelectPanel, ...
            'units', 'normalized', 'position', clust_chk_pos(cj), ...
            'string', sprintf('%d (%d)', clustIds(cj), clustNSpikes(cj)), ...
            'callback', @clustSelectChkbox_callback, 'userdata', cj, 'value', 1);
    end
    
    
    % CROSS CORRELOGRAMS
    figure(ccgsFig); clf; set(ccgsFig, 'Name', 'CCGs', 'NumberTitle', 'off');
        
%             rangeOptions = arrayfun(@(rng_ms) sprintf('%d ms', rng_ms), allRanges_ms, 'un', 0);
%             nbinOptions = arrayfun(@(nbin) sprintf('%d bins', nbin), allNbins, 'un', 0);
    
    h_ccg_range = uicontrol('style', 'popupmenu', 'string', rangeOptions, 'parent', ccgsFig, 'value', 5, ...
                          'tag', 'range_ms', 'units', 'normalized', 'position', [.05 .95, .2, .04], 'callback', @updateCCGs);
    h_ccg_nbins = uicontrol('style', 'popupmenu', 'string', nbinOptions, 'parent', ccgsFig, 'value', 3, ...
                          'tag', 'nbins', 'units', 'normalized', 'position', [.3 .95, .2, .04], 'callback', @updateCCGs);
    h_ccg_xtick = uicontrol('style', 'checkbox', 'string', 'X', 'parent', ccgsFig, 'value', 1, ...
                          'tag', 'nbins', 'units', 'normalized', 'position', [.55 .95, .2, .04], 'callback', @updateCCGs);
    h_ccg_ytick = uicontrol('style', 'checkbox', 'string', 'Y', 'parent', ccgsFig, 'value', 1, ...
                          'tag', 'nbins', 'units', 'normalized', 'position', [.65 .95, .2, .04], 'callback', @updateCCGs);

%     ccg_mbEdges = SavedData.ccgData.mbEdges;
%     allCCGs_ms = SavedData.ccgData.allCCGs_ms(idx_clustsUse, idx_clustsUse);
    


    ccgFigMargin = [.01 .01 .01 .1];            
    nCCGs = nClust-1;
    [h_ccg_ax, h_ccg_bar, h_ccg_txt, h_ccg_line] = deal( zeros(nClust, nClust) );    
    for ci = 1:nClust
        for cj = 1:ci-1
            h_ccg_ax(ci, cj) = mySubplot(nCCGs, nCCGs, cj, ci-1, 0, ccgFigMargin, 'parent', ccgsFig);
            h_ccg_bar(ci, cj) = bar([1 2], [0 0], 1, 'visible', 'off');
                        
            h_ccg_txt(ci, cj) = text(0,1, ' ', 'parent', h_ccg_ax(ci, cj), 'horiz', 'center', 'vert', 'top');            
            h_ccg_line(ci, cj) = line([0 0], [1 1], 'color', 'k', 'linestyle', ':');
            set(h_ccg_ax(ci, cj), 'xtick', [], 'ytick', []);
        end
    end
    3;

     % SPIKE RASTERS
    figure(rastersFig); clf; set(rastersFig, 'name', 'Spike Rasters', 'NumberTitle', 'off');
    rastersFigMargin = [0 0 0 .07];

    rasterNBinOptions = [10, 40, 100, 400, 1000];
    corrNBin = 40;
    nnbins = length(rasterNBinOptions);
    
    rasterNBinOptions_str = arrayfun(@(nbin) sprintf('%d bins', nbin), rasterNBinOptions, 'un', 0);
    h_raster_nbins = uicontrol('style', 'popupmenu', 'string', rasterNBinOptions_str, 'parent', rastersFig, 'value', 3, ...
                               'tag', 'nbins', 'units', 'normalized', 'position', [.2 .95, .2, .04], 'callback', @updateRasterNBins);
    
    h_raster_ax = zeros(1,nClust);
    h_raster_bar = zeros(1,nClust);
    h_raster_ylab = zeros(1,nClust);
    
    raster_binVals = cell(nClust, nnbins);
    raster_binE = cell(1, nnbins);
    raster_binC = cell(1, nnbins);
    
    expEnd_ms = expDuration_sec*1000;
    
    for nbi = 1:nnbins
        raster_binE{nbi} = linspace(0, expEnd_ms, rasterNBinOptions(nbi)+1);
        raster_binC{nbi} = binEdge2cent(raster_binE{nbi});
        
        for cj = 1:nClust
            raster_binVals{cj, nbi} = uint16( histcnt(spikeTimes_ms(clustSpkIdxs{cj}), raster_binE{nbi}) );
        end
            
        if rasterNBinOptions(nbi) == corrNBin
            allRasters_i = double([raster_binVals{:, nbi}]');
            pairRasterCorr = 1-squareform( pdist(allRasters_i, 'correlation') );
        end
        3;
    end
                
    for cj = 1:nClust
        h_raster_ax(cj) = mySubplot(nClust, 1, cj, 1, 0, rastersFigMargin);         
        h_raster_ylab(cj) = ylabel(num2str(clustIds(cj)), 'parent', h_raster_ax(cj));
        h_raster_bar(cj) = bar(0, 0, 1, 'linewidth', 1);        
    end    
    set(h_raster_ax, 'xlim', [0 expEnd_ms]);  
    
    
          
    % FIGURE WITH MANUAL CLUSTER MERGING     
    function positionCellSelectCheckboxes
        nClusterGroups = length(absClusterGrouping_disp);
%         numSizCol2str = @(n,col) sprintf('\\color[rgb]{%d %d %d}%d', col(1), col(2), col(3), n);
        numTf2str = @(n,tf) iff(tf, sprintf('%d', n), sprintf('[%d]', n));
         
        nPerGrp = cellfun(@length, absClusterGrouping_disp);
        curChkBox_pos = @(i) [ ((sum(nPerGrp(1:(i-1))))+.1 ) / (sum(nPerGrp)+nExtraButtonColumns), chkBox_Bot, .95*nPerGrp(i)/(sum(nPerGrp)+nExtraButtonColumns), chkBox_Hgt];
        curChkBoxTxt_pos = @(i) curChkBox_pos(i) + chkTxt_diff;
%         curChkBox_pos = @(i) [ (i) / (nClusterGroups+3), .005, 1/(nClusterGroups+3), .04];
        for clust_grp_i = 1:nClusterGroups
            clustIdxs_inGrp = sort( absClusterGrouping_disp{clust_grp_i} );
            clustIds_inGrp = clustIds_orig( clustIdxs_inGrp );
            
            useClust = absUseClusters_disp(absClusterGrouping_disp{clust_grp_i});
            if (length(clustIds_inGrp) == 1) || (length(unique(useClust)) == 1)
                col = iff(useClust(1), 'black', 'red');
                clustList_str = cellstr2csslist(arrayfun(@num2str, clustIds_inGrp, 'un', 0));                    
            else
                col = 'black';
                clustList_str = cellstr2csslist(arrayfun(numTf2str, clustIds_inGrp, useClust, 'un', 0));    
            end
%             clustIds_cols = arrayfun(@(tf) iff(tf, 'black', 'red'), absUseClusters_disp(absClusterGrouping_disp{clust_grp_i}), 'un', 0);
%             clustList_str = 
%             clustId1 = find(clustIds == absClusterGrouping_disp{clust_grp_i}(1));

            posChk = curChkBox_pos(clust_grp_i);            
            set( h_clustMergeChkbox(clustIdxs_inGrp(1)), 'position', posChk, 'ForegroundColor', col, ...
                'string', iff(useChkboxTxt, '', clustList_str), 'visible', 'on');%, 'userdata', clustIds_inGrp);  
            set( h_clustMergeChkbox(clustIdxs_inGrp(2:end)), 'visible', 'off');%, 'userdata', []);
            
            if useChkboxTxt
                set( h_clustMergeChkboxTxt(clustIdxs_inGrp(1)), 'position', curChkBoxTxt_pos(clust_grp_i), 'Color', col, 'string', clustList_str);                
            end
            
            
            
        end
%         set(h_clustMergeChkbox(nClusterGroups+1:end), 'visible', 'off');
        
        
        absClusterGrouping_tmp = cellfun(@(x) x(absUseClusters(x)), absClusterGrouping_disp, 'un', 0);
        absClusterGrouping_tmp = absClusterGrouping_tmp( ~cellfun(@isempty, absClusterGrouping_tmp));
                
        if ~isequal(absClusterGrouping_tmp, absClusterGrouping) || ~isequal(absUseClusters_disp, absUseClusters)
            set(h_calcButton, 'backgroundcolor', [.9 .0 0])
        else
            set(h_calcButton, 'backgroundcolor', [.94 .94 .94])
        end
                
% % %         % this is the 'absClusterGrouping' that will be used in calculations.  
        % make sure that any invisible checkboxes are not selected (from a load/save)
        idx_invisible = strcmp(get(h_clustMergeChkbox, 'visible'), 'off');                 
        set(h_clustMergeChkbox(idx_invisible), 'value', 0);                
%         mrgChkboxesVisible = find(~idx_invisible);
        
    end

    function clustCheckBoxSelect(src,~)
        switch get(src, 'tag')
            case 'txt', ci = get(src, 'userdata'); 
                h_chkbox = h_clustMergeChkbox(ci);
                toggleValue(h_chkbox);                
                
            case 'chk', h_chkbox = src;
        end            
        
        onOff = get(h_chkbox, 'value');        
        col = iff(onOff, selectedColor, unselectedColor);
%         set(src, 
        set(h_chkbox, 'BackgroundColor', col);

    end

    function deleteButton_callBack(~, ~)        
        mrgChkboxesVisible = cellfun(@(x) x(1), absClusterGrouping_disp);
        
        idx_selected = find( cell2mat( get(h_clustMergeChkbox(mrgChkboxesVisible), 'value')) );                
        if isempty(idx_selected)
            return;
        end           
        
        selected_clustIdxs = [absClusterGrouping_disp{idx_selected}];
%         clustIds(selected_clustIdxs);
        
        isDeleted = ~absUseClusters_disp(selected_clustIdxs(1));
        idx_mrg_selected = mrgChkboxesVisible(idx_selected);
        if ~isDeleted  % delete now
            set(h_clustMergeChkbox(idx_mrg_selected), 'ForegroundColor', 'r');
            if useChkboxTxt
                set(h_clustMergeChkboxTxt(idx_mrg_selected), 'Color', 'r');
            end
            
            absUseClusters_disp(selected_clustIdxs) = 0;            
        else % undelete
            set(h_clustMergeChkbox(idx_mrg_selected), 'ForegroundColor', 'k');
            if useChkboxTxt
                set(h_clustMergeChkboxTxt(idx_mrg_selected), 'Color', 'k');
            end
            absUseClusters_disp(selected_clustIdxs) = 1;
        end
        positionCellSelectCheckboxes;        
                    
    end


    function mergeButton_callBack(~, ~)        
        mrgChkboxesVisible = cellfun(@(x) x(1), absClusterGrouping_disp);
        idx_selected = find( cell2mat( get(h_clustMergeChkbox(mrgChkboxesVisible), 'value')) );                
        if length(idx_selected) > 1
            set(h_clustMergeChkbox(mrgChkboxesVisible(idx_selected(1:end))), 'value', 0, 'BackgroundColor', unselectedColor);  % unselect before merging.                  
            
            mergeClustGroups(idx_selected);
            % un-select all except #1;            
        end
            
    end

    function mergeClustGroups(clustGroupIds_toMerge)
        id1 = clustGroupIds_toMerge(1);
        allCellIdsInNewGroup = [absClusterGrouping_disp{clustGroupIds_toMerge}];
        absClusterGrouping_disp{id1} = allCellIdsInNewGroup;
        absClusterGrouping_disp(clustGroupIds_toMerge(2:end)) = [];
        
        positionCellSelectCheckboxes;        
%         updateCellWaveforms;
    end



    function splitButton_callBack(~, ~)
        mrgChkboxesVisible = cellfun(@(x) x(1), absClusterGrouping_disp);
        
        idx_selected = find( cell2mat( get(h_clustMergeChkbox(mrgChkboxesVisible), 'value')) );                
        if ~isempty(idx_selected)
            idx_grpsWithMultipleCells = idx_selected(cellfun(@length, absClusterGrouping_disp(idx_selected)) > 1);
            if ~isempty(idx_grpsWithMultipleCells)
                clustIds_toBeSelectedAfter = [absClusterGrouping_disp{idx_selected}];
                
                splitClustGroups(idx_grpsWithMultipleCells);            
                
                idx_chkBoxToBeSelected = find(cellfun(@(cids) any(clustIds_toBeSelectedAfter == cids(1)),  absClusterGrouping_disp));
%                 set(h_clustMergeChkbox, 'value',0)
                set(h_clustMergeChkbox(idx_chkBoxToBeSelected), 'value',1, 'backgroundcolor', selectedColor)                 %#ok<FNDSB>
            end
        end
         
    end

    function splitClustGroups(clustGroupIds_toSplit)
        recoveredClustIds = num2cell([absClusterGrouping_disp{clustGroupIds_toSplit}]);        
        absClusterGrouping_disp(clustGroupIds_toSplit) = [];
        absClusterGrouping_disp = [absClusterGrouping_disp, recoveredClustIds];
        ascend_idx = ord( cellfun(@(x) x(1), absClusterGrouping_disp), 'ascend');
        absClusterGrouping_disp = absClusterGrouping_disp(ascend_idx);        
        
        positionCellSelectCheckboxes;                
%         updateCellWaveforms;
%         updateDendrogram;
    end


     function revertButton_callback(~,~)

         absClusterGrouping_disp = absClusterGrouping_disp_saved;
         absUseClusters_disp = absUseClusters_saved;

         positionCellSelectCheckboxes;
     end

     function clearAllSelectedButton_callBack(~,~)
        mrgChkboxesVisible = cellfun(@(x) x(1), absClusterGrouping_disp);
         
        set(h_clustMergeChkbox(mrgChkboxesVisible), 'value', 0);
        for i = mrgChkboxesVisible
            clustCheckBoxSelect(h_clustMergeChkbox(i));
        end         

     end

 
    function recalculateAfterClustMerges(~,~)
        3;
        % 1. Save the current grouping & absClustersUse (in case want to 'revert' later on) 
        absClusterGrouping_disp_saved = absClusterGrouping_disp;
        absUseClusters_saved = absUseClusters_disp;

        % 2. Compute the new data.
        if ~any(absUseClusters_disp)
            error('No clusters are being used');
        end
        inputs = struct('absClustClustering', {absClusterGrouping_disp}, 'absUseClusters', {absUseClusters_disp}, ...
            'clustMeanWaveforms_ccw', clustMeanWaveforms_ccw);        
        newData = getClusterData(Gid, 2, 2, inputs);        
        
        % 3. Compute the 'absClusterGrouping' that will be used in calculations (it has the deleted clusters removed).          
        absClusterGrouping = cellfun(@(x) x(absUseClusters_disp(x)), absClusterGrouping_disp, 'un', 0);
        absClusterGrouping = absClusterGrouping( ~cellfun(@isempty, absClusterGrouping));                
        clustIdxs = sort(cellfun(@(x) x(1), absClusterGrouping)); % a list of all the clusters in use.        
                                    
        absUseClusters = absUseClusters_disp;
        clustIds = clustIds_orig( clustIdxs );
        clustsAvailable = false(1,nClust_orig); clustsAvailable(clustIdxs) = true;        
        inactiveClustIdxs = find(~clustsAvailable);        
        nClust = length(clustIdxs);
        
        clustRefrPeriod_ms(:) = nan;
        clustRefrPeriod_ms(clustIdxs) = cellfun(@(x) min(clustRefrPeriod_ms_orig(x)), absClusterGrouping);
        
        % recompute clustSpkIdxs (list of which spikes belong to which cluster)
        clustSpkIdxs(clustIdxs) = cellfun(@(lst) sortCat( clustSpkIdxs_orig(lst) ), absClusterGrouping, 'un', 0);
        [clustSpkIdxs{inactiveClustIdxs}] = deal([]);        
        if all(absUseClusters_disp)
            assert(sum( cellfun(@length, clustSpkIdxs)) == sum( cellfun(@length, clustSpkIdxs_orig)));        
        end
        
        % recompute # spikes per cluster, and firing rates.
        clustNSpikes = cellfun(@length, clustSpkIdxs);
        clustNSpikes(inactiveClustIdxs) = nan;         
        clustFiringRates = clustNSpikes / expDuration_sec;        

        % recompute mean spike waveforms
        clustMeanWaveforms_raw_plot = updateMeanWaveforms(clustMeanWaveforms_raw_plot_orig, clustNSpikes_orig, absClusterGrouping );
        clustMeanWaveforms_scl_plot = updateMeanWaveforms(clustMeanWaveforms_scl_plot_orig, clustNSpikes_orig, absClusterGrouping );
        clustMeanWaveforms_ccw_plot = updateMeanWaveforms(clustMeanWaveforms_ccw_plot_orig, clustNSpikes_orig, absClusterGrouping );
        clustMeanWaveforms_ccw      = updateMeanWaveforms(clustMeanWaveforms_ccw_orig,      clustNSpikes_orig, absClusterGrouping );

        
        % Apply the Saved Data (after recomputing the mean spike waveforms)
        applySavedData(newData);
                
        % update cluster spread                
        coefSpreadData = getClustSpreadData(Gid, spreadOptions, clustSpkIdxs, clustIdxs);
        set([h_clust_spread_txt(inactiveClustIdxs), h_clust_spread_elps(inactiveClustIdxs), h_clust_spread_circ(inactiveClustIdxs)], 'visible', 'off');
        set([h_clust_spread_txt(clustIdxs), h_clust_spread_elps(clustIdxs), h_clust_spread_circ(clustIdxs)], 'visible', 'on');        
%         updateClusterSpreadMeasure;
        updateClusterCircles;
                
        % update histograms
        [prevMeasureName, prevCoordsName] = deal('');
        updateMeasureType;
        updateCellGrouping;
        set(h_qual_ax, 'xlim', [.9, length(clustIdxs)+.1]);
        
        % update waveforms
        allHnds = [h_wvfm; h_clustSelect]; % h_wvfm_std; 
        if multipleOris
            allHnds = [allHnds; h_otc];
        end
        if multipleSpfs
            allHnds = [allHnds; h_spf];
        end            
        set(allHnds(:,inactiveClustIdxs), 'visible', 'off');  % [h_wvfm(inactiveClustIdxs), h_wvfm_std(inactiveClustIdxs), h_otc(inactiveClustIdxs), h_spf(inactiveClustIdxs), h_clustSelect(inactiveClustIdxs)
        set(allHnds(:,inactiveClustIdxs), 'visible', 'on');
        updateWaveformYvals;
        updateWvfms;
                
        % update densities
        for b_i = 1:nSets            
            spkDensities{b_i} = updateDensities(spkDensities_orig{b_i}, absClusterGrouping);
        end

        % update Physiological data (oriTCs, spfTCs, and colors)
        physData_std = updatePhysData(physData_std_orig, absClusterGrouping, clustNSpikes_orig, 'std', oris_rad);
        if isOptPhysDataAvailable
            physData_opt = updatePhysData(physData_opt_orig, absClusterGrouping, clustNSpikes_orig, 'opt', oris_rad);
        end
        updateTuningCurves;
        
        % update pairwise-data : cross-pruning data
        pairPruningMatrix = updatePruningMatrix(pairPruningMatrix_orig, absClusterGrouping, 'plus');
        pairRefrPeriods_ms = updatePruningMatrix(pairRefrPeriods_ms_orig, absClusterGrouping, 'min');
        
        positionCellSelectCheckboxes;
%         updateCellClustersUsage_clustOnOff;
        3;
        updatePairImage_clustsActive;
    end
    
                
    
    function applySavedData(data)
        % ccgs
        ccg_mbs(clustIdxs,clustIdxs) = data.ccgData.ccgs_mb;
        ccg_mbIdxs = data.ccgData.mbIdxs;    

        allRefr_isi_mb(clustIdxs,clustIdxs) = data.isiData.isis_mb;        
        
        % kmeans:
        S_kmean = data.KmeansData;
        if max((max(abs(clustMeanWaveforms_ccw(clustIdxs,:)- S_kmean.clustMeanWaveforms_ccw(:,:))) > 1e-4)) || ~isequal(clustIds(:), S_kmean.clustIds(:))
            error('Data was computed with different starting data - please recalculate');
        end
            
        allCellGroupings = S_kmean.allCellGroupings_wvfm;        
    
        
        % with updated mean spike waveforms, recompute pairwise cluster CCs and EDs, and the
        % histograms of CCs and EDs, and the best N cells for each scale.
                
        negAmps = clustMeanWaveforms_scl_plot(:,idxNegAmps);
        coordOptionVals{1} = clustMeanWaveforms_ccw; 
        coordOptionVals{2} = negAmps;                        
        
        KM = getKMeansStruct(coordOptionNames, coordOptionVals, measureOptionNames, measureOptions_dists, measuresScale0, NScl);        
        
        nGroupings = length(allCellGroupings);
        
        for co_i = 1:length(coordOptionNames)            
            
            dWBn = nan(1,nGroupings);   
            for m_i = 1:length(measureOptionNames)
                
                s = KM.(coordOptionNames{co_i}).(measureOptionNames{m_i});

                s.dLogBWdist = nan(1,nGroupings);
                for ni = 1:nGroupings
                    cell_clustIds = allCellGroupings{ni}; % calculate Quality of this clustering                    
                    [s.dLogBWdist(ni), dWBn(ni)] = getInterIntraDiffs(cell_clustIds, s.pairDists);
                end                                                
                s.bestNCells_vs_scl = getNumClustVsScale(s.dLogBWdist, dWBn, s.allScales);
        
                KM.(coordOptionNames{co_i}).(measureOptionNames{m_i}) = s;
            end
        end
        
        % pairdata                
        pairErrorMtx(clustIdxs,clustIdxs) = log10(data.pairData.errorMtx);
        pairGaussOverlaps(clustIdxs,clustIdxs) = -data.pairData.overlaps;        
%         pairGaussOverlaps(eye(nClust)==1) = nan;
        pairWvfmCCs   (clustIdxs,clustIdxs) = 1-KM.Waveform.CD.pairDists(clustIdxs,clustIdxs);        
        pairNegAmpsCCs(clustIdxs,clustIdxs) = 1-KM.NegAmps.CD.pairDists(clustIdxs,clustIdxs);        
        
%         pairNegAmpsCCs = 1-squareform( pdist(negAmps, 'correlation') );
    end

    

        
%     curDistMeasure = 'cc';                

%     updateMeasureType;
    
    
    return;


 end