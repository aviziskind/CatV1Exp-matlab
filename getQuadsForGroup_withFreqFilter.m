function [dists_pca, dists_neg, meanWaveforms, wvfmChannel_ccs, pca_M, pca_Cov, allPcaCoeffs_C, allSpkWaveformsRaw_C, cellSpikeIdxs_C] = getQuadsForGroup_withFreqFilter(Gid, curFilter, defaultFilterMode)

    persistent allFreqData saveCount    
    
    freqDatafile = [CatV1Path 'MatlabDB_avi' filesep 'FreqDataFile.mat'];
    redo_all = 0;
    redo_current = 0;
    saveCountSpacing = 10;    
    
    getNegDistsIfMissing = 0;
    getMeanWvfmsIfMissing = 0;
    
    if isempty(allFreqData)        
        if exist(freqDatafile, 'file') && ~redo_all
            allFreqData = load(freqDatafile);        
        else
            allFreqData = struct;
        end
        saveCount = 0;
    end
    
    ccw_tf = 1;
    nComp = 8;    
    
    fld_name = sprintf('Gid_%d_%s', Gid, getVoltageFilterName(curFilter, 'measure', Gid));
    normCh = 0; whitenCh = 1; matchDB = 0;
    sort_opt = struct('voltageFilterMode_str', getVoltageFilterName(defaultFilterMode, 'measure', Gid));

    [uCellIds, cellSpikeIdx] = getCellSorting(Gid, 'cells', [], sort_opt);
    cellIds = uCellIds(uCellIds > 0);
    cellIdxs = binarySearch(uCellIds, cellIds);

    nCells = length(cellIds);
    t_ms = [-.9:.05:1.2]; t_idx0 = find(t_ms == 0, 1);
    if ((~isfield( allFreqData, fld_name) || (~isfield(allFreqData.(fld_name), 'dists_neg') && getNegDistsIfMissing) ...
                                          || (~isfield(allFreqData.(fld_name), 'meanWaveforms') && getMeanWvfmsIfMissing)) )...
        || redo_current 
    %%                
%         copyRawdatafiles(Gid);
    
        [allSpkWaveforms_raw] = getSpkWaveformsForFilter(Gid, curFilter, defaultFilterMode);
        [allPcaCoeffs] = getGroupWaveformCoefficients('PCA', Gid, nComp, 'concat', normCh, whitenCh, matchDB, curFilter);        
        
        allNegAmps_raw = squeeze( allSpkWaveforms_raw(t_idx0,:,:) )';
        
%         [spkPCA] = doPCAonSpikes(allSpkWaveforms, nComp);
        %%
                
        cellSpikeIdx_use = cellSpikeIdx(cellIdxs);

        pca_M = cellfun(@(i) mean(allPcaCoeffs(i, :), 1 ), cellSpikeIdx_use, 'un', 0);
        pca_Cov = cellfun(@(i) cov(allPcaCoeffs(i, :)), cellSpikeIdx_use, 'un', 0);
            
        neg_M   = cellfun(@(i) mean(allNegAmps_raw(i, :), 1 ), cellSpikeIdx_use, 'un', 0);
        neg_Cov = cellfun(@(i) cov(allNegAmps_raw(i, :), 1 ), cellSpikeIdx_use, 'un', 0);
        
        %%
        cell_pairs = nchoosek(1:nCells, 2);
        nPairs = size(cell_pairs, 1);

        dists_pca = zeros(nPairs, 1);
        dists_neg = zeros(nPairs, 1);

        for pr_i = 1:nPairs
            i1 = cell_pairs(pr_i, 1);
            i2 = cell_pairs(pr_i, 2);        

            M1_pca = pca_M{i1};    M2_pca = pca_M{i2};
            C1_pca = pca_Cov{i1};  C2_pca = pca_Cov{i2};

            M1_neg = neg_M{i1};    M2_neg = neg_M{i2};
            C1_neg = neg_Cov{i1};  C2_neg = neg_Cov{i2};
            
            dists_pca(pr_i) = -quadProdGaussians(M1_pca, C1_pca, M2_pca, C2_pca, 'log', [], [], 1);
            dists_neg(pr_i) = -quadProdGaussians(M1_neg, C1_neg, M2_neg, C2_neg, 'log', [], [], 1);            
        end

        meanWaveforms = cellfun(@(idx) mean(allSpkWaveforms_raw(:,:,idx), 3), cellSpikeIdx_use, 'un', 0);
        meanNegAmps = cellfun(@(idx) mean(allNegAmps_raw(idx, :), 1), cellSpikeIdx_use, 'un', 0);                
        
        s = struct('Gid', Gid, 'cellIds', cellIds, 'pca_M', {pca_M}, 'pca_Cov', {pca_Cov}, 'dists_pca', dists_pca, 'dists_neg', dists_neg, ...
            'meanWaveforms', meanWaveforms, 'meanNegAmps', meanNegAmps);
        
        allFreqData.(fld_name) = s;
        saveCount = saveCount + 1;
        
        if saveCount >= saveCountSpacing
            %%
            save(freqDatafile, '-struct', 'allFreqData');
            saveCount = 0;
        end
        
    else        
        s = allFreqData.(fld_name);        
    end
    
    if ~isfield(s, 'wvfmChannel_ccs')
        %%
        
        nC = length(s.cellIds);
        ch_idxs_best = cellfun(@(mwvfm)  indmin( min(mwvfm, [], 1) ),  s.meanWaveforms, 'un', 0);
        meanWvfms_best_C = cellfun(@(mwvfm, ch_i) mwvfm(:,ch_i), s.meanWaveforms,ch_idxs_best, 'un', 0);
        meanWvfms_best = [meanWvfms_best_C{:}];
        
        ccs_all = corr(meanWvfms_best);
        idx_lower = find(tril(ones(nC), -1));
        wvfmChannel_ccs = ccs_all(idx_lower); %#ok<FNDSB>
        s.wvfmChannel_ccs = wvfmChannel_ccs;                
        
    end
    
    if nargout >= 1
        dists_pca = s.dists_pca;
    end
    if nargout >= 2
        dists_neg = s.dists_neg;
    end
    if nargout >= 3
        meanWaveforms = s.meanWaveforms;
    end    
    if nargout >= 4
        wvfmChannel_ccs = s.wvfmChannel_ccs;
    end
    if nargout >= 5
        pca_M = s.pca_M;
    end
    if nargout >= 6
        pca_Cov = s.pca_Cov;
    end
    
    ccw_tf = 0;  % if retrieve waveforms - is for plotting.
    skip_wvfms_if_dont_have = 1;
    if nargout >= 7
        if ~exist('allSpkWaveforms_raw', 'var')
            [allSpkWaveforms_raw, ~, t_ms] = getSpkWaveformsForFilter(Gid, curFilter, defaultFilterMode, skip_wvfms_if_dont_have);
            if isempty(allSpkWaveforms_raw)
                [allPcaCoeffs_C, allSpkWaveformsRaw_C] = deal([]);
                return;
            end
                
            [allPcaCoeffs] = getGroupWaveformCoefficients('PCA', Gid, nComp, 'concat', normCh, whitenCh, matchDB, curFilter);
            [uCellIds, cellSpikeIdx] = getCellSorting(Gid, 'cells', [], sort_opt);
            cellIds = uCellIds(uCellIds > 0);

            cellIdxs = binarySearch(uCellIds, cellIds);
            cellSpikeIdx_use = cellSpikeIdx(cellIdxs);
        end
        
        allPcaCoeffs_C       = cellfun(@(idx) allPcaCoeffs(idx,:), cellSpikeIdx_use, 'un', 0);
        allSpkWaveformsRaw_C = cellfun(@(idx) allSpkWaveforms_raw(:,:,idx), cellSpikeIdx_use, 'un', 0);        
        
    end
        3;

    
end



function [spikeWaveforms_raw, spikeWaveforms_ccw, t_ms] = getSpkWaveformsForFilter(Gid, filterMode, defaultFilterMode, skip_flag)    
    
    recreateFilteredDataFiles = 0;
    redetectSpikes = 0;
    
    redoMeasurements = 0;
    redoBefore = 735492.913674;  % sprintf('%.6f', now);
    
    skip_if_dont_have_files = exist('skip_flag', 'var') && isequal(skip_flag, 1);
    
%     dfname = filteredDatafile('create', Gid, recreateFilteredDataFiles, filterMode);                
%         idx = strfind(fname_i, dfInfo.dataFileName);
%         file_ext = fname_i( idx + length(dfInfo.dataFileName) : end );         
%         file_hnd_filt(fi) = getOpenFileHandle(dfInfo, struct('fileExt', file_ext, 'precision', outputClass));        
        
    curVoltageFilter(filterMode);
%     voltageModeStr = getVoltageFilterName(filterMode, 'measure', Gid);
    
    if redetectSpikes
%         file_opt_detect = struct('voltageFilterMode_str', getVoltageFilterName(defaultFilterMode, 'detect', Gid));
%         spikesStartEndFile  = getFileName('startEnd', Gid, 0, file_opt_detect);
%         if ~exist(spikesStartEndFile, 'file') || fileOlderThan(spikesStartEndFile, redoBefore)
%             detectSpikes(Gid, defaultFilterMode);
%         end
    end

   allFetNames = {'doMahalanobisPeak', 'doNegAmps', 'doPosAmps', 'doPtPWidths', 'doFWHM', 'doChannelDelays', 'doFullWaveform'};
   shortFetNames = {'doNegAmps', 'doMahalanobisPeak', 'doFullWaveform'};
    
    
    file_opt_measure = struct('voltageFilterMode_str', getVoltageFilterName(filterMode, 'measure', Gid));
    file_opt_default = struct('voltageFilterMode_str', '');
    
    
    spikeProperties_default = load(getFileName('properties', Gid, 1, file_opt_default));
    if ~isfield(spikeProperties_default, 'keepSpikes')
        measureSpikeParameters(Gid, allFetNames, {defaultFilterMode, defaultFilterMode});
        spikeProperties_default = load(getFileName('properties', Gid, 1, file_opt_default));
        assert(isfield(spikeProperties_default, 'keepSpikes'));
    end
%     spikeWaveforms_default = load(getFileName('waveforms', Gid, 1, file_opt_default));
        
    spikePropertiesFile = getFileName('properties', Gid, 1, file_opt_measure);
    spikeWaveformsFile  = getFileName('waveforms', Gid, 1, file_opt_measure);
    
    isDefaultFilter = all(cellfun(@(fn) isequal(filterMode.(fn), defaultFilterMode.(fn)), fieldnames(filterMode)));    
    
    
    haveFiles = exist(spikePropertiesFile, 'file') && ~fileOlderThan(spikePropertiesFile, redoBefore)  ...
             && exist(spikeWaveformsFile, 'file') && ~fileOlderThan(spikeWaveformsFile, redoBefore);
    
    
    if ~haveFiles && skip_if_dont_have_files
        [spikeWaveforms_raw, spikeWaveforms_ccw, t_ms] = deal([]);
        return;        
    end
    
    if ~haveFiles || redoMeasurements
           if isDefaultFilter
               fetNames = allFetNames;
           else
               fetNames = shortFetNames;
           end
%            measureSpikeParameters(Gid, fetNames, {defaultFilterMode, defaultFilterMode}, spikeProperties_default, spikeWaveforms_default);
           measureSpikeParameters(Gid, fetNames, {defaultFilterMode, filterMode}, spikeProperties_default);
    end
    
    S_prop = load(spikePropertiesFile);
    
    if nargout >= 1
%     S_wav = load(spikeWaveformsFile);
%     allSpkWaveforms = S_wav.spikeWaveforms;
        [spikeWaveforms_raw, t_ms] = getSpikeWaveforms(Gid, [], 0, 0, 0, filterMode);
    end
%     assert(isequal(allSpkWaveforms, spikeWaveforms2_raw));
       
    if nargout >= 2
        [spikeWaveforms_ccw] = getSpikeWaveforms(Gid, [], 0, 1, 0, filterMode);
    end
        
%     spkT = S_prop.position;
    
%     t_ms = S_wav.t_ms;
    3;
    
end

% function [coeff, PCA_comps] = doPCAonSpikes(allSpkWaveforms, nComp)
%     [L, nChannels, nSpk] = size(allSpkWaveforms);
%     allSpkWaveforms = reshape(allSpkWaveforms, [L*nChannels, nSpk])';
%     
%     [coeff, PCA_comps] = doPCA(allSpkWaveforms, nComp);
%     coeff = coeff';            
% end

% function cellSpikeIdx_mod = adjustSpikeIdx(cellSpikeIdx, newIdx, cellIdxs)
%     nCells = length(cellIdxs);
%     cellSpikeIdx_mod = cell(1, nCells);
% 
%     for ci = 1:nCells
%         cell_spk_idx_orig = cellSpikeIdx{cellIdxs(ci)};
%         if ~isempty(newIdx)
%             cell_spk_idx = nonzeros(  newIdx(cell_spk_idx_orig)  );
% %                  cell_spk_idx2 = newIdx(cell_spk_idx);
%         else
%             cell_spk_idx = cell_spk_idx_orig;
%         end
%         cellSpikeIdx_mod{ci} = cell_spk_idx;
%     end
% 
% end



%{
option 1

    just skip the unnecessary spikes - have clusters of different sizes
        + easy
        - small spikes neglected

    
    if spike uncertain - skip in all files
        - need to maintain global record of which spikes to skip, & update
        with each new filter

    re-measure
        + correct,
        - complicated to implement, almost from scratch

%}

%{
%     default_spks = getSpikes(Gid, [], [], [], 0, struct('voltageFilterMode_str', ''));
%     default_spkTs = default_spks(:,1);
%     if ~isDefaultFilter
%         %%        
% %         idx_have = binarySearch(default_spkTs, S_prop.position);
%         assert(isequal(idx_have, 1:length(idx_have)) );
%     else
%         assert(length(default_spkTs) == S_prop.nSpikes)        
%     end


    
        if isequal(idx_have, 1:length(idx_have))
            newIdx = [];
        else
            % could have some extra spikes (discard), and some missing
            % spikes (mark which ones are missing)
            
            if any(diff(idx_have) == 0)
                [uIdx, idx_list, nInList] = uniqueList(idx_have);                
                idx0 = find (diff(idx_have) == 0);
                t_diffs = abs(default_spkTs(idx_have) - S_prop.position);
            end
            
            newIdx = zeros(1, length(default_spkTs));
            newIdx(idx_have) = 1:length(idx_have);        
        end


%}