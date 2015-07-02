function compareFrequencyFilters(Gid, cellIds)    
    
%     
%     if length(Gid) > 1
%         for gid_i = 1:length(Gid)
%             fprintf('\n-----%d/%d/------\n', gid_i, length(Gid));
%             compareFrequencyFilters(Gid(gid_i))
%         end
%         return;
%     end
% 
%     
    if nargin < 1
        %%
        Gids_all = getAllGids('f');
%         Gids_all = [4570, 2597, 4254, 1648];

%         maxGid = 3041;        
%         Gids = Gids(Gids <= maxGid);
               

        sd = siteDataFor('Gid', Gids_all, 1);
        nCells_eachGroup = arrayfun(@(s) nnz(s.cellIds), sd);

        minNCells = 6;
%         minNCells = 2;
        Gids = Gids_all(nCells_eachGroup >= minNCells);    
        Gids = Gids(1:5);
%        
%         Gids = Gids(1:20);
%         Gids = Gids(1);
%         Gids = 2288;
    end

    testMode = 'filter';
%     testMode = 'interp';
    
    t_ms = [-.9:.05:1.2]; idx_t0 = find(t_ms==0,1);

    switch testMode
        case 'filter'                        
            freqs = [25, 35, 50, 75, 105, 150, 210, 300, 425, 600, 850, 1200, 1700, 2400];
            for freq_i = 1:length(freqs)
                filters_S.(sprintf('butter_%d', freqs(freq_i))) = struct('highPass_freq_Hz', freqs(freq_i), 'filterOrder', 1,   'filterName', 'butter');
            end
%             for freq_i = 1:length(freqs)                
%                 filters_S.(sprintf('median_%d', freqs(freq_i))) = struct('highPass_freq_Hz', freqs(freq_i), 'filterOrder', 1,   'filterName', 'median');
%             end                        
            
            allFilterNames = fieldnames(filters_S);
            for i = 1:length(allFilterNames)
                filters_S.(allFilterNames{i}).interpMethod = 'sinc';
                filters_S.(allFilterNames{i}).interpN = 16;
            end            
            allFilters = struct2array(filters_S);
            nFilters = length(allFilters);
            
        case 'interp'            
%             filters_S.noItp  = struct('interpMethod', 'sinc',  'interpN', 1);
            nInterps = [1, 2, 3, 4, 6, 8, 16];  % 1 == no interp
            for itp_i = 1:length(nInterps)
                filters_S.(sprintf('sinc_%d', nInterps(itp_i))) = struct('interpMethod', 'sinc', 'interpN', nInterps(itp_i));
            end
            for itp_i = 1:length(nInterps)
                filters_S.(sprintf('spline_%d', nInterps(itp_i))) = struct('interpMethod', 'spline', 'interpN', nInterps(itp_i));
            end            
            
            allFilters = struct2array(filters_S);
            allFilterNames = fieldnames(filters_S);
            for i = 1:length(allFilterNames)
                filters_S.(allFilterNames{i}).filterName = 'butter';
                filters_S.(allFilterNames{i}).highPass_freq_Hz = 300;
                filters_S.(allFilterNames{i}).filterOrder = 1;
            end                        
            
%             filter_butter_600 = struct('highPass_freq_Hz', 600,   'filterOrder', 1,   'filterName', 'butter');
%             filter_med_300    = struct('highPass_freq_Hz', 300,   'filterOrder', 1,   'filterName', 'median');
                    
            nFilters = length(allFilters);
            [allFilters.highPass_freq_Hz] = deal(300);
            [allFilters.filterOrder] = deal(1);
            [allFilters.filterName] = deal('butter');
            
    end
    
    defaultFilterMode = curVoltageFilter('default');
    nGids = length(Gids);
    [dists_pca_grp_flt_C, dists_neg_grp_flt_C, mWvfms, wvfm_ccs, pcaM, pcaCov, cell_pca_C, cell_wvfms_C] =  deal( cell(nGids, nFilters) );
    dists_flt_S = struct;
    
    getAllWvfms = 1;
    
%     progressBar('init-', nGids);
    filt_file = [CatV1Path 'MatlabDB_avi' filesep 'freq_data.mat'];   redo_file = 0;
    if exist(filt_file, 'file') && ~redo_file
        tic; filt_file_S = load(filt_file); toc;
    else
        filt_file_S = struct;
    end
    saveCount = 0;
    
    for gi = 1:nGids
        Gid = Gids(gi);        
        fprintf('\n----- (Gid = %d) %d/%d ------\n', Gid, gi, nGids);        
                
        
        for fi = 1:nFilters                    
            flt_field = sprintf('Gid_%d_%s', Gid, getVoltageFilterName(allFilters(fi), 'measure', 'Gid') );
            if getAllWvfms
                if isfield(filt_file_S, flt_field)
                    s = filt_file_S.(flt_field);
                    [dists_pca_grp_flt_i, dists_neg_grp_flt_i, mWvfms_i, wvfm_ccs_i, pcaM_i, pcaCov_i, cell_pca_i, cell_wvfms_i] = deal(...
                        s.dists_pca_grp_flt, s.dists_neg_grp_flt, s.mWvfms, s.wvfm_ccs, s.pcaM, s.pcaCov, s.cell_pca, s.cell_wvfms);                            
                else                    
                    [dists_pca_grp_flt_i, dists_neg_grp_flt_i, mWvfms_i, wvfm_ccs_i, pcaM_i, pcaCov_i, cell_pca_i, cell_wvfms_i] = ...
                        getQuadsForGroup_withFreqFilter(Gid, allFilters(fi), defaultFilterMode);
%%
                    filt_file_S.(flt_field) = struct('dists_pca_grp_flt', dists_pca_grp_flt_i, ...
                        'dists_neg_grp_flt', dists_neg_grp_flt_i, 'mWvfms', {mWvfms_i}, 'wvfm_ccs', wvfm_ccs_i, ...
                        'pcaM', {pcaM_i}, 'pcaCov', {pcaCov_i}, 'cell_pca', {cell_pca_i}, 'cell_wvfms', {cell_wvfms_i});
                    saveCount = saveCount+1;
                end

            elseif ~getAllWvfms
                [dists_pca_grp_flt_i, dists_neg_grp_flt_i, mWvfms_i, wvfm_ccs_i, pcaM_i, pcaCov_i] = ...
                    getQuadsForGroup_withFreqFilter(Gid, allFilters(fi), defaultFilterMode);
                [cell_pca_i, cell_wvfms_i] = deal([]);
            end

            3;
            if getAllWvfms
                fprintf('[%s]', allFilterNames{fi});
                if mod(fi,9) == 0
                    fprintf('\n');
                end
            end
            
            [dists_pca_grp_flt_C{gi, fi}, dists_neg_grp_flt_C{gi, fi}, mWvfms{gi, fi}, wvfm_ccs{gi,fi}, pcaM{gi, fi}, pcaCov{gi, fi}, ...
                cell_pca_C{gi,fi}, cell_wvfms_C{gi,fi}] = deal( ...
                dists_pca_grp_flt_i, dists_neg_grp_flt_i, mWvfms_i, wvfm_ccs_i, pcaM_i, pcaCov_i, cell_pca_i, cell_wvfms_i);

            3;
        end
            
            
%         progressBar(gi);
    end
    if saveCount > 0
        %%
        save(filt_file, '-struct', 'filt_file_S', '-v6');
    end
    
    for fi = 1:nFilters
        dists_flt_S.(allFilterNames{fi}) = cat(1, dists_pca_grp_flt_C{:,fi});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    if strcmp(testMode, 'filter')
%         flt_type = 'median';
        flt_type = 'butter';
        
        X = freqs;
        fld_names = arrayfun(@(f) sprintf('%s_%d', flt_type, f), freqs, 'un', 0);
        
    elseif strcmp(testMode, 'interp')
%         interp_method = 'sinc';
        interp_method = 'spline';
        
        X = nInterps;
        fld_names = arrayfun(@(n) sprintf('%s_%d', interp_method, n), nInterps, 'un', 0);
                
    end
    
    nEachGrp = cellfun(@(x) size(x,1), dists_pca_grp_flt_C(:,1));
    fld_idxs = cellfun(@(fld_nm) find( strcmp(allFilterNames, fld_nm), 1), fld_names);
    nFilt = length(fld_names);
    
    dists_pca_C = cell(1, nFilt);
    for fi = 1:nFilt
        dists_pca_C{fi} = cat(1, dists_pca_grp_flt_C{:, fld_idxs(fi)});
    end
    dists_pca_all = [dists_pca_C{:}];
    %%
    [max_dist_pca, ind_max] = max(dists_pca_all, [], 2);
    dist_pca_nrm_max = bsxfun(@rdivide, dists_pca_all, max_dist_pca);
    ratio_lo_hi = dists_pca_all(:,1) ./ dists_pca_all(:,end);
    
    %%
    dists_pca_all_sm = gaussSmooth(dist_pca_nrm_max, 4, 2);
    [max_dist_pca_sm, ind_max_sm] = max(dists_pca_all_sm, [], 2);
    [uIdx, idx_indices] = uniqueList(ind_max_sm);    
        
    wccs = cat(1, wvfm_ccs{:,1});
    
    m = cellfun(@(i) mean(wccs(i)), idx_indices)
    s = cellfun(@(i) std(wccs(i)), idx_indices);
    
    
  %%
    imax1 = find(ind_max == 1)';
    imax14 = find(ind_max == 14)';
    
    i_best_lo = find(ratio_lo_hi > 1);
    i_best_hi = find(ratio_lo_hi < 1);
    
    idx_all = imax14(max_dist_pca(imax14)>20);
    
    cumNGroups = [0; cumsum(nEachGrp)]+1;
%%
    idx_glob = idx_all(2);
    
    Gid_idx = find(idx_glob >= cumNGroups, 1, 'last');
    pr_idx = idx_glob-cumNGroups(Gid_idx)+1;
    
%     Gid_idx = 1;
%     pr_idx = 13;
    
    %             for freq_i = 1:length(freqs)    
    
%     dists_pca = zeros(1, nFilt);
%     dists_neg = zeros(1, nFilt);
    wvfm_ccs = zeros(1, nFilt);
    neg_ccs = zeros(1, nFilt);
    
    nCells = nCells_eachGroup ( Gids(Gid_idx) == Gids_all );
    combos = nchoosek(1:nCells, 2);
    
    figure(1); clf; 
    ax1 = subplotGap(1,2,1,1); hold on; box on; title('Cell 1 Waveform');
    ax2 = subplotGap(1,2,1,2); hold on; box on; title('Cell 2 Waveform');

    
    %             figure(5); clf; h_pca(1) = plot(0,0, 'bo'); hold on; h_pca(2) = plot(0,0, 'ro');
    %             set(h_pca, 'markersize', 3);
    %             h_ax5 = gca;
    
    figure(6); clf; h_pcapca(1) = plot(0,0, 'bo'); hold on; h_pcapca(2) = plot(0,0, 'ro');
    set(h_pcapca, 'markersize', 3);
    h_ax6 = gca;
    
    figure(7); clf; h_proj = bar(zeros(2), 1, 'stacked'); hold on; h_gs = plot([0 0], zeros(2));
    set(h_proj(1), 'facecolor', 'b'); set(h_proj(2), 'facecolor', 'r');
    set(h_gs(1), 'color', 'b'); set(h_gs(2), 'color', 'r');
    bar_ax = gca;
    
    figure(15); clf; 
    h_ba(1) = subplot(1,2,1);
    h_ba(2) = subplot(1,2,2);
    set(h_ba, 'nextplot', 'add');
    
    
    
    gs_func = @(beta, x) gaussian(x, beta(1), beta(2));    
    
    use_vecOfFirst = 0;
    pca_comps = [];
    normalizeWvfmsToNeg = 1;
    plotConcatWvfm = 0;
    usePCAforProj = 1;
        pca_comps = [];
    
    L = [];
       
%     for fi = 1:nFilt
    
%         dists_pca(fi) = dists_pca_grp_flt_C{Gid_idx, fld_idxs(fi)}(pr_idx);
%         dists_neg(fi) = dists_neg_grp_flt_C{Gid_idx, fld_idxs(fi)}(pr_idx);
    col1 = [0 0 1];
    col2 = [0 .75 0];
    
    h_wvfm_plot1 = [];
    h_wvfm_plot2 = [];
    
    
    blues = bsxfun(@times, [0 0 1], [1:nFilt]'/nFilt);
    greens = bsxfun(@times, [0 1 0], [1:nFilt]'/nFilt);

    
    dists_pca = arrayfun(@(fidx) dists_pca_grp_flt_C{Gid_idx, fidx}(pr_idx), fld_idxs);        
    dists_neg = arrayfun(@(fidx) dists_neg_grp_flt_C{Gid_idx, fidx}(pr_idx), fld_idxs);
    
    figure(3); clf; hold on; box on;
    plot(X, dists_pca, 'ko-');
    plot(X, dists_neg, 'rs-');
    ylims_dist = lims([dists_pca, dists_neg], .05);
%     plot(X, dists_chk, 'r.-');
%     set(gca, 'xscale', 'log', 'xlim', lims(X, .1, [], 1));
    set(gca, 'xscale', 'log')
    set(gca, 'xlim', lims(X, .1, [], 1), 'ylim', [0, ylims_dist(2)]);
    
    h_cur_dist_wvfm = plot(gca,0,0, 'ko', 'visible', 'on');
    h_cur_dist_neg = plot(gca,0,0, 'rs', 'visible', 'on');
    xlabel('Highpass filter cutoff frequency (Hz)');
    ylabel('Overlap');
    legend({'Neg Amps', 'Waveform'}, 'location', 'bestOutside')
    
    if plotConcatWvfm
        set([ax1, ax2], 'xtick', []);
    end
    
    %%
    for fi = [1:nFilt]
        %%
%         fi = nFilt;
        
        
        i1 = combos(pr_idx, 1); i2 = combos(pr_idx, 2);
        if getAllWvfms
            all_wvfms_1 = cell_wvfms_C{Gid_idx,fld_idxs(fi)}{i1};
            all_wvfms_2 = cell_wvfms_C{Gid_idx,fld_idxs(fi)}{i2};
        
            mean_wvfms_1 = mean(all_wvfms_1, 3); ch1 = indmax(abs(mean_wvfms_1(idx_t0,:)));
            mean_wvfms_2 = mean(all_wvfms_2, 3); ch2 = indmax(abs(mean_wvfms_2(idx_t0,:)));
        else
            mean_wvfms_1 = mWvfms{Gid_idx,fld_idxs(fi)}{i1}; ch1 = indmin(mean_wvfms_1(idx_t0,:));            
            mean_wvfms_2 = mWvfms{Gid_idx,fld_idxs(fi)}{i1}; ch2 = indmin(mean_wvfms_2(idx_t0,:));
        end     
        mean_negAmps_1 = mean_wvfms_1(idx_t0,:);
        mean_negAmps_2 = mean_wvfms_2(idx_t0,:);
        
        if plotConcatWvfm
            dt = diff(t_ms(1:2));
            t_ms_plot = t_ms(1) + [0: (length(t_ms)*4 - 1)]*dt;
            mean_wvfms_1_ch = mean_wvfms_1(:);
            mean_wvfms_2_ch = mean_wvfms_2(:);
            
            if normalizeWvfmsToNeg
                mean_wvfms_1_ch = mean_wvfms_1(:) / max(abs(mean_wvfms_1_ch(:)));
                mean_wvfms_2_ch = mean_wvfms_2(:) / max(abs(mean_wvfms_2_ch(:)));                                
            end
            
        else
            t_ms_plot = t_ms;
            mean_wvfms_1_ch = mean_wvfms_1(:,ch1);
            mean_wvfms_2_ch = mean_wvfms_2(:,ch2);
            
            if normalizeWvfmsToNeg
                mean_wvfms_1_ch = mean_wvfms_1_ch / abs(mean_wvfms_1_ch(idx_t0));
                mean_wvfms_2_ch = mean_wvfms_2_ch / abs(mean_wvfms_2_ch(idx_t0));
            end
        end
        
        wvfm_ccs(fi) = pearsonR(mean_wvfms_1_ch, mean_wvfms_2_ch);
        neg_ccs(fi) = pearsonR(mean_negAmps_1, mean_negAmps_2);
        
        line_w = iff(fi == 1 || fi == nFilt, 2, 1);
        
        
        
        h_wvfm_plot1(fi) = plot(ax1, t_ms_plot, mean_wvfms_1_ch, 'color', 'k', 'linewidth', 2);
        h_wvfm_plot2(fi) = plot(ax2, t_ms_plot, mean_wvfms_2_ch, 'color', 'k', 'linewidth', 2);
        for fj = 1:fi-1
            set(h_wvfm_plot1(fj), 'color', blues(fj,:), 'linewidth', 1)
            set(h_wvfm_plot2(fj), 'color', greens(fj,:), 'linewidth', 1)
        end
        
        set(h_cur_dist_wvfm, 'xdata', X(fi), 'ydata', dists_pca(fi), 'marker', 'o', 'markersize', 10, 'markerfacecolor', 'k');
        set(h_cur_dist_neg, 'xdata', X(fi), 'ydata', dists_neg(fi), 'marker', 's', 'markersize', 10, 'markerfacecolor', 'r');    
        
        
        set([ax1, ax2], 'xlim', lims(t_ms_plot));
        %                 plot(ax1, t_ms, mean_wvfms_1, ':');
        %                 plot(ax2, t_ms, mean_wvfms_2, ':');
        
        if fi == 1 || fi == nFilt
            if fi == 1
                plot_i = 1; 
            elseif fi == nFilt
                plot_i = 2;
            end
            
            plot(h_ba(plot_i), t_ms_plot, mean_wvfms_1_ch, 'color', col1);
            plot(h_ba(plot_i), t_ms_plot, mean_wvfms_2_ch, 'color', col2);            
            set(h_ba(plot_i), 'xlim', lims(t_ms_plot));
            title(h_ba(plot_i), sprintf('cc = %.2f', wvfm_ccs(fi)));            
            
        end
            
            
        
        if getAllWvfms
            if use_vecOfFirst
                if fi == 1
                    %%
                    allSpks_used = cat(3, cell_wvfms_C{Gid_idx,fld_idxs(fi)}{:});

                    [pca_coeff_all, PCA_comps] = doPCAonSpikes(allSpks_used, 8);
%                     all_pca_saved = cell_pca_C{Gid_idx,fld_idxs(fi)};
%                         assert(isequal(pca_comp_all, all_pca_saved));    

                end
                all_pca_1 = doPCAonSpikes(all_wvfms_1, 8, PCA_comps)';
                all_pca_2 = doPCAonSpikes(all_wvfms_2, 8, PCA_comps)';
                                

            else
                all_pca_1 = cell_pca_C{Gid_idx,fld_idxs(fi)}{i1};
                all_pca_2 = cell_pca_C{Gid_idx,fld_idxs(fi)}{i2};
            end
            idx_dim_use = 1:8;
            
            M1 = mean(all_pca_1(:,idx_dim_use)); C1 = cov(all_pca_1(:,idx_dim_use));
            M2 = mean(all_pca_2(:,idx_dim_use)); C2 = cov(all_pca_2(:,idx_dim_use));
        
            pca_1 = all_pca_1(:,1:2);
            pca_2 = all_pca_2(:,1:2);
        
        else
            
            

%             dists_chk(fi) = -quadProdGaussians(M1, C1, M2, C2, 'log', [], [], 1);
        end
        
        %                 set(h_pca(1), 'xdata', pca_1(:,1), 'ydata', pca_1(:,2));
        %                 set(h_pca(2), 'xdata', pca_2(:,1), 'ydata', pca_2(:,2));
        
        [pca_1_proj, pca_2_proj, pca_1_pca, pca_2_pca, pca_comps] = projectToConnectingLine(all_pca_1, all_pca_2, usePCAforProj, pca_comps);
        
        [x_el1, y_el1] = ellipsoidFromCov(pca_1_pca(:,1:2), [], 3, 100);
        [x_el2, y_el2] = ellipsoidFromCov(pca_2_pca(:,1:2), [], 3, 100);
        h_elps1(fi) = plot(h_ax6, x_el1, y_el1, color_s(fi), 'k', 'linewidth', 2);
        h_elps2(fi) = plot(h_ax6, x_el2, y_el2, color_s(fi), 'k', 'linewidth', 2);
        for fj = 1:fi-1
            set(h_elps1(fj), 'color', blues(fj,:), 'linewidth', 1)
            set(h_elps2(fj), 'color', greens(fj,:), 'linewidth', 1)
        end
        
        %                 mPCA_1 = mean(pca_1, 1); mPCA_2 = mean(pca_1, 1);
%         if fi <= 10
            L = lims([pca_1_proj(:); pca_2_proj(:)], .1);
            L = [-.75, 1.75];
            binE = linspace(L(1), L(2), 50);
            binC = binEdge2cent(binE);
            
            xL1 = lims([pca_1(:,1); pca_2(:,1)], .05);
            yL1 = lims([pca_1(:,2); pca_2(:,2)], .05);
            
            %                     axis(h_ax5, [xL1; yL1]);
            
            xL2 = lims([pca_1_pca(:,1); pca_2_pca(:,1)], .05);
            yL2 = lims([pca_1_pca(:,2); pca_2_pca(:,2)], .05);
            
            axis(h_ax6, [xL2; yL2]);
                        
%         end
        set(bar_ax, 'xlim', binE([1, end]), 'xtickmode', 'auto');
        dBin = diff(binE(1:2));
        n1 = histcnt(pca_1_proj, binE);
        n2 = histcnt(pca_2_proj, binE);
        
        m1 = mean(pca_1_proj); s1 = std(pca_1_proj);
        m2 = mean(pca_2_proj); s2 = std(pca_2_proj);
        binC_fine = binC(1):dBin/10:binC(end);
        gs1 = gaussian(binC_fine, m1, s1);
        gs2 = gaussian(binC_fine, m2, s2);
        
        normHists = 1;
        if normHists
            n1 = n1/(sum(n1)*dBin);
            n2 = n2/(sum(n2)*dBin);
        end
        
%                         set(h_proj(1), 'xdata', binC, 'ydata', n1);
%                         set(h_proj(2), 'xdata', binC, 'ydata', n2);
        %                 plot(bar_ax, binC, gs1
        set(h_gs(1), 'xdata', binC_fine, 'ydata', gs1/max(gs1), 'color', col1, 'linewidth', 2)
        set(h_gs(2), 'xdata', binC_fine, 'ydata', gs2/max(gs2), 'color', col2, 'linewidth', 2);
        
                
        
        %                 set(h_pca(1), 'xdata', pca_1(:,1), 'ydata', pca_1(:,2));
        %                 set(h_pca(2), 'xdata', pca_2(:,1), 'ydata', pca_2(:,2));
        set(h_pcapca(1), 'xdata', pca_1_pca(:,1), 'ydata', pca_1_pca(:,2), 'color', col1);
        set(h_pcapca(2), 'xdata', pca_2_pca(:,1), 'ydata', pca_2_pca(:,2), 'color', col2);
        
        drawnow;
        pause(.1);
        
        
    end
    
    %%
    
    3;
            
            
%     figure(8); clf; hold on; box on;
%     wvfm_ccs_inv = 1-wvfm_ccs;
%     plot(X, wvfm_ccs_inv, 'bo-'); 
%     neg_ccs_inv = 1-neg_ccs;
%     plot(X, neg_ccs_inv/neg_ccs_inv(1)*wvfm_ccs_inv(1), 'r.-'); 
%     set(gca, 'xscale', 'log');
    
    
end


% function [M, C] = getMeanCov(X)
%     M = mean(X);
%     C = cov(X);
% end
    

function [coeff, PCA_comps] = doPCAonSpikes(spks, nComp, PCA_comps_prev)
    [L, nChannels, nSpk] = size(spks);
    spks = reshape(spks, [L*nChannels, nSpk])';
    
    if (nargin >= 3) && ~isempty(PCA_comps_prev)
        PCA_comps = PCA_comps_prev;
        spks_msub = bsxfun(@minus, spks, mean(spks, 1));
        coeff = [PCA_comps' * spks_msub'];        
        
    else
        [coeff, PCA_comps] = doPCA(spks, nComp);
        coeff = coeff';
        
%         spks_msub = bsxfun(@minus, spks, mean(spks, 1));
        spks_msub = spks; 
        coeff2 = [PCA_comps' * spks_msub'];        
%         assert(isequal(coeff, coeff2'));
        
    end
end



%{
            figure(55); clf;    
        %     plot(dists_med_300, dists_butter_300, '.'); axis square; hold on; fplot(@(x) x, xlim, 'k:')
        %     xlabel('Median (300 Hz)'); ylabel('Butter (300 Hz)');
            diffs = (dists_med_300 - dists_butter_300)./dists_med_300 * 100;
            hist(diffs, 30);
            xlabel('med 300 - butter 300');
            drawVerticalLine(mean(diffs), 'color', 'r');
            title(sprintf('%.2f', mean(diffs)));

            figure(56); clf;    
        %     plot(dists_butter_300, dists_butter_600, '.'); axis square; hold on; fplot(@(x) x, xlim, 'k:')
            diffs2 = (dists_butter_300 - dists_butter_600)./dists_butter_300 * 100;
            hist(diffs2, 30);
            xlabel('butter 300 - butter 600');
            drawVerticalLine(mean(diffs2), 'color', 'r');
            title(sprintf('%.2f', mean(diffs2)));
        %     xlabel('Butter (300 Hz)'); ylabel('Butter (600 Hz)');
        3;
    
        case 'interp'
%             allFilters = [filter_noItp,   filter_sinc_4, filter_sinc_16, filter_spline_4, filter_spline_16];
            
            dists_noItp     = dists_flt_S.noItp;
            dists_sinc_4    = dists_flt_S.sinc_4;
            dists_sinc_16   = dists_flt_S.sinc_16;
            dists_spline_4  = dists_flt_S.spline_4;
            dists_spline_16 = dists_flt_S.spline_16;
                        
            %%            
            figure(155); clf;    
        %     plot(dists_med_300, dists_butter_300, '.'); axis square; hold on; fplot(@(x) x, xlim, 'k:')
        %     xlabel('Median (300 Hz)'); ylabel('Butter (300 Hz)');
            diffs = (dists_sinc_4 - dists_spline_4)./dists_spline_4 * 100;
            hist(diffs, 30);
            xlabel('sinc 4 - spline 4');
            drawVerticalLine(mean(diffs), 'color', 'r');
            title(sprintf('%.2f', mean(diffs)));

            figure(156); clf;    
        %     plot(dists_butter_300, dists_butter_600, '.'); axis square; hold on; fplot(@(x) x, xlim, 'k:')
            diffs2 = (dists_sinc_16 - dists_spline_16)./dists_spline_16 * 100;
            hist(diffs2, 30);
            xlabel('sinc_16 - spline_16');
            drawVerticalLine(mean(diffs2), 'color', 'r');
            title(sprintf('%.2f', mean(diffs2)));
        %     xlabel('Butter (300 Hz)'); ylabel('Butter (600 Hz)');
            
            
            figure(157); clf;    
        %     plot(dists_med_300, dists_butter_300, '.'); axis square; hold on; fplot(@(x) x, xlim, 'k:')
        %     xlabel('Median (300 Hz)'); ylabel('Butter (300 Hz)');
            medians = [median(dists_noItp), median(dists_sinc_4), median(dists_sinc_16)];
            means = [mean(dists_noItp), mean(dists_sinc_4), mean(dists_sinc_16)];
            x = [0, 4, 16];
            plot(x, medians, 'bo-', x, means, 'gs:');            
            3;
            
    end
%}

%{

%%
%     file_hnd_raw = getOpenFileHandle(dfInfo, struct('fileExt', '', 'precision', outputClass));
%     %%
    
    
%     nInterp = 10;
% %     paramOpt.interpMethod = 'fourier'; % 'fourier' 'sinc', or 'spline'    
%     interpMethod = 'sinc'; 
    
%      allFieldNames = {'highPass_freq_Hz', 'filterOrder', 'filterName', ...
%             'window_ms', 'useVoltageRMS', 'useVoltageNEO', 'bartlettN', 'threshold'};
    
    %%
%     filtName = 'butter'; highPass_freq_Hz = 300;    
%     filtName = 'median'; highPass_freq_Hz = 300;    
    

%     voltageFilterMode = curVoltageFilter('filterName', filtName, 'highPass_freq_Hz', highPass_freq_Hz, 'itpMethod', interpMethod, 'interpN', nInterp);
%     voltageFilter_str = getVoltageFilterName(voltageFilterMode, 1, Gid);    
%     %%
    
%     spkWvforms = getSpkWaveformsForFilter(Gid, voltageFilterMode);
%         voltageFilter_str = getVoltageFilterName(voltageFilterMode, 1);    
%         spikeWaveformsFile  = getFileName('waveforms', Gid, 1, struct('voltageFilter_str', voltageFilter_str));    
%     
%     [allSpikeParams, allSpkWaveforms] = measureSpikeParameters(Gid, {'doFullWaveform'}, voltageFilterMode, paramOpt);
    
    3;
    %%
    %%
    return;
    
    nSpksEachCell = cellfun(@length, spkTimes_tk);

    
    t_ms_rep = repmat(t_ms', 1, 4);
    nSpk = size(allSpkWaveforms_C{1}, 3);
    %%
    
    eq12 = zeros(1, nSpk);
    for i = 1:nSpk
        eq12(i) = isequal(allSpkWaveforms_C{1}(:,:,i), allSpkWaveforms_C{2}(:,:,i));
    end
    
    %%
    eq12 = zeros(1, nSpk);
    for i = 1:nSpk
        eq12(i) = isequal(spikeWaveforms(:,:,i), S_wav.spikeWaveforms(:,:,i));
    end    
    %%
    idx = 15413;    
    figure(3); clf;
            
    spks_ci1 = double( allSpkWaveforms_C{1}(:,:,idx) );
    spks_ci2 = double( allSpkWaveforms_C{2}(:,:,idx) );
    plot(t_ms_rep, spks_ci1); hold on;
    plot(t_ms_rep-0, spks_ci2, ':');
%     subplot(1,3, 3);  plot(t_ms_rep, spks_ci1-spks_ci2);
        
        
    %%
    
    for fi = 1:nFilters
        figure(fi); clf;
        for ci = 1:nCells        
            subplot(1,nCells, ci);
            spks_ci = double( allSpkWaveforms_C{fi}(:,:,cellSpikeIdx{cellIdxs(ci)}) );
            errorbar(t_ms_rep, mean(spks_ci, 3), std(spks_ci, [], 3));
%             plot(t_ms_rep, mean(spks_ci, 3));

        end
    end
    %%
        
    for ci = 1:nCells
        for spk_i = 1:nSpksEachCell(ci)
            %%
            tk = spkTimes_tk{ci}(spk_i);
            
            samp_start = rectified( round(tk - nSampAroundSpike) ); 
            N_retrieve = nSampAroundSpike*2+1;
            data_0 = readSamples(file_hnd_raw, dfInfo.channelIds, samp_start, N_retrieve);            
            
            data_1   = readSamples(file_hnd_filt(1), dfInfo.channelIds, samp_start, N_retrieve);            
            data_2   = readSamples(file_hnd_filt(2), dfInfo.channelIds, samp_start, N_retrieve);            
                        
            if 1
                %%
                tk_idx = [-nSampAroundSpike:nSampAroundSpike];
                figure(4); clf
                data_0_ms = bsxfun(@minus, data_0, mean(data_0, 2));
                [mx, idx_largest] = max( max(abs(data_0_ms), [], 2) );
                plot(tk_idx, data_0_ms(idx_largest, :), 'k.-', tk_idx, data_1(idx_largest, :), 'b.-', tk_idx, data_2(idx_largest, :), 'r.-');
%                 plot(tk_idx, data_0_ms(idx_largest, :), 'k.-');
                
                
            end
            
            
            
            3;
        end        
    end
    
    
    for fi = 1:nFilters
        fclose(file_info(fi).fid);    
    end
    
  %}

