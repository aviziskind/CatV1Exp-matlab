function slides5


    % examples of cells at a site

    use_sm = 1;
    
    doFlashed = 1;
    doDrifting_spf = 1;
    doDrifting_ori = 1;
    
    
    
    % f1/dc plots:
    %%
      
    if doFlashed
        S_f = load('flashedGratingCells_GLFcuw8_degree_ori');
        %%
        idx_cells = [S_f.allCells.cellId] > 0;
        allCells_f = S_f.allCells(idx_cells);

        binE = 0:.2:2;
        use_sm = 1;
        nPh_use = [8];
        cell_nPh = arrayfun(@(s) length(s.ph), allCells_f);
        idx_use = find(cell_nPh == nPh_use);
    %     idx_use = 1:length(cell_nPh);

        if ~use_sm
            f_F1oDC_all = [allCells_f.F1oDC_maxR_avP_sm];        
        else
            f_F1oDC_all = [allCells_f.F1oDC_maxR_avP];
        end
        f_F1oDC = f_F1oDC_all(idx_use);
        nSimple = nnz(f_F1oDC_all>=1);
        nComplex = nnz(f_F1oDC_all<1);
        nTotal = length(f_F1oDC_all);
        fprintf('Flashed: nSimple:%d. nComplex:%d. Ntotal = %d\n', nSimple, nComplex, nTotal)


        plotF1oDCs(f_F1oDC, binE, 11)
        title('F1/DC ratios for responses to flashed gratings');
    %%
    

    %     idx_use = 1:length(cell_nPh);

        ori_rep_f = nestedFields(allCells_f(idx_use), 'stats', 'tuningStats', 'oriStats_si', 'ori_rep_pval');    
        good_rep_f = ori_rep_f < 1e-10;

        idx_simp_f = ibetween(f_F1oDC, [1.3, 1.7]);
        idx_comp_f = ibetween(f_F1oDC, [0.3, 0.4]);
        idx_good_simple_f = find(good_rep_f & idx_simp_f);
        idx_good_complex_f = find(good_rep_f & idx_comp_f);

        i = 0;    
        %%

        % nice simple: idx=56, i=1, shift =0
        % nice complex: idx=143, shift =0
    %     figure(60);  

    %     clf;
        figure(61); clf;
    %     i = i+1;
        idx = idx_use(idx_good_simple_f(i));
    %     idx = idx_use(idx_good_complex_f(i));
        idx = 143;

        plotPTC(S_f, idx, 1, 0);
        ylims = ylim;
        ylim([0, ylims(2)]);
        3;
    %     oriSp_maxR_avP = cdata.oriSp_maxR_avP;    
    %     r = squeeze(cdata.R(oriSp_maxR_avP(1), oriSp_maxR_avP(2), :));
    %     PSTH = cdata.PSTH;
    %     LR_bins = PSTH.timeWindow_bins;
    %     R_full = getOspDataForPsthWindow(cdata.Gid, cdata.cellId, [], [], LR_bins(1), LR_bins(2), PSTH.windowProfile, 'osp_full');
    % %     R_os = mean(mean(R_full,3),4);
    % %     [~, idx_max] = maxElement(R_os);
    %     [nOri, nSp, nPh, nTrials] = size(R_full);
    %     R_ph_trials = reshape(R_full(oriSp_maxR_avP(1), oriSp_maxR_avP(2), :, :), [nPh, nTrials]);
    % 
    % %     r = gaussSmooth(r, 1.5);
    %     ph = cdata.ph;   
    %     plot(ph, R_ph_trials', '.-');    
    %     3;


    end
    
    if doDrifting_spf
        %%
        S_d_s = load('driftingGratingCells_GLFcuw8_degree_spf');    
        idx_cells = [S_d_s.allCells.cellId] > 0;
        allCells_ds = S_d_s.allCells(idx_cells);

        use_sm = 1;
        if ~use_sm
            ds_F1oDC = [allCells_ds.F1oDC_maxR_avP_sm];
        else
            ds_F1oDC = [allCells_ds.F1oDC_maxR_avP];
        end

        nSimple_d = nnz(ds_F1oDC>=1);
        nComplex_d = nnz(ds_F1oDC<1);
        nTotal_d = length(ds_F1oDC);
        fprintf('Drifting (spf): nSimple:%d. nComplex:%d. Ntotal = %d\n', nSimple_d, nComplex_d, nTotal_d)    

        %%
        binE = 0:.2:2;
        plotF1oDCs(ds_F1oDC, binE, 12);
        title('F1/DC ratios for responses to drifting gratings');
        3;
        %%
        spf_stats = nestedFields(allCells_ds, 'stats', 'tuningStats', 'spfStats_si');        
        spf_rep_pval = arrayfun(@(s) s.spf_rep_pval(1), spf_stats);
        good_rep = spf_rep_pval < 1e-5;

        idx_simp = ibetween(ds_F1oDC, [1.3, 1.7]);
        idx_good = find(good_rep & idx_simp);

        i = 0;
        %%
        % good simple: idx = 55;
        % good complex: idx = 16;

        figure(61); clf; 
        i = i+1;
        idx = idx_good(i);
    % idx = 55;
        plotPTC(S_d_s, idx, 1, 0);


        ylims = ylim;
        ylim([0, ylims(2)]);    
        3;
    %     oriSp_maxR_avP = allCells_ds(idx_good(i)).oriSp_maxR_avP;
    %     
    %     r = squeeze(allCells_ds(idx_good(i)).R(oriSp_maxR_avP(1), oriSp_maxR_avP(2), :));    
    %     r = gaussSmooth(r, 1.5);
    %     ph = allCells_ds(idx_good(i)).ph;   
    %     plot(ph, r, '.-');



    
    end
    
    if doDrifting_ori
    
        %%
        S_d_o = load('driftingGratingCells_GLFcuw8_degree_ori');    
        idx_cells = [S_d_o.allCells.cellId] > 0;
        allCells_do = S_d_o.allCells(idx_cells);

        use_sm = 1;
        if ~use_sm
            do_F1oDC = [allCells_do.F1oDC_maxR_avP_sm];
        else
            do_F1oDC = [allCells_do.F1oDC_maxR_avP];
        end

        nSimple_d = nnz(do_F1oDC>=1);
        nComplex_d = nnz(do_F1oDC<1);
        nTotal_d = length(do_F1oDC);
        fprintf('Drifting (ori): nSimple:%d. nComplex:%d. Ntotal = %d\n', nSimple_d, nComplex_d, nTotal_d)    

        %%
        binE = 0:.2:2;
        plotF1oDCs(do_F1oDC, binE, 13);
        title('F1/DC ratios for responses to drifting gratings');
    
    end    
end

function plotF1oDCs(F1oDCs, binE, figId)
    idx_simple = F1oDCs > 1;        
    binC = binEdge2cent(binE);    
    
    binV_simp = histcnt(F1oDCs(idx_simple), binE); 
    binV_comp = histcnt(F1oDCs(~idx_simple), binE); 
    figure(figId);
    h = bar(binC, [binV_comp(:), binV_simp(:)], 1, 'stacked');    
    set(h(1), 'facecolor', 'r');
    set(h(2), 'facecolor', 'b');
    drawVerticalLine(1, 'linestyle', ':');
    xlabel('F1/DC'); ylabel('Number of Cells');
    legend('Complex Cells', 'Simple Cells');

end

function plotPTC(S, idx, std_flag, nshift, plotMeanR)

    if nargin < 3
        std_flag = 1;
    end
    if nargin < 4
        nshift = 0;
    end
    if nargin < 5
        plotMeanR = 0;
    end

    cdata = S.allCells(idx);
    Gid = cdata.Gid;
    cellId = cdata.cellId;
%     idx = find(([S.allCells.Gid] == Gid)& ([S.allCells.cellId] == cellId),1);
    ph = cdata.ph;
    gratingType = flashedOrDrifting(Gid);
    if gratingType == 1
                
        PSTH = cdata.PSTH;
        LR_bins = PSTH.timeWindow_bins;        
        %     R_os = mean(mean(R_full,3),4);
        %     [~, idx_max] = maxElement(R_os);
        l_bin = LR_bins(1);
        r_bin = LR_bins(2);
        windowProfile = PSTH.windowProfile;
        sm = 0;
    else
        l_bin = 1;
        r_bin = 1;
        windowProfile = [];
        sm = 1;
    end
    oriSp_maxR_avP = cdata.oriSp_maxR_avP;
    R = cdata.R;
    [R_full, meanRate] = getOspDataForPsthWindow(Gid, cellId, [], [], l_bin, r_bin, windowProfile, {'osp_full', 'meanRate'});
    [nOri, nSp, nPh, nTrials] = size(R_full);
    R = R/mean(R(:))*meanRate;
    R_full = R_full/mean(R_full(:))*meanRate;    

    if nshift > 0
        idx_shift = [nPh-nshift+1:nPh, 1:nPh-nshift];    
        R = R_full(:,:,idx_shift);
        R_full = R_full(:,:,idx_shift,:);
    end    
    
    r = squeeze(cdata.R(oriSp_maxR_avP(1), oriSp_maxR_avP(2), :));    
    
    R_ph_trials = reshape(R_full(oriSp_maxR_avP(1), oriSp_maxR_avP(2), :, :), [nPh, nTrials]);
    
    figure(63);
    subplot(1,2,1); imagesc(mean(mean(R_full,3),4)); colorbar;
    subplot(1,2,2); imagesc(mean(R,3)); colorbar;
    
        
    
    
    if sm
        R_ph_trials = gaussSmooth(R_ph_trials, 1.5, 1, 1);
        r = gaussSmooth(r(:), 1.5, 1, 1);
    end
            
    ph_ext = [ph(:)', 360];
    R_ph_trials_ext = [R_ph_trials; R_ph_trials(1,:)];
    r_ext = [r(:);r(1)];
        %%
    figure(60); clf;
    plot(ph_ext, R_ph_trials_ext', '.-');       
    figure(61); clf; hold on; box on;
    y_mean = mean(R_ph_trials_ext,2);
    if nargin < 1 || std_flag == 1
        y_std = stderr(R_ph_trials_ext, 2);
    elseif std_flag == 2
        y_std = std(R_ph_trials_ext, [], 2);
    end

    doErrorBars = gratingType == 1;
    
    if doErrorBars
        errorbar(ph_ext, y_mean, y_std, 'bo-'); 
    else
        plot(ph_ext, y_mean, '.-'); hold on;
        plot(ph_ext, y_mean+y_std, ':');
        plot(ph_ext, y_mean-y_std, ':');
    end
    set(gca, 'xtick', [0:90:360]);
    if gratingType == 1
        xlabel('Spatial Phase');
    else
        xlabel('Temporal Phase');
    end
    ylabel('Firing Rate (Hz)');
    f1odc = getF1oDC(ph, y_mean(1:nPh));
    title(sprintf('F1/DC = %.1f', f1odc));
    
    xlim(lims([0, 360], .02));
    if plotMeanR
        hold on;
        plot(ph_ext, r_ext, 'rs-');
    end
    
    
%     if sm == 1
%         %%
%         gt = getGratingStimType(Gid);
%         nUniqueSeq = gt.nUniqueSeq;
%         nRep = gt.nSeqRep;
%         
%         R_full5 = reshape(R_full, [nOri, nSp, nPh, nUniqueSeq, nRep]);
%         
%         
%     end
    
    3;
    
end