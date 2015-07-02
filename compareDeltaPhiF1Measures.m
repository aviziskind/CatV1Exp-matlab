function compareDeltaPhiF1Measures
    
    plotCC_vsItpAndHarmonics = 0;
    showControlPanel         = 0;
    plotInterpolatedPeaks    = 0;
    barPlotDistMeansStd      = 1;

    

    gratingType = curGratingType('');
%     [pt_ids, pairTypes] = curPairTypes;
%     pairTypes_str = [pairTypes{:}];
    pairTypes_str = 'Wcc';
    
    limToComplexPairs = false;
    nPhSelect = 4;
    
    
    cmpDatafile = [CatV1Path gratingType 'GratingComparisonData_DB_tuning_' pairTypes_str '.mat'];         
    S_cmp = load(cmpDatafile);
    pairData = S_cmp.pairData;

    idx = true(1,length(pairData));
    if limToComplexPairs 
        f1odc_max = 0.5;
        f1odc_min = 0;
        
        pair_minF1oDCs = [pairData.loc_minF1oDC_cmp];
        pair_maxF1oDCs = [pairData.loc_maxF1oDC_cmp];    
        idx = (pair_minF1oDCs >= f1odc_min) & (pair_maxF1oDCs <= f1odc_max);
        if strcmp(gratingType, 'flashed') && ~isempty(nPhSelect)
            nPh = [pairData.n_phases];         
%             idx = idx & (nPh==nPhSelect);
            idx = (nPh==nPhSelect);
        end    
    end
    
    idx_dphi = find(strcmp(S_cmp.measures, 'dphi'), 1);
    idx_dF1 = find(strcmp(S_cmp.measures, 'dF1'), 1);
    orig_dphis = S_cmp.allStatsC{idx_dphi}.val(idx);
    orig_dF1s = S_cmp.allStatsC{idx_dF1}.val(idx);
    orig_cc = corr(orig_dphis, orig_dF1s, 'rows', 'complete');

    pairData = pairData(idx);

    compDphiFile = [CatV1Path gratingType 'GratingDphiComparison_' pairTypes_str '.mat'];
    redo = false;
    ignoreFile = false;
    
    if plotInterpolatedPeaks   
        nHrm = 4;
        tcs = cat(3, pairData.phaseTCs);
        f1odcs = cat(1, pairData.loc_F1oDCs_cmp);
        tc_peaks = zeros(size(f1odcs));
        nItp = 20;
        phs = linspace(0, 360, nPhSelect*nItp+1); phs = phs(1:end-1);
        progressBar('init-', size(tcs, 3));
        for i = 1:size(tcs, 3)
            progressBar;
            tc = tcs(:,1,i);
            [tmp, tmp2, tc_itp] = fourierInterp(tc, nHrm, nItp, 'spline');
            tc_peaks(i) = phs(indmax(tc_itp));
        end
        figure(8);
%         unique(
        
        plot(tc_peaks, f1odcs, 'o', 'markersize', 3)
        xlim([-1 361]); set(gca, 'xtick', linspace(0, 360, nPhSelect+1));
        xlabel('peak of interpolated smoothing curve');
        ylabel('F1/DC')
        nHrmTxt = iff(isempty(nHrm), 'all harmonics', sprintf('%d harmonics', nHrm) );
        title(sprintf('n = %d phases. %s', nPhSelect, nHrmTxt));
        
        figure(9);
        
            n_dPhiBins = 33;
            L = (360/(n_dPhiBins-1))/2;
            xBins_dphi = linspace(0-L, 360+L, n_dPhiBins+1);
        n = histcnt(tc_peaks, xBins_dphi);
        bar(binEdge2cent(xBins_dphi), n, 1);
        set(gca, 'xtick', linspace(0, 360, nPhSelect+1));
        xlim([xBins_dphi(1), xBins_dphi(end)])
        xlabel('peak of interpolated smoothing curve');
        title(sprintf('n = %d phases. %s', nPhSelect, nHrmTxt));
        3;
    end
    
    
    if strcmp(gratingType, 'flashed')        
%         interps = [1, 2, 5, 10, 100];
%         harmonics = [1 2 3 5 10];
%         interps = [1, 2, 3, 5, 10 30];

        interps = [1 2 5 10 200];
        harmonics = [1 2 3 4];

    elseif strcmp(gratingType, 'drifting')
%         interps = [1 2 3 4 5 20];
%         harmonics = [50];
        interps = [1 10];
        harmonics = [1 2 3 4 5 30];
        
    end
        
    nInterps = length(interps);
    nHarmonics = length(harmonics);    

    nPairs = length(pairData);
%     nPairs = 100;

    
    % CALCULATE OR LOAD

    if exist(compDphiFile, 'file') && ~redo && ~ignoreFile
        S = load(compDphiFile);
    end
    
    if exist('S', 'var') && isequal(S.interps, interps) && isequal(S.harmonics, harmonics) && isequal(S.nPairs, nPairs) 
        allDphs = S.allDphs;

    else
                    
        allDphs = zeros(nInterps, nHarmonics, nPairs);        
        tic;
        progressBar('init-', nPairs, 30);
        for p_i = 1:nPairs
            tc1 = pairData(p_i).phaseTCs(:,1,1);
            tc2 = pairData(p_i).phaseTCs(:,1,2);

            for itp_i = 1:nInterps   
                for hrm_i = 1:nHarmonics
                    allDphs(itp_i, hrm_i, p_i) = abs( deltaPhiStar(tc1, tc2, interps(itp_i), harmonics(hrm_i)) );
                end
            end
            progressBar;
        end
        progressBar('done');
        toc;
        if ~ignoreFile
            save(compDphiFile, 'interps', 'harmonics', 'nPairs', 'allDphs');
        end
        
    end

    
    if barPlotDistMeansStd    
        3;
        dPhs_mean = nanmean(allDphs, 3);
        dPhs_std = nanstderr(allDphs, 3);
        
        figure(745); clf;
%         h = bar(dPhs_mean);
        h1 = errorb(dPhs_mean', dPhs_std');        
        set(h1.bars, 'basevalue', 90)
        xlabel('# harmonics')
        legend(legendarray('nItp = ', interps), 'location', 'bestoutside');
        set(h1.bars, 'linewidth', 1)
        ylabel('mean \Delta\phi')        
        set(gca, 'xtick', 1:nHarmonics, 'xticklabel', num2str(harmonics'))
        
        for i = 1:nInterps
            for j = 1:nHarmonics
                [h, pvals(i,j)] = ttest(squeeze(allDphs(i,j,:)), 90);
            end
        end
        figure(746); clf;
        h = bar(-log10(pvals'));
        xlabel('# harmonics')
        set(gca, 'xtick', 1:nHarmonics, 'xticklabel', num2str(harmonics'))
        legend(legendarray('nItp = ', interps), 'location', 'bestoutside');
        ylabel('-log_{10} p-value (ttest)')
        
        
    end
    
    if plotCC_vsItpAndHarmonics 
        graphMode = 'plot';
        
        dim = iff((nInterps > 1) && (nHarmonics > 1), 3, 2);
        idx_Dphi = sub2ind([nInterps, nHarmonics], 1, nHarmonics);
        idx_dF1  = sub2ind([nInterps, nHarmonics], nInterps, 1);

        idx_selects = [idx_dF1, idx_Dphi, ];%(nHarmonics-1)*nInterps+1;
        select_x = [1, nHarmonics];    
        lbls = {'dF1', '\Delta\phi'};
        
%     sDp = 'r^';
%     sdp = 'ro';
%     sDF1 = 'b^';
%     sdF1 = 'bo';

    
        for i = 1:length(idx_selects)
            idx_select = idx_selects(i);
            figure(120+i); clf; hold on;
            allDphs_cols = reshape(permute(allDphs, [3, 1, 2]), [nPairs, nInterps*nHarmonics]);

            allCorrs_cols = corr(allDphs_cols(:,idx_select), allDphs_cols, 'rows', 'complete');
            allCorrs = reshape(allCorrs_cols, [nInterps, nHarmonics]);

            if (dim == 3) && strcmp(graphMode, 'bar')

                rng = [min(allCorrs(:)), max(allCorrs(:))];
                rng(1) = roundToNearest(rng(1), .25, 'down');

                bar3(allCorrs - rng(1));
                zlim([rng - rng(1)])
                zticks = get(gca, 'ztick');
                zticks = zticks + rng(1);
                set(gca, 'zticklabel', num2str(zticks'));

                xlabel('# harmonics');
                ylabel('# interps');
                set(gca, 'xtick', 1:nHarmonics, 'xticklabel', num2str(harmonics'))
                set(gca, 'ytick', 1:nInterps, 'yticklabel', num2str(interps'));

            elseif (dim == 2) || strcmp(graphMode, 'plot')

                if dim == 2
                    if nInterps > 1 
                        X = interps;
                        Xname = '# interp';
                    else
                        X = harmonics;
                        Xname = '# harmonics';
                    end
                elseif dim == 3
                    X = harmonics;
                    Xname = '# harmonics';
                    Y = interps;
                    Yname = 'nInt';                                

                end

                plot(1:length(X), allCorrs', 's-'); hold on;                

                plot(select_x(i), allCorrs(idx_select), 'r*');%'bs', 'markerfacecolor', 'g')
                set(gca, 'xtick', 1:length(X), 'xticklabel', num2str(X'));
                ylim([0 1])
                xlim([.8 length(X)+.2]);

                xlabel(Xname);
        %         ylabel('cc')

                ylabel(['cc with ' lbls{i}]);
    %             if idx_select > 1
    %                 ylabel('cc with \DeltaFn')
    %             else
    %                 ylabel('cc with \DeltaF1')
    %             end
                drawHorizontalLine(orig_cc, 'linestyle', ':', 'color', 'r');

                if dim == 3
                    legend(legendarray([Yname ' = '], Y), 'location', 'SE');                
                end

            end

        end
        
    end
    
    
    
    
    if showControlPanel
        
        msVsMsFig = 50;
        [int_grid, harm_grid] = meshgrid(interps, harmonics);
        wrp = @(x) [x(1:end); x(1)];
        ext = @(x) [x(1:end); x(end)+diff(x(1:2))];
        shft = @(x,n) x([end-n+1:end, 1:end-n]);
        
        [nHarm, nItp, dPhi_names, showTC1] = deal(0);
        
        allDPhiNames = arrayfun(@(hrm, int) sprintf('dF%d_itp%d', hrm, int), harm_grid, int_grid, 'un', 0);
        allDPhiNames = allDPhiNames(:);    

        figure(msVsMsFig); clf; 
        h_dphi12_L = plot(0,0, 'bo', 'markersize', 4);
%         h_dphi12_ax = gca;
        h_dphi12_xlab = xlabel(' ');
        h_dphi12_ylab = ylabel(' ');
        h_dphi12_tit = title(' ');
        d = 2;
        xlim([-d 180+d]);
        ylim([-d 180+d]);
        set(gca, 'xtick', [0:45:180], 'ytick', [0:45:180]);
        set(msVsMsFig, 'windowButtonDownFcn', {@selectPointsInFigure, @updateTuningCurvePlots});

        n_dPhiBins = 13;
        binE = linspace(0, 180, n_dPhiBins+1);
        binC = binEdge2cent(binE);
        
        figure(msVsMsFig+1); clf;         
        for i = 1:2
        
            h_tc_ax{i} = subplot(2,3, [3*i-2 : 3*i-1] ); %#ok<AGROW>
            h_tc{i} = plot(0,0, 'bo-', 0,0, 'b.:',  0,0, 'go-', 0,0, 'g.:', 0,0, 'bs'); %#ok<AGROW>
            h_tc_tit{i} = title(' '); %#ok<AGROW>
            ylim([-0.1,  1.1]);
            
            set(h_tc{i}(1), 'linewidth', 1, 'markerfacecolor', 'b')
            set(h_tc{i}(3), 'linewidth', 1, 'markerfacecolor', 'g')
            set(h_tc{i}(5), 'linewidth', 1, 'markersize', 3);
            
            h_cc_ax{i} = subplot(2,3, 3*i ); %#ok<AGROW>
            h_cc{i} = plot(0,0, 'k.-', 0,0, 'rs'); %#ok<AGROW>
            h_cc_tit{i} = title(' '); %#ok<AGROW>
            ylim([-1 1]);
            %         set(h_tc(6:7), 'markerfacecolor', 'auto');        
        end                
        figure(msVsMsFig+2); clf;         
        for i = 1:2        
            h_hist_ax{i} = subplot(2,1,i); %#ok<AGROW>
            h_hists{i} = bar(0,0,1); %#ok<AGROW>
            h_hists_tit{i} = title(' '); %#ok<AGROW>
        end
        
        xticks = [0:90:360];
        set([h_tc_ax{:}, h_cc_ax{:}], 'xlim', [0 360], 'xtick', xticks);
        set([h_cc_ax{:}], 'xticklabel', num2str(min(xticks, 360-xticks)'))
    
        updateComparisonPlot(allDPhiNames{1}, allDPhiNames{2}, false)

        args = { {'dPhi1', allDPhiNames, allDPhiNames{1}}, {'dPhi2', allDPhiNames, allDPhiNames{2}}, {'showTC1', [false, true]} };
        manipulate(@updateComparisonPlot, args, 'FigId', 149);

    end
    
        
%     zlim( rng + [-1, 1]*diff(rng)/10 );

    
    function updateComparisonPlot(dPhi1_name, dPhi2_name, showTC1_arg)
        
        dPhi_names = {dPhi1_name, dPhi2_name};
        showTC1 = showTC1_arg;
        
        [nHarm(1), nItp(1)] = dealV( sscanf(dPhi1_name, 'dF%d_itp%d') );
        [nHarm(2), nItp(2)] = dealV( sscanf(dPhi2_name, 'dF%d_itp%d') );
                
        nHarm_idx(1) = find(harmonics == nHarm(1), 1);
        nHarm_idx(2) = find(harmonics == nHarm(2), 1);
        nItp_idx(1) = find(interps == nItp(1), 1);
        nItp_idx(2) = find(interps == nItp(2), 1);
        
        X = allDphs(nItp_idx(1), nHarm_idx(1), :);
        Y = allDphs(nItp_idx(2), nHarm_idx(2), :);
        
        corr_r = corr(X(:), Y(:), 'type', 'pearson', 'rows', 'complete');
        corr_rho = corr(X(:), Y(:), 'type', 'spearman', 'rows', 'complete');                
        
        if isempty(h_hists_tit)
            h_hists_tit{1} = get(h_hist_ax{1}, 'title');
            h_hists_tit{2} = get(h_hist_ax{2}, 'title');
        end
        
        set(h_dphi12_L, 'xdata', X, 'ydata', Y);
        set([h_dphi12_xlab, h_hists_tit{1}], 'string', dPhi1_name, 'interpreter', 'none');
        set([h_dphi12_ylab, h_hists_tit{2}], 'string', dPhi2_name, 'interpreter', 'none');
        set(h_dphi12_tit, 'string', sprintf('r = %.3f, \\rho = %.3f (N = %d)', corr_r, corr_rho, length(X)));
                    
        nX = histcnt(X, binE);
        nY = histcnt(Y, binE);
        set(h_hists{1}, 'xdata', binC, 'ydata', nX);
        set(h_hists{2}, 'xdata', binC, 'ydata', nY);
        set([h_hist_ax{:}], 'xlim', [binE(1), binE(end)], 'xtick', [0:45:180]);
    end


    function updateTuningCurvePlots(glob_id, grp_id, loc_id, pt_dphiX, pt_dphiY) %#ok<INUSL>
                
        pt_vals = [pt_dphiX, pt_dphiY];
        pData = pairData( loc_id );        
%         assert(pt_dphi == S{idx_dphi}.val(pair_idx));
%         assert(pt_dF1  == S{idx_dF1}.val(pair_idx));

%         Gid = pData.Gids(1);
%         cellIds = pData.cellIds;
        
        tc1 = pData.phaseTCs(:,1,1);
        tc2 = pData.phaseTCs(:,1,2);
        nPh = length(tc1);
        ph_ext = linspace(0,360, nPh+1);  ph = ph_ext(1:nPh);
%         dph = round(360/nPh);
                      
        scl1 = 1/max(tc1);
        scl2 = 1/max(tc2);

        doTC1shift = true;
                
        for dph_i = 1:2

            [deltaPhi_sgn, phs{dph_i}, full_tc1{dph_i}, full_tc2{dph_i}] = deltaPhiStar(tc1, tc2, nItp(dph_i), nHarm(dph_i)); %#ok<AGROW>            
                        
            dphi{dph_i} = abs(deltaPhi_sgn); %#ok<AGROW>
            assert( dphi{dph_i} == pt_vals(dph_i) );
            
            ph_step = diff(phs{dph_i}(1:2));
            nPh_itp = 360/ph_step;
%             dphi_bin = deltaPhi(phs{dph_i}, full_tc1{dph_i}, full_tc2{dph_i}, 'cross-correlation') / ph_step;
            dphi_deg = deltaPhi(phs{dph_i}, full_tc1{dph_i}, full_tc2{dph_i}, 'cross-correlation');
            dphi_deg = mod(dphi_deg, 360);
            
            if doTC1shift
                phFull_shft = mod(phs{dph_i} + dphi_deg, 360);
                ph_shft = mod(ph + dphi_deg, 360);
            else
                phFull_shft = phs{dph_i};
                ph_shft = ph;
            end
            idx_shft_full = ord(phFull_shft);
            [ph_shft_srt, idx_shft] = sort(ph_shft);
            
            allDphiCCs = arrayfun(@(n) corr(shft(full_tc1{dph_i}, n), full_tc2{dph_i}), 1:nPh_itp);
            [max_cc, idx_max] = max(allDphiCCs);
                        
            set(h_cc{dph_i}(1), 'xdata', phs{dph_i}, 'ydata', allDphiCCs);
            set(h_cc{dph_i}(2), 'xdata', phs{dph_i}(idx_max), 'ydata', max_cc);
            
            
%             
%             tc1_shifted = shft(full_tc1{dph_i}, dphi_bin);

    %         [phi1_rad] = getF1phase( ph, tc1, 360);
    %         [phi2_rad] = getF1phase( ph, tc2, 360);
    %         phi1_deg = rad2deg(phi1_rad); phi1_y = f_cos1(indmin(abs(t_f1-phi1_deg)));
    %         phi2_deg = rad2deg(phi2_rad); phi2_y = f_cos2(indmin(abs(t_f2-phi2_deg)));                

%             f1_show = 'off';
            set(h_tc{dph_i}(1), 'xdata', ph_shft_srt, 'ydata', tc1(idx_shft)*scl1, 'linestyle', 'none', 'linewidth', 1);
            set(h_tc{dph_i}(3), 'xdata', ph_ext, 'ydata', wrp(tc2)*scl2, 'linestyle', 'none', 'linewidth', 1);        

            set(h_tc{dph_i}(2), 'xdata', ext(phs{dph_i}(:)), 'ydata', wrp(full_tc1{dph_i}(idx_shft_full))*scl1, 'linestyle', '-');
            set(h_tc{dph_i}(4), 'xdata', ext(phs{dph_i}(:)), 'ydata', wrp(full_tc2{dph_i})*scl2, 'linestyle', '-');
            
            set(h_tc{dph_i}(5), 'xdata', phs{dph_i}, 'ydata', full_tc1{dph_i}*scl1, 'linestyle', ':', 'marker', 'none', 'linewidth', 1, 'visible', iff(showTC1, 'on', 'off'));
            
%             y_ax2 = 
            ydata_pts = get(get(h_tc_ax{dph_i}, 'children'), 'ydata');
            ydata{dph_i} = [ydata_pts{:}]; %#ok<AGROW>
            
            nm = strrep(dPhi_names{dph_i}, '_', '\_');
            s = sprintf('%s.   \\Delta\\phi = %.1f', nm, deltaPhi_sgn);
            set(h_tc_tit{dph_i}, 'string', s);
        end
        ylim1 = min([ydata{:}, 0])-.1; 
        ylim2 = max([ydata{:}, 1])+.1;
        set([h_tc_ax{:}], 'ylim', [ylim1 ylim2]);
        
%         set(h_tc(3), 'xdata', ph_ext, 'ydata', wrp(tc2));        
%         set(h_tc(2), 'xdata', ph_ext, 'ydata', wrp(tc1_shft), 'visible', 'on');
%         set(h_tc(4), 'xdata', t_f1, 'ydata', f_cos1, 'visible', f1_show);
%         set(h_tc(5), 'xdata', t_f2, 'ydata', f_cos2, 'visible', f1_show);
%         set(h_tc(6), 'xdata', phi1_deg, 'ydata', phi1_y, 'visible', f1_show);
%         set(h_tc(7), 'xdata', phi2_deg, 'ydata', phi2_y, 'visible', f1_show);
%         set(h_tc(6), 'markerfacecolor', 'auto', 'marker', '*', 'color', 'r');%[.2 .2 1]);        
%         set(h_tc(7), 'markerfacecolor', 'auto', 'marker', '*', 'color', 'r');%[.2 1 .2]);        
%         ax = get(h_ms12_tit, 'parent');
%         legend off;
%         legend(ax, h_tc([1,3]), {['Cell ' num2str(cellIds(1))], ['Cell ' num2str(cellIds(2))]}, 'location', 'best')
%         tc_cc = pearsonR(tc1, tc2);
%         tc_rho = spearmanRho(tc1, tc2);
%         tc_dphi = deltaPhi(ph, tc1, tc2, 'cross-correlation');
%         tc_dF1 = deltaPhi(ph, tc1, tc2, 'angle');
%         [cc, rho, dp, df1] = deal(S{1}.val(pair_idx), S{2}.val(pair_idx), S{3}.val(pair_idx), S{4}.val(pair_idx) );
%         s1 = sprintf('Site # %d. Cells %d, %d.', Gid, cellIds(1), cellIds(2));        
%         s2 = sprintf('r = %.2f, \\rho = %.2f, \\Delta\\phi = %.0f, \\DeltaF1 = %.0f', cc, rho, dp, df1);
%         set(h_ms12_tit, 'string', {s1, s2}, 'fontsize', 10);%[.2 1 .2]);        
        
        
    end


    
end





function [dphi, phis, r1, r2] = deltaPhiStar(r1_orig, r2_orig, nInterp, nHarmonics)
    itpMethod = 'harmonics';
%     itpMethod = 'linear';
    
    showWorking = false;
    n = length(r1_orig);
    phis = linspace(0, 360, n*nInterp+1); phis = phis(1:end-1);
        
%     if nInterp > 1    
        [r1_hrm, phs_tmp, r1] = fourierInterp(r1_orig, nHarmonics, nInterp, itpMethod);
        [r2_hrm, phs_tmp, r2] = fourierInterp(r2_orig, nHarmonics, nInterp, itpMethod);                
%     elseif nInterp == 1
%         [r1] = fourierInterp(r1_orig, nHarmonics);
%         [r2] = fourierInterp(r2_orig, nHarmonics);
%     end
    
    if showWorking
        figure;
        phis_orig = linspace(0, 360, n+1); %phis_orig = phis_orig(1:end-1);
        wrp = @(x) [x(:); x(1)];
        figure(15); clf; plot(phis_orig, wrp(r1_orig), 'bs-', phis_orig, wrp(r2_orig), 'rs-',...
                              phis, r1, 'b:', phis, r2, 'r:');
    end
        
    dphi = deltaPhi(phis, r1, r2, 'cross-correlation');            
end


