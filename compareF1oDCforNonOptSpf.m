function compareF1oDCforNonOptSpf
%%    
3;
    curGratingType(1);

    
    S_d = load('driftingGratingCells_GLFcuw8_degree');
    idx_ori = find(strncmp({S_d.allCells.stimType}, 'Grating:Orient', 9) & [S_d.allCells.cellId] > 0);
    
    dpp = arrayfun(@(s) s.degPerPix, S_d.allCells(idx_ori));
    sp  = arrayfun(@(s) s.sp, S_d.allCells(idx_ori));
    
    spfs_cpd_dg = 1./(sp .* dpp);        
    
    3;
    
    S = load('flashedGratingCells_GLFcuw8_degree');
    nph = arrayfun(@(s) length(s.ph), [S.allCells]);
    allCells = S.allCells(nph == 8);
        
    Ncl = length(allCells);
    
    allF1oDCs = zeros(Ncl, 10);
    optF1oDC  = zeros(Ncl, 1);
    F1oDCat05  = zeros(Ncl, 1);
    
    idx_closeTo05 = zeros(Ncl, 1);
    all_spfs_cpd = zeros(Ncl, 10);
    
    for cell_i = 1:Ncl
        
        cell_data = allCells(cell_i);
        R = cell_data.R;
        
        R_os_smoothed = mean(R,3);
        
        dw_ori = 1;         
        R_os_smoothed = gaussSmooth(R_os_smoothed, dw_ori, 1, 1);
        
        log_spf = log(cell_data.sp);
        dw_spf = mean(diff(log_spf));
        R_os_smoothed = gaussSmooth_nu(log_spf, R_os_smoothed, dw_spf, 2);    
        
        ori_pref_idx = indmax( sum(R_os_smoothed, 2) );
        spf_pref_idx = indmax( R_os_smoothed(ori_pref_idx(1), :) );

        % calculate F1/DC at opt.
        nPh = size(R, 3); phs = linspace(0, 360, nPh+1); phs = phs(1:nPh);        
        ori_pref_idx = indmax( sum(R_os_smoothed, 2) );
        for spf_j = 1:10
            phaseTC_i = R(ori_pref_idx, spf_j, :);
            allF1oDCs(cell_i,spf_j) = getF1oDC(phs, phaseTC_i(:), 360);                        
        end
        F1oDC_atSpfPref = allF1oDCs(cell_i, spf_pref_idx);
        assert( abs(F1oDC_atSpfPref - cell_data.stats.tuningStats.oriStats_si.F1oDC) < 1e-5 );        
        optF1oDC(cell_i) = F1oDC_atSpfPref;
                        
        spPeriod_pix = cell_data.sp;
        degPerPix = cell_data.degPerPix;
        spfs_cpd = 1./(spPeriod_pix * degPerPix);        
        
        all_spfs_cpd(cell_i,:) = spfs_cpd;
        
        idx_closeTo05(cell_i) = indmin( abs( log2(spfs_cpd) - log2(0.5)) );
        F1oDCat05(cell_i) = allF1oDCs(cell_i, idx_closeTo05(cell_i));
    end
    
    %%
    3;
    figure(55);
    plot(optF1oDC, F1oDCat05, 'o', 'markersize', 2);
    L = lims([0 2], .015);
    
    axis equal square;
    axis([L; L]);
    drawHorizontalLine(1);
    drawVerticalLine(1);
    xlabel('F1/DC at optimal Spatial Frequency'); ylabel('F1/DC at 0.5 cyc/deg')
    
%     for j = 1:
    %%
    [cc, err_frac, n_err] = getF1oDC_ErrorCC(optF1oDC, F1oDCat05);
    
    title(sprintf('cc: %.2f. S/C error: %d/%d (%.2g%%)', cc, n_err, Ncl, err_frac*100))
    %%
%     spfs_cpd_test = logspace( [0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2];
    
    spfs_cpd_test = logspace(log10(0.1), log10(2.3), 20);
    n_test = length(spfs_cpd_test);
    err = zeros(n_test, 1);
    ccs = zeros(n_test, 1);
    F1oDCs_test = zeros(Ncl, 1);
    
    for j = 1:n_test
        spf_cur = spfs_cpd_test(j);
        for ci = 1:Ncl
            idx = indmin( abs( log2(spfs_cpd) - log2(spf_cur)) );
            F1oDCs_test(ci) = allF1oDCs(ci, idx);
        end
        
        [ccs(j), err(j)] = getF1oDC_ErrorCC(optF1oDC, F1oDCs_test);
        
    end
    %%
    figure(56);
    plot(spfs_cpd_test, ccs, 'bo-', spfs_cpd_test, err, 'rs-');
    legend('cc', 'error')
    xlabel('spatial frequency (cyc/deg)');
    3;
    



end

function [cc, err_frac, n_err] = getF1oDC_ErrorCC(optF1oDC, testF1oDC)

    nCC = nnz(optF1oDC <  1 & testF1oDC <  1);
    nSC = nnz(optF1oDC >= 1 & testF1oDC <  1);
    nCS = nnz(optF1oDC <  1 & testF1oDC >= 1);
    nSS = nnz(optF1oDC >= 1 & testF1oDC >= 1);
    cc = corr(optF1oDC, testF1oDC, 'rows', 'complete');
    n_err = nSC+nCS;
    err_frac = n_err/length(optF1oDC);

end