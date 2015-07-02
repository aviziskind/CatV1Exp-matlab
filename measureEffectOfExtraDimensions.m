function measureEffectOfExtraDimensions

    M1 = -1;
    M2 = 1;
    
    N1 = 1000;
    N2 = 1500;
    

    randn('state', 1);
    showClusts = 0;

    
    all_sig1s = [.3:.1:.6];
    nSig1s = length(all_sig1s);
    all_nExtraDims = [0:30];
    nNExtras = length(all_nExtraDims);
    

    nTrials = 50;
    
    
    allErrs = nan(nSig1s, nNExtras, nTrials);
    [overlaps_pred, overlaps_pred_norm, overlaps_ac, overlaps_ac_norm] = deal( nan(nSig1s, nNExtras, nTrials) );    
    
    progressBar('init-', nSig1s * nNExtras * nTrials);

    
    for ti = 1:nTrials         
        
        rand('state', ti);
        X1_base = randn(N1, 1);
        X2_base = randn(N2, 1);
        
        for i_s1 = 1:nSig1s
            
            sig1 = all_sig1s(i_s1);
            sig2 = .45;

            X1 = sig1* X1_base +M1;
            X2 = sig1* X2_base +M2;
            
            
            for di = 1:nNExtras
                nExtraDims = all_nExtraDims(di);

        
        %         sig1 = .5;
        %         sig2 = .2;


%                 nExtraDim = 1;
                X1_full = [X1, randn(N1, nExtraDims)*sig2];
                X2_full = [X2, randn(N1, nExtraDims)*sig2];

%                 X_dim1 = [X1; X2];
                
%                 X = [X_dim1, randn(N1+N2, nExtraDims)*sig2];
                X = [X1_full; X2_full];
                
                % Calculate Predicted overlap
                m1_pred = [M1; zeros(nExtraDims, 1)];
                m2_pred = [M2; zeros(nExtraDims, 1)];
                c1_pred = diag([sig1^2; ones(nExtraDims, 1)*(sig2^2)]);
                c2_pred = c1_pred;
                
                ovlp_pred = quadProdGaussians(m1_pred, c1_pred, m2_pred, c2_pred, 'log');
                ovlp_pred_factor = quadProdGaussians(m1_pred, c1_pred, m1_pred, c2_pred, 'log');
                
                overlaps_pred(i_s1, di, ti) = ovlp_pred;
                overlaps_pred_norm(i_s1, di, ti) = ovlp_pred - ovlp_pred_factor;
                
                % Calculate Actual overlap
                m1_ac = mean(X1_full, 1);
                m2_ac = mean(X2_full, 1);
                c1_ac = cov(X1_full);
                c2_ac = cov(X2_full);
                                
                ovlp_ac = quadProdGaussians(m1_ac, c1_ac, m2_ac, c2_ac, 'log');
                ovlp_ac_factor = quadProdGaussians(m1_ac, c1_ac, m1_ac, c2_ac, 'log');
                
                overlaps_ac(i_s1, di, ti) = ovlp_ac;
                overlaps_ac_norm(i_s1, di, ti) = ovlp_ac - ovlp_ac_factor;
                
                

                id_orig = [ones(N1, 1); ones(N2,1)*2];

                %     ids_kk = klustakwik(X);

                gmm = gmdistribution.fit(X, 2, 'start', id_orig);
                ids_gmm = cluster(gmm, X);

                ids_gmm2 = renumberClusters(ids_gmm, id_orig);

                allErrs(i_s1, di, ti) = errRate(id_orig, ids_gmm2);

                progressBar;

                %%
                if showClusts

                    L = max(abs(lims([X1, X2], .05)));
                    nbin = 30;
                    binE = linspace(-L, L, nbin+1);
                    binC = binEdge2cent(binE);
                    binV_1 = histcnt(X1, binE);
                    binV_2 = histcnt(X2, binE);

                    figure(1); clf;
                    plot(binC, binV_1, 'b'); hold on;
                    plot(binC, binV_2, 'g'); 


                    idx1 = 1:N1;
                    idx2 = N1+[1:N2];

                    figure(2); clf; hold on;
                    plot(X(idx1, 1), X(idx1, 2), 'b.');
                    plot(X(idx2, 1), X(idx2, 2), 'g.');    

                    figure(3); clf; hold on;
                    [uId, clustIdxs] = uniqueList(ids_gmm2);
                    for i = 1:length(uId)
                        plot(X( clustIdxs{i}, 1), X( clustIdxs{i}, 2), [color_s(i) '.']);
                    end
                    title(sprintf('nClust = %d', length(unique(uId))));

                    figure(10); clf;
                    plot(all_sig1s, allErrs, 'o-');
                end


            end
        end
        3;
    end
    progressBar('done');
    3;

    %%
    3;
    figure(10); clf; hold on;
    allErrs_trialsAv = mean(allErrs, 3);
    plot(all_sig1s, allErrs_trialsAv, 'o-');
    if nTrials > 1
        all_sig1s_M = repmat(all_sig1s(:), [1, nNExtras]);
        errorbar(all_sig1s_M, allErrs_trialsAv, stderr(allErrs, 3), '.');
    end
    xlim(lims(all_sig1s, .05));
    
    legend(legendarray('nExtra = ', all_nExtraDims), 'location', 'NW');
    
    %%
    figure(11); clf; hold on
    allErrs_trialsAv_norm0Extra = bsxfun(@rdivide, allErrs_trialsAv, allErrs_trialsAv(:,1));
    plot(all_nExtraDims, allErrs_trialsAv_norm0Extra);
    rel_sigmas = all_sig1s / sig2;
    legend(legendarray('\sigma1/\sigma2 = ', rel_sigmas), 'location', 'NW')
    3;

    overlaps_pred_av = mean(overlaps_pred,3);
    overlaps_pred_norm_av = mean(overlaps_pred_norm,3);
    overlaps_ac_av = mean(overlaps_ac,3);
    overlaps_ac_norm_av = mean(overlaps_ac_norm,3);
    
    %%
    figure(12); clf; hold on
    overlaps_pred_norm_av_norm0Extra = bsxfun(@rdivide, overlaps_pred_norm_av, overlaps_pred_norm_av(:,1));
    overlaps_ac_norm_av_norm0Extra = bsxfun(@rdivide, overlaps_ac_norm_av, overlaps_ac_norm_av(:,1));
    plot(all_nExtraDims, overlaps_pred_norm_av_norm0Extra', ':');
    plot(all_nExtraDims, overlaps_ac_norm_av_norm0Extra', '.-');
    
    %%
    overlaps_ac_norm_norm0Extra = bsxfun(@rdivide, overlaps_ac_norm, overlaps_ac_norm(:,1,:));    
    all_nExtraDims_M = repmat(all_nExtraDims(:)', [nSig1s, 1]);
    figure(13); clf; hold on        
    errorbar(all_nExtraDims_M', mean(overlaps_ac_norm_norm0Extra, 3)', stderr(overlaps_ac_norm_norm0Extra, 3)', '.-');
    legend(legendarray('\sigma1/\sigma2 = ', rel_sigmas), 'location', 'NW')
    %     rel_sigmas = all_sig1s / sig2;
%     legend(legendarray('\sigma1/\sigma2 = ', rel_sigmas), 'location', 'NW')
    3;
    

end

function r = errRate(ids_orig, ids_calc, eps)
%     if ~exist('eps', 'var')
%         eps = 0.5;
%     end
    r = nnz(ids_orig ~= ids_calc)/length(ids_orig);
    
    
end

function ids_renum = renumberClusters(ids, ids_ref)
    %%
    Nclust = length(unique(ids_ref));
    combos = [ids(:), ids_ref(:)];
    [uRows, comboCount] = uniqueCount(combos, 'rows');
    
    A = zeros(Nclust, Nclust);
    for i = 1:size(uRows, 1)
        A(uRows(i,1), uRows(i,2)) = comboCount(i);
    end
    
    [~, idx_max] = max(A, [], 1);
    
    ids_renum = zeros(size(ids));
    for j = 1:Nclust
        ids_renum(ids == j) = idx_max(j);
    end
    
    3;


end
