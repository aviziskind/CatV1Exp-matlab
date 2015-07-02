function exploreBestNumberOfDimensions

    doIDcalc = 0;
    doIndivQuadCalc = 1;
    doCumQuadCalc = 0;


    featureType = 'GLF';
    channelAction = 'concat';
    normChannels = 0;
    whitenChannels = 1;
    
    maxNDim = 88;

    opt.doOtherCellSpk = 0;
    opt.doID_radius = 1;
    
    Gid = 6001;
    allCoeffs = getGroupWaveformCoefficients(featureType, Gid, maxNDim, channelAction, normChannels, whitenChannels);
    
    %%    
    [stats, idxs] = identifyCellfromIC(Gid, struct('clustGrouping', 'clusters'), [], 1);
    IC_idx = idxs.correctDetections_EC_spkIdx;
    figure(55); clf;
    plot(allCoeffs(:,1), allCoeffs(:,2), 'bo', 'markersize', 2); hold on;
    plot(allCoeffs(IC_idx,1), allCoeffs(IC_idx,2), 'ro', 'markersize', 2);
    figure(56); clf;
    plot(allCoeffs(:,3), allCoeffs(:,4), 'bo', 'markersize', 2); hold on;
    plot(allCoeffs(IC_idx,3), allCoeffs(IC_idx,4), 'ro', 'markersize', 2);

    3;
    %%
        
    nSpikes = size(allCoeffs, 1);
    
    if doIDcalc
    
        allIDs = zeros(1, maxNDim);
        allLs = zeros(1, maxNDim);
        allID_radii = zeros(1, maxNDim);
        progressBar('init-', maxNDim)

        for i = 1:maxNDim

            fet = allCoeffs(:,1:i);
            [L_ratio, L, ID, frac, ID_radius] = calculateIsolationDistance(fet, IC_idx);
            allIDs(i) = ID;
            allID_radii(i) = ID_radius;
            allLs(i) = L;

            progressBar(i)
        end

        3;
        figure(93); plot(allIDs ./allID_radii, '.-')
    %     profile viewer;

        figure(94); plot(allLs, '.-')
        3;
    end

    if doIndivQuadCalc
        
        allQuads = zeros(1, maxNDim);
        
        ks_stats = zeros(1, maxNDim);
        ks_pvals = zeros(1, maxNDim);
        allLinearErrs = zeros(1, maxNDim);
        noise_idx = 1:nSpikes; noise_idx(IC_idx) = [];
        for dim_i = 2:maxNDim
            clust_fet_i = allCoeffs(IC_idx,dim_i);
            noise_fet_i = allCoeffs(noise_idx,dim_i);
            M1 = mean(clust_fet_i);
            C1 = var(clust_fet_i);
            M2 = mean(noise_fet_i);
            C2 = var(noise_fet_i);
            
            [h, ks_pvals(dim_i), ks_stats(dim_i)] = kstest2(clust_fet_i, noise_fet_i);
            
            allQuads(dim_i) = quadProdGaussians(M1, C1, M2, C2, 'log') - quadProdGaussians(M1, C1, M1, C2, 'log');
            [th, bestErr] = optimalLinearThreshold(clust_fet_i, noise_fet_i, .5);
            allLinearErrs(dim_i) = bestErr;
            
        end

        figure(104); plot(allQuads, '.-')
        figure(105); plot(ks_stats, '.-');
        figure(106); plot(allLinearErrs, '.-');
    %     profile viewer;
        3;
        
        
        
        
    end

    if doCumQuadCalc
        allQuads = zeros(1, maxNDim);                
        
        noise_idx = 1:nSpikes; noise_idx(IC_idx) = [];
        for dim_i = 1:maxNDim
            clust_fet_i = allCoeffs(IC_idx,1:dim_i);
            noise_fet_i = allCoeffs(noise_idx,1:dim_i);
            M1 = mean(clust_fet_i, 1);
            C1 = cov(clust_fet_i);
            M2 = mean(noise_fet_i, 1);
            C2 = cov(noise_fet_i);
            
            allQuads(dim_i) = quadProdGaussians(M1, C1, M2, C2, 'log') - quadProdGaussians(M1, C1, M1, C2, 'log');
%             [th, err] = optimalLinearThreshold(clust_fet_i, noise_fet_i, .5);
            
        end

        figure(204); plot(allQuads, '.-')
    %     profile viewer;
        3;
    end



end