function simulateDiffOriPairs

    
%     runmodes = {'confint', 'plot_v_sigma', 'diff_hist'};
%     runmode = 'confint';
%     runmode = 'plot_v_sigma';
    runmode = 'diff_hist';

%     median_dOriMU_observed = 8.7;
%     median_dOriMU_observed = 9.2;
    saveSampleOriMapFig = 1;  % A
    saveDiffHistFig = 0;  % B
    saveMedOriPctOutlierVsSigmaFigs = 1; % C,D
    saveDistSigmaOutlierFigs = 0; % E,F
    saveMedianVsDistFigs = 0; % G, H

    
     nCentPoints_side_set = 4; 
    doLargeMapForDiffHist = 1;
    if doLargeMapForDiffHist && any(strcmp(runmode, {'diff_hist', 'plot_v_sigma'}))
%         nCentPoints_side_set = 15;
        nCentPoints_side_set = 4;
    end

    median_dOriMU_observed = 8.72;
    pctOutliers_observed = 3.7;
    
    
%     curEstimateOfGaussSigma = 92.3;
    curEstimateOfGaussSigma = 84.51;
    
    % step 1: run with runmode='plot_v_sigma', to get a *rough* estimate of
        % the correct value of sigma.
    % step 2: run with runmode='confint', (set the appropriate values in
        % dist_vals and do 500 reps to get a good estimate of sigma.
        % and to find the index of the rep that gave the closest to the
        % mean.  [--> get plots with conf intervals]
        
        idx_bestRep = 497;
        
    % step 3: run with runmode='plot_v_sigma', with the rep_idx set to the
        % rep_idx found in step 2. [--> get both plots v sigma]
        % also, get map_idx that corresponds to best match to median
        % observed.
%         idx_bestMap = 8; % ncent_side = 4;
        idx_bestMap = 17;  % ncent_side = 15;
        
    % step 4: run with runmode='diff_hist' with the correct rep_idx and
        % map_idx to get the difference histogram.
     
        
    

    % in drifting grating ori batches, used 265 recording sites, 17 animals
    

    nGlobalRep = 500;
    nMaps = 17;
    
    um_spacing = 10; % for plotting nice graph

%     nCentPoints_side_set = 10; 
   
    
    rangeIntegrate_um = 600;
    rangeSpacing_um = 1000;

    
    
%     n_wav = 8;
    lambda_mm = 1.2;
    
    

    
    nModesForOriMaps = 50;
    
    
    
    make3Dmap = 0;
%   
%     r50_vals = [10 : 5 : 100, 110:10: 500 ];
%     r50_vals = [5 : 5 : 130]; 
% %         r50_vals = [5 : 5 : 130]; 
%     r50_vals = [5, 10]; 
%     dist_vals = [80, 100, 110, 114.3, 120, 150];
    figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];    
    
    figureSize = [370, 250];
    figureLB = [1300, 500];
    figureLB_offset = [20,-20];
    
   

    
     median_glob_estimates = zeros(1, nGlobalRep);
    pct_outlier_estimates = zeros(1, nGlobalRep);
    
%     nCentPoints_side = 15;
%     nCentPoints = nCentPoints_side^2;
    
    
%     seeing_distance_func = 'exponential';
    seeing_distance_func = 'gaussian';

    exp_r50_factors   = [0.6931, 1.6785,  2.6741];  % for D = 1, D = 2, D = 3.
    gauss_r50_factors = [0.6745, 1.17741, 1.538172];
    gauss_sig_to_R50_3D_factor = 1.538172;
    
    useR50inCalc = 0;
    if useR50inCalc
        dist_vals = [80, 100, 110, 114.3, 120, 150];
    else
        switch seeing_distance_func
%             case 'gaussian',    dist_vals = [50, 60, 70, 75, 80, 90];
%             case 'gaussian',    dist_vals = [60, 75, 90, 100];
                case 'gaussian',    
                    switch runmode, 
                        case 'plot_v_sigma', dist_vals = [60:5:110];  % (2) multiple distances - for plot vs sigma.
                        case 'confint',      dist_vals = curEstimateOfGaussSigma + [-25, 0, 25];  % (1) for multiple estimates, to generate confidence intervals 
                        case 'diff_hist',    dist_vals = [curEstimateOfGaussSigma];  % (3) close to actual value - for histogram of pairwise differences
                    end
%                     case 'gaussian',    dist_vals = [60];
%                 case 'gaussian',    dist_vals = [55:10:85];
            case 'exponential', dist_vals = [30:5:50];
        end
        
    end
%     um_spacing = 10; % 50?    
%     um_spacing = 30;

  
    
    
    
    nDists = length(dist_vals);
    
%     fracFromEdge = 1/5;    
    showMaps = saveSampleOriMapFig && 1;
    
    nDims = iff(make3Dmap, 3, 2);

    
    

    
    lambda_um = lambda_mm * 1000;
    margin = 0;
    %         point_spacing_microns = 50;
%     um_spacing = 50; % 50?
%     points_per_mm = 1000/um_spacing;
    column_height_um = 500;
%     nPointsPerWavelength = ceil( points_per_mm / (lambda/1000) );
    nPointsPerWavelength = ceil( lambda_um / um_spacing );
    

    
    
    n_wav = (nCentPoints_side_set+1)*rangeSpacing_um*(1.005)  / (lambda_um);


    
    rng_cent = [rangeSpacing_um, n_wav*lambda_um-rangeSpacing_um];

    xs_cent = rng_cent(1) : rangeSpacing_um : rng_cent(2);
    ys_cent = xs_cent;
    nCentPoints_side = length(xs_cent);
    nCentPoints = nCentPoints_side^2;
    
    
    median_values = zeros( nDists, nMaps);
%     median_values_eachCentCell = zeros(nR, nCentPoints, nMaps);
%     median_values_no3Dwgt = zeros(nDists, nMaps);
%     median_values_unconv = zeros(nDists, nMaps);
    
%     median_values_unconv_orig = zeros( nDists, nMaps);
    
    pct_outliers = zeros( nDists, nMaps);
%     pct_outliers_no3Dwgt = zeros( nDists, nMaps);
    
    dists_fromCent = cell(1,nCentPoints);
    idx_inRange_fromCent = cell(1,nCentPoints);
    diffOri_fromCent = cell(1,nCentPoints);
    
    firstTime = 1;
    rng('default');
    progressBar('init-', nGlobalRep * nMaps * nDists * nCentPoints);    

%     rand_seeds = reshape(1:nGlobalRep*nMaps, [nMaps, nGlobalRep])';
    
    switch runmode
        case 'confint',      rep_idxs = 1:nGlobalRep;
        case 'plot_v_sigma', rep_idxs = idx_bestRep;   % once you know the best value of sigma from running confint, use the rep_idx that gets the closest value (set here)
        case 'diff_hist',    rep_idxs = idx_bestRep; %177;
    end
    

    for rep_i = rep_idxs
        rng(rep_i);

        switch runmode
            case 'confint',      map_idxs = 1:nMaps;
            case 'plot_v_sigma', map_idxs = 1:nMaps;
            case 'diff_hist',    map_idxs = idx_bestMap;
        end
        
        if length(map_idxs) == 1 && strcmp(runmode, 'diff_hist');
            for map_i_tmp = 1:map_idxs-1
                [Z_map, x_map, y_map] = generateOrientationMap(nModesForOriMaps, n_wav, lambda_um, margin, nPointsPerWavelength);
            end
        end
        
        for map_i = map_idxs

%             rand_seed = rand_seeds(rep_i, map_i);
            
    %         dens_range = [3, 3.3];
            %% generate ori map
            plotInMM = 1;
            plot_scale_factor = iff(plotInMM, 1/1000, 1);

            [Z_map, x_map, y_map] = generateOrientationMap(nModesForOriMaps, n_wav, lambda_um, margin, nPointsPerWavelength);

            
%             [Z_map, x_map, y_map] = getOrientationMap(nModes, n_wav, lambda_um, nPointsPerWavelength, rand_seed);
            ori_map = angle(Z_map)/2 + pi/2;
%             ori_map = angle(Z_map)/2 + pi/2;
%             assert(isequal(Z_map, Z_map2))
%             assert(isequal(x_map, x_map2))
%             assert(isequal(y_map, y_map2))
            3;
            
%             [x_pin, y_pin] = findPinwheels(x_map,y_map,Z_map);
            x_pin = [0]; y_pin = [0];


            if firstTime
                 [x_cent, y_cent] = meshgrid(xs_cent,ys_cent);
                 z_cent = column_height_um(ones(numel(x_cent), 1)) /2;
                 if make3Dmap
                     Xcent = [x_cent(:)'; y_cent(:)'; z_cent(:)' ];
                 else
                     Xcent = [x_cent(:)'; y_cent(:)' ];
                 end
                 
                 idx_x_cent = binarySearch(x_map, Xcent(1,:)); idx_y_cent = binarySearch(x_map, Xcent(2,:));
                 idx_cent = sub2indV(size(ori_map), [idx_y_cent(:), idx_x_cent(:)]);
                 
                 
                 if make3Dmap
                     nZ = round(column_height_um / um_spacing);
                     z_map = linspace(0, column_height_um, nZ);
                     [x_samp, y_samp, z_samp] = meshgrid(x_map,y_map,z_map);
                     Xsamp = [x_samp(:)'; y_samp(:)'; z_samp(:)'];
                 else
                     [x_samp, y_samp] = meshgrid(x_map,y_map);
                     Xsamp = [x_samp(:)'; y_samp(:)'];
                 end
                 nSampPoints = size(Xsamp, 2);
                 
                 
                 idx_x_samp = binarySearch(x_map, Xsamp(1,:)); idx_y_samp = binarySearch(x_map, Xsamp(2,:));
                 idx_samp = sub2indV(size(ori_map), [idx_y_samp(:), idx_x_samp(:)]);
                 
                 
                 for ci = 1:nCentPoints
                     % pick a point close to the center
                                              
                     dists_ci = normV( bsxfun(@minus, Xsamp, Xcent(:,ci)), 1);
                     idx_inRange = find( dists_ci < rangeIntegrate_um );
                     
                     idx_inRange_fromCent{ci} = idx_inRange;
                     dists_fromCent{ci} = dists_ci(idx_inRange);
                 end

                 firstTime = 0;
            end
                    

            %% plot map
            if showMaps
                %%
                fig_id = 1;
                figure(fig_id); clf; % phase angle of Z (orientation map).
                imagesc(x_map*plot_scale_factor, y_map*plot_scale_factor, rad2deg(ori_map)); 

                colormap(hsv(250));
                h_col = colorbar; 
                set(h_col, 'ytick', [0 : 45 : 180], 'ylim', [0, 179.9]);
                
                axis xy square;
                hold on;
                dens = length(x_pin) / n_wav^2;
    %         fprintf('pinwheel density = %.2f\n', dens);
            end
            3;
        
            
            x_cent_ori = ori_map(idx_cent);
            if showMaps
%                 plot(Xcent(1,:), Xcent(2,:), 'w.');
%                 plot(x_pin, y_pin, 'ko');
                plot(Xcent(1,:)*plot_scale_factor, Xcent(2,:)*plot_scale_factor, 'k.');

            end
            %%
            Xpin = [x_pin; y_pin];

            dists_to_pinwheels = zeros(1, nCentPoints);
            for zi = 1:nCentPoints
                dists_to_all_pinwheels = normV(bsxfun(@minus, Xcent(1:2,zi), Xpin),1);
                dists_to_pinwheels(zi) = min(dists_to_all_pinwheels);        
            end

            %%
        %     rng_full = [0, lambda*n_wav];
        %     Xsamp = rng_full(1) + rand(2, nSampPointsPerCent)*(rng_full(2)-rng_full(1));
              
            x_samp_ori = ori_map(idx_samp);

            if ~make3Dmap  % check indexing is correct
                assert(isequal(x_samp_ori(:), ori_map(:)));            
            end

            if showMaps
        %     plot(Xsamp(1,:), Xsamp(2,:), 'w.', 'markersize', 1);
        
                    set(fig_id, 'windowstyle', 'normal');
                    set(fig_id, 'color', 'w', 'position', [figureLB + figureLB_offset*1, figureSize]);
                    
                    set(gca, 'xtick', [0:5], 'ytick', [0:5]);
                    %%
                    if saveSampleOriMapFig
                        hgsave(fig_id, sprintf('%sFig13A_simOri_sampleOriMap', figureFolder))
                    end
                    3;

            end


            showSingleElectrodeHistograms = 0;
            showFullMapHistograms = strcmp(runmode, 'diff_hist');
%             showHistograms = (map_i == 10);
            plotProportionOfOutliers = strcmp(runmode, 'diff_hist');
            %%
            
            %%



            switch seeing_distance_func
                case 'exponential', 
                    if useR50inCalc
                        r50_factor = exp_r50_factors(nDims);
                    else
                        r50_factor = 1;
                    end
                    seeing_dist_func = @(r,a) (1./a)*exp(- (r * r50_factor ) ./a);

                case 'gaussian',
                    if useR50inCalc
                        r50_factor = gauss_r50_factors(nDims);
                    else
                        r50_factor = 1;
                    end
                    seeing_dist_func = @(r,a) gaussian(r, 0, a/r50_factor);
            end

    
            ori_bins_E = linspace(0, 90, 50);
            binC = binEdge2cent(ori_bins_E);

            gauss_std_dev = 5 / 0.6745;

    %         dists = zeros(nCentPoints, nSampPoints);        
    %         diff_ori = zeros(nCentPoints, nSampPoints);
    %         idx_inRange = cell(1, nCentPoints);

            
            for ci = 1:nCentPoints                

                idx_inRange = idx_inRange_fromCent{ci};                        
                diffOri_fromCent{ci} = rad2deg( circDist( x_cent_ori(ci), x_samp_ori(idx_inRange)', pi) );
                    3;
            end
            


            for dist_i = 1:nDists;
%%
                all_dOris_C = cell(1, nCentPoints);
                all_wgts_C = cell(1, nCentPoints);
                all_wgt_dist_C = cell(1, nCentPoints);
                all_dists_C = cell(1, nCentPoints);
                for ci = 1:nCentPoints
                    progressBar;
                    % pick a point close to the center            

                    %%
                    dists_ci = dists_fromCent{ci};
                    diff_ori_ci = diffOri_fromCent{ci};


                    wgts = seeing_dist_func(dists_ci, dist_vals(dist_i));
                    wgts = wgts/sum(wgts);

                    wgts_min_keep = 1e-7;

        %                 figure(55); clf;
        %                 hist(log(wgts), 40);
                    idx_keep = wgts > wgts_min_keep;

                    %% check distance
                                
                    if showSingleElectrodeHistograms
                        %%

                        v = dists_ci(idx_keep);
    %                     v = diff_ori_ci(idx_keep);
                        ori_bins_E = linspace(0, max(v), 50);
                        binC = binEdge2cent(ori_bins_E);

                        w = wgts(idx_keep);
    %                     [wgt_md, wgt_md_itp] = weightedMedian(diff_ori, wgts);
                        [wgt_md, wgt_md_itp] = weightedMedian(v(:), w(:));


                        n = histcnt(v, ori_bins_E, w);
                        figure(55); clf;
                        h55 = bar(binC, n, 1);
                        xlim(ori_bins_E([1, end]));                                
                        drawVerticalLine(wgt_md_itp, 'color', 'k');
                        
                    end

        %             [wgt_md, wgt_md_itp] = weightedMedian(diff_ori, wgts);
        %             [wgt_md2, wgt_md_itp2] = weightedMedian(diff_ori(idx_keep), wgts(idx_keep));                

                    all_dOris_C{ci} = diff_ori_ci(idx_keep);
                    all_wgts_C{ci} = wgts(idx_keep);
                    all_dists_C{ci} = dists_ci(idx_keep);

                    wgt_md_dist = weightedMedian(dists_ci(idx_keep), wgts(idx_keep));
                    all_wgt_dist_C{ci} = wgt_md_dist;

                end
                %%
                all_dOris = [all_dOris_C{:}];
                all_wgts = [all_wgts_C{:}];
                all_dists = [all_dists_C{:}];

                % 1. naive weighted median.
                wgt_med_all = weightedMedian(all_dOris, all_wgts);
                
                % 2. weighted median, convolving orientation with gaussian.                
                nOriBins = 90;
                ori_bins_E = linspace(0, 90, nOriBins+1);
                ori_bins_C = binEdge2cent(ori_bins_E);

                nDistBins = 50;
                dist_bins_E = linspace(min(all_dists), max(all_dists), nDistBins+1); 
                dist_bins_C = binEdge2cent(dist_bins_E);                                                                                        
                
                [wgt_med_conv, wgt_med_unconv, pct_outlr, n_conv] = getMedian_convolvedWithGaussian(all_dOris, all_wgts, gauss_std_dev, ori_bins_C);
                n_conv_origBin = n_conv;
%                 dev = abs(wgt_med_all - wgt_med_unconv) / wgt_med_all;
    %                 assert( dev < .02);

    %                 for ci = 1:nCentPoints
    %                     median_values_eachCentCell(dist_i, ci, map_i) = ...
    %                         getMedian_convolvedWithGaussian(all_dOris_C{ci}, all_wgts_C{ci}, gauss_std_dev);                    
    %                 end
                
                median_values(dist_i, map_i) = wgt_med_conv;
                pct_outliers(dist_i, map_i) = pct_outlr;                
%                 median_values_unconv(dist_i, map_i) = wgt_med_unconv;
    
                % 3. weighted median, correcting for 3D distance, & convolving orientation with gaussian 
                calculate2Dhist = plotProportionOfOutliers;
                if calculate2Dhist
                    applySeeingDistanceWeightsWhenHistogramming = 0;
                    applySeeingDistanceWeightsAtEnd = 1;

                    do3DdistanceReweighting = 1;
                    doOriGaussianConv = 1;

    %                 ori_bins_E = reBinE; ori_bins_C = binEdge2cent(ori_bins_E);
    %                 nOriBins = length(ori_bins_C);                                                

                    if applySeeingDistanceWeightsWhenHistogramming
                        hist_wgts = all_wgts;
                    else
                        hist_wgts = [];
                    end

                    binV_xy = hist2d(all_dOris, all_dists, ori_bins_E, dist_bins_E, hist_wgts);                
%                     binV_xy1 = binV_xy;

                    if do3DdistanceReweighting
                        binV_xy = reweightWith3Dweighting(  binV_xy, dist_bins_C );
                    end
    %                 binV_xy2 = binV_xy;

                    if doOriGaussianConv
                        binV_xy = convolveBinsWithGaussian(ori_bins_C, binV_xy, gauss_std_dev);
                    end
    %                 binV_xy3 = binV_xy;

                    seeing_distance_wgts = seeing_dist_func(dist_bins_C, dist_vals(dist_i));

                    if applySeeingDistanceWeightsAtEnd
                        binV_xy = bsxfun(@times, binV_xy, seeing_distance_wgts');
                    end
                    
                end
                
                
%                 binV_diffOri1 = normalized(sum(binV_xy1, 2));
%                 binV_diffOri2 = normalized(sum(binV_xy2, 2));
%                 binV_diffOri3 = normalized(sum(binV_xy3, 2));
%                 binV_diffOri = normalized(sum(binV_xy, 2));
%                                 
%                 binV_diffOri = sum(binV_xy, 2);
%                 binV_dist = sum(binV_xy, 1);
%                 
%                 n_norm_conv = n_norm_conv/sum(n_norm_conv);
%                 binV_diffOri = binV_diffOri/sum(binV_diffOri);
%                 figure(44); clf;
%                 
%                 plot(ori_bins_C, n_norm_conv, 'b.-', ...
%                      ori_bins_C, binV_diffOri, 'ro-' ...
%                  );

%                 plot(ori_bins_C, binV_diffOri1, 'go-', ...
%                      ori_bins_C, binV_diffOri2, 'rs-', ...
%                      ori_bins_C, binV_diffOri3, 'm*-', ...
%                      ori_bins_C, binV_diffOri,  'k.-' ...                     
%                  );
                
%                 plot(ori_bins_C, binV_diffOri1, 'b.-', ...
%                      ori_bins_C, binV_diffOri2, 'ro-' ...                     
%                  );
                
%                 med_diffOri = getMedianFromBins(ori_bins_C, binV_diffOri);
                
%                 idx_ori_outlier = ori_bins_C > 45;
%                 pct_outlr = sum(binV_diffOri(idx_ori_outlier)) / sum(binV_diffOri);
                
%                 median_values(dist_i, map_i) = med_diffOri;
%                 pct_outliers(dist_i, map_i) = pct_outlr;

                  %%              
                
                
                %%
%                 n_norm_conv_orig = n_norm_conv;
                %%
                if showFullMapHistograms
                    %%
                    
                    alsoShowUnconvolvedDistribution = 0;
                    showUnconvolvedInstead = 0;
                    showTxtPctOutliers = 1;
                    
                    bar_norm_col = [.5, .5, .5]; bar_out_col = [.8, .8, .8];
                    
                    n_unconv = histcnt(all_dOris, ori_bins_E, all_wgts);                     
                    n_unconv = n_unconv / sum(n_unconv) / diff(ori_bins_C(1:2));
                    
                    n_conv = n_conv / sum(n_conv) / diff(ori_bins_C(1:2));
                                        
                    ori_bins_plot_E = 0 : 2 : 90;
                    ori_bins_plot_C = binEdge2cent(ori_bins_plot_E);
                    
                    closest_bin_idxs = arrayfun(@(origBinC) binarySearch(ori_bins_plot_C, origBinC), ori_bins_C);
                    [~, binIdxList] = uniqueList( closest_bin_idxs );                    
                    
%                     ori_bins_E = binCent2edge(binC_conv);
%                     binE = linspace(0, 90, 5);
                    
                    
%                     med = getMedianFromBins(ori_bins_C, n);
                    
%                     wgt_med2 =  interp1(c, binE, 0.5);

                    n_conv = n_conv_origBin;
                    n_unconv_origBin = n_unconv;
                                        
                    n_conv   = cellfun(@(lst) sum(n_unconv_origBin(lst)), binIdxList);
                    n_unconv = cellfun(@(lst) sum(n_conv_origBin(lst)), binIdxList);                    
                    
                    idx_typical = find(ori_bins_plot_C < 45); idx_outlr = find(ori_bins_plot_C  >= 45);
                    n_conv_split   = [n_conv(idx_typical),   zeros(1, length(idx_outlr)); zeros(1, length(idx_typical)), n_conv(idx_outlr);]';
                    n_unconv_split = [n_unconv(idx_typical), zeros(1, length(idx_outlr)); zeros(1, length(idx_typical)), n_unconv(idx_outlr);]';
                    
                    if showUnconvolvedInstead
                        binVals = n_unconv_split;
                        wgt_md_use = wgt_med_unconv;
                    else
                        binVals = n_conv_split;
                        wgt_md_use = wgt_med_conv;
                    end
                    %%
                    fig_id = 56;
                    figure(fig_id); clf; 
                    set(fig_id, 'windowstyle', 'normal');
                    set(fig_id, 'color', 'w', 'position', [figureLB + figureLB_offset*1, figureSize]);
                    
                    if binVals(5,1) < binVals(6,1)
                        [binVals(5,1), binVals(6,1)] = deal(binVals(6,1), binVals(5,1));
                    end
                    
                    if alsoShowUnconvolvedDistribution
                        bar(ori_bins_plot_C, n_unconv_split, 1);                        

%                     bar(binC_conv, n_conv_norm, 1);
%                     drawVerticalLine(wgt_med_all, 'color', 'g');
                    
                        plot(ori_bins_plot_C, n_conv_split, 'r--', 'linewidth', 2);                    
%                     plot(binC, n_norm, 'r.-');                    
                    
                    else
                        
                        h_bar = bar(ori_bins_plot_C, binVals, 1, 'stacked');                        
                        set(h_bar(1), 'facecolor', bar_norm_col);
                        set(h_bar(2), 'facecolor', bar_out_col);                        
                        
                    end

                    
                    
                    xlim(ori_bins_E([1, end]));     set(gca, 'xtick', [0:15:90]);
%                     title(sprintf('Distribution of differences ( \\sigma = %.1f µm)', dist_vals(dist_i)));
%                     xlim([0, 45]); set(gca, 'xtick', [0:5:90]);
                    
%                     drawVerticalLine(wgt_med2, 'color', 'g', 'linestyle', ':');
                    xlabel('Difference in orientation'); ylabel('(weighted) Count');
                    drawVerticalLine(wgt_md_use, 'color', 'k', 'linewidth', 2);
                    
                    if alsoShowUnconvolvedDistribution
                        legend('distribution', 'conv(dist)', 'median');
                    else
                        legend('Differences < 45°', 'Differences > 45°', 'median');
                    end
%                     title(sprintf('Seeing distance = %.1f \\mum. Median Diff = %.2f\\circ', dist_vals(dist_i), wgt_med_all));
                    
%                     title(sprintf('Distribution of differences (\\sigma = %.1f \\mum) Conv Wgt Median Diff = %.2f\\circ', dist_vals(dist_i), wgt_med_conv));
%                     
%                     drawnow;
                    
                    ylims = ylim;
                    text(wgt_md_use, ylims(2)*[0.5], sprintf('\\leftarrow median = %.1f°', wgt_md_use+.1));
                    
                    if showTxtPctOutliers
                        %%
                        s = sum(binVals);
                        pct_out = s(2)/sum(s);
                        text(mean([45, 90]), ylims(2)*[.1], sprintf('Outliers : %.1f %%', pct_out*100), 'horiz', 'center')                        
                    end
                    set(fig_id, 'color', 'w', 'position', [figureLB + figureLB_offset*2, figureSize]);

                    if saveDiffHistFig
                        hgsave(fig_id, sprintf('%sFig13B_simOri_diffOriHist', figureFolder))
                    end
                    
                    
                    show2Dhist = 0;
%%
                    if show2Dhist
                        %%
%                         xy = [all_dOris(:), all_dists(:)];
                        
                        figure(77); clf;
                        imagesc(ori_bins_C, dist_bins_C, binV_xy'); axis xy; colorbar;
                        xlabel('Difference in orientation'); set(gca, 'xtick', 0:15:90); ylabel('Distance (\mum)');                        
                        
%                         figure(78); clf;
%                         imagesc(ori_bins_C, dist_bins_C, binV_xy_rewgt'); axis xy; colorbar;
%                         xlabel('Difference in orientation'); set(gca, 'xtick', 0:15:90); ylabel('Distance (\mum)');                        
                        
                        %}  
                        ori_diff_dist = sum(binV_xy, 2);
                        ori_diff_dist = ori_diff_dist/sum(ori_diff_dist);% *diff(ori_bins_C(1:2));
%                         ori_bins_C ( find(cumsum(ori_diff_dist) > 0.5, 1))
                        figure(78); clf;                        
                        plot( ori_bins_C, ori_diff_dist, '.-' )
                        
                        distance_dist = sum(binV_xy, 1);
                        distance_dist = distance_dist / sum(distance_dist);% * diff(dist_bins_C(1:2));
%                         dist_bins_C ( find(cumsum(distance_dist) > 0.5, 1))
                        figure(79); clf;
                        plot( dist_bins_C, distance_dist, '.-' )                                                
                        
                        %%
                        % proportion of outliers vs distance.
                        idx_outlier = ori_bins_C > 45;
                        nOutliers_v_dist = sum(binV_xy(idx_outlier, :), 1);
                        nTot_v_dist = sum(binV_xy, 1);
                        pct_outliers_v_dist = nOutliers_v_dist ./ nTot_v_dist;

                        figure(80); clf;
                        plot( dist_bins_C, pct_outliers_v_dist, '.-' );
                        xlabel('Distance from electrode (\mum)'); ylabel('Proportion of outliers');
                        
                        
                    end
                    
                    
                    keyboard;
                end
                3;

                
                if plotProportionOfOutliers
                    %%
                    doNaevePlot = 0;
                    doPrecisePlot = 1;
                    if doNaevePlot
                        %%
                        dOri_typ = all_dOris <= 45;
                        dOri_out = all_dOris > 45;
                        bL = lims(all_dists, .01); bE = linspace(bL(1), bL(2), 30); bC = binEdge2cent(bE);
    %                     bV_dist_typ = histcnt(all_dists(dOri_typ), bE);  
    %                     bV_dist_out = histcnt(all_dists(dOri_out), bE);  
                        bV_dist_typ = histcnt(all_dists(dOri_typ), bE, all_wgts(dOri_typ));  
                        bV_dist_out = histcnt(all_dists(dOri_out), bE, all_wgts(dOri_out));  
                        bV_dist_typ_nrm = bV_dist_typ / sum(bV_dist_typ);
                        bV_dist_out_nrm = bV_dist_out / sum(bV_dist_out);
                        wgt_med_dist_typ = weightedMedian(all_dists(dOri_typ), all_wgts(dOri_typ));
                        wgt_med_dist_out = weightedMedian(all_dists(dOri_out), all_wgts(dOri_out));
                        figure(42); clf; hold on; box on;
                        stairs(bC, bV_dist_typ_nrm, 'b-')
                        stairs(bC, bV_dist_out_nrm, 'r-');
                        ylabel('Proportion of cells');
                        xlabel('distance between cell & electrode (\mum)');
                        legend('Typical cells', 'outliers', 'location', 'NE');
                        drawVerticalLine(wgt_med_dist_typ, 'color', 'b', 'linestyle', ':', 'linewidth', 2);
                        drawVerticalLine(wgt_med_dist_out, 'color', 'r', 'linestyle', ':', 'linewidth', 2);
                        title(sprintf('Median Distance: Typical: %.2f \\mum.  Outliers: %.2f \\mum', wgt_med_dist_typ, wgt_med_dist_out));

                        figure(43); clf; hold on; box on;
                        prop_out = bV_dist_out ./ (bV_dist_out + bV_dist_typ);
                        prop_typ = bV_dist_typ ./ (bV_dist_out + bV_dist_typ);
%                         plot(bC, prop_typ, 'bo-');
                        plot(bC, prop_out, 'b.-');
%                         legend('Fraction typical', 'Fraction outliers', 'location', 'E');
                        xlabel('distance between cell & electrode (\mum)');
                        title('Proportion of typical/outlier cells vs distance');
                        drawHorizontalLine([.1, .2, .3], 'linestyle', ':', 'color', 'k')
                    end
                    
                    if doPrecisePlot
                        %%
                        for j = 1:2
                            useMedianPlot = j == 1;

                            idx_outlier = ori_bins_C > 45;
                            nOutliers_v_dist = sum(binV_xy(idx_outlier, :), 1);
                            nTot_v_dist = sum(binV_xy, 1);
                            pct_outliers_v_dist = nOutliers_v_dist ./ nTot_v_dist * 100;                        
                            for i = 1:nDistBins
                                median_dOri_v_dist(i) = getMedianFromBins(ori_bins_C, binV_xy(:,i));
                            end

                            prop_cells_v_dist = sum(binV_xy, 1);
                            prop_cells_v_dist = prop_cells_v_dist/max(prop_cells_v_dist);

                            if useMedianPlot
                                y1 = median_dOri_v_dist;
                                dy_tick = 5;
                                ylab = 'Median ori difference';
                                fig_id = 81;
                                fig_name = 'Fig13G_simOri_medianDiff_vs_dist';
                            else
                                y1 = pct_outliers_v_dist;
                                dy_tick = 5;
                                ylab = 'Percentage of outliers';
                                fig_id = 82;
                                fig_name = 'Fig13H_simOri_propOutliers_vs_dist';
                            end
                            nskp = 2;    
                            idx = 1:nskp:nDistBins;

                            line_w = 1.5;
                            col1 = [1 1 1]*0;
                            col2 = [1 1 1]*.5;

                            figure(fig_id); clf;
                            set(fig_id, 'color', 'w', 'windowstyle', 'normal')
                            set(fig_id, 'position', [figureLB + figureLB_offset*(6+j), figureSize]);
                            [h_ax, h1, h2] = plotyy(dist_bins_C(idx), y1(idx), dist_bins_C(idx), prop_cells_v_dist(idx));
    %                         plot( dist_bins_C, prop_outliers_v_dist, '.-' );

                            set(h_ax, 'outerposition', getNormPosition(1,1,1,1, [.03, 0, .08], [.09, 0, .09]));
                            ylim1 = [0 max(y1)*1.1];
                            set(h_ax, 'box', 'off', 'nextplot', 'add');
                            set(h_ax(2), 'XTickLabel','','XAxisLocation','Top');
                            

                            set(h_ax(1), 'ycolor', col1, 'ylim', ylim1, 'ytick', [0 : dy_tick : ylim1(2)]);
                            set(h_ax(2), 'ycolor', col2);
                            set(h1, 'color', col1, 'marker', 's', 'linewidth', line_w);
                            set(h2, 'linestyle', 'none', 'marker', 'o', 'linewidth', line_w, 'color', col2);                        

                            r = dist_bins_C; sig = dist_vals(dist_i);
                            prop_cells_exact = r.^2 .* exp(-r.^2/(2*sig.^2));
                            prop_cells_exact = prop_cells_exact /max(prop_cells_exact);

                            h3 = plot(r, prop_cells_exact, ':', 'linewidth', line_w, 'parent', h_ax(2), 'color', col2);

                            ylabel(ylab, 'parent', h_ax(1));
                            ylabel('Relative proportion of cells', 'parent', h_ax(2));                       
                            xlabel('Distance from electrode (\mum)'); 
                        
                            set(h_ax, 'Xlim', [0 400]);
                            if saveMedianVsDistFigs
                                hgsave(fig_id, sprintf('%s%s', figureFolder, fig_name))
                            end
%                             xlim([0 400])
                            

                        end
                            3;
                        keyboard;
                    end
                    
                end
                
                3; 


            end   
        %%

        end
        %%
        
        idx_map_closest_to_target = indmin(abs(median_values - median_dOriMU_observed));
        
        doPlots = strcmp(runmode, 'plot_v_sigma');
        if doPlots
            plots_to_do = 1;
        else
            plots_to_do = [];
        end
        
        
        useUnconvMedians = 0;
        
        doOutliersPlot = 1;
        outputText = 0;
    %     plots_to_do = 1:3;
        
        for ii = plots_to_do
            %%
            median_diff_fig = 50+ii;
            pct_outliers_fig = 150+ii;
            if doPlots
                figure(median_diff_fig); clf;
                set(median_diff_fig, 'windowstyle', 'normal', 'color', 'w', 'position', [figureLB + figureLB_offset*(3), figureSize])
                figure(pct_outliers_fig); clf;
                set(pct_outliers_fig, 'windowstyle', 'normal', 'color', 'w', 'position', [figureLB + figureLB_offset*(4), figureSize])
            end
%%
            if useUnconvMedians
                median_values_use = median_values_unconv;
            else
                median_values_use = median_values;
            end
            
            fig_ids = [median_diff_fig, pct_outliers_fig];
            y_vals = {median_values_use, pct_outliers_fig};
            y_labels = {'Median ori difference', 'Percentage of outliers'};
            
%             for fig_i = 1:length(fig_ids)
            
            switch ii
                case 1, median_values_use = median_values;
                case 2, median_values_use = median_values_unconv;
%                 case 3, median_values_use = median_values_unconv_orig;
            end

%             r50_vals_microns = dist_vals;

            dist_vals_X = dist_vals(:)';

            median_values_M = mean(median_values_use, 2);            
            pct_outliers_M = mean(pct_outliers, 2);
            
            median_values_S = std(median_values_use, [], 2);            
            pct_outliers_S = std(pct_outliers, [], 2);
            
            if doPlots
                if nMaps > 1                
                
                    figure(median_diff_fig);
                    h1 = errorbar(dist_vals_X, median_values_M, median_values_S, 'o-'); hold on;                    
                    
                    if doOutliersPlot
                        figure(pct_outliers_fig);
                        h2 = errorbar(dist_vals_X, pct_outliers_M, pct_outliers_S, 'o-'); hold on;
                    end
                    
                else                                
                    figure(median_diff_fig);
                    h1 = plot(dist_vals_X', median_values_M, 's-'); hold on;
                    
                    if doOutliersPlot
                        figure(pct_outliers_fig);
                        h2 = errorbar(dist_vals_X, pct_outliers_M, pct_outliers_S, 'o-'); hold on;
                    end
                end
                set([h1, h2], 'color', 'k', 'marker', '.')           
            end
            
                        %%
            med_targets = median_dOriMU_observed;  nTargets = length(med_targets);
            
            a_targets = zeros(1, nTargets);
            pct_outlier_targets = zeros(1, nTargets);
            
            if outputText
                fprintf('\n\n');
            end
            
            for ti = 1:nTargets
                if outputText
                    fprintf('For median difference of %.1f degrees: vs seeing distance (D):\n', med_targets(ti));
                end
                                
                a_targets(ti) = interp1(median_values_M, dist_vals, med_targets(ti));
                pct_outlier_targets(ti) = interp1(dist_vals, pct_outliers_M, a_targets(ti));
                    
                if isnan(a_targets(ti))
                    beep;
                    fprintf('nan at %d    ', rep_i);
                    keyboard;
                end
                
                if doPlots
                    figure(median_diff_fig);
                    h3 = plot(a_targets(ti), med_targets(ti), 'ko', 'markerfacecolor', [.7 .7 .7]);                
                    ylims = ylim; xlims = xlim;
                    line([xlims(1), a_targets(ti)], med_targets(ti)*[1, 1], 'linestyle', '-', 'color', 'k');
                    line(a_targets(ti)*[1, 1], [ylims(1), med_targets(ti)], 'linestyle', ':', 'color', 'k');
                
                    figure(pct_outliers_fig);
                    h4 = plot(a_targets(ti), pct_outlier_targets(ti), 'ko', 'markerfacecolor', [.7 .7 .7]);
                    ylims = ylim; xlims = xlim;
                    line(a_targets(ti)*[1, 1], [ylims(1), pct_outlier_targets(ti)], 'linestyle', '-', 'color', 'k');
                    line([xlims(1), a_targets(ti)], pct_outlier_targets(ti)*[1, 1], 'linestyle', ':', 'color', 'k');
                end
                if outputText
                    fprintf(' %s : D= %.1f.  ', Qstr, a_targets(ti));
                end
                
                if outputText
                    fprintf('\n');
                end

            end
            
            if doPlots
                
                if useR50inCalc
                    dist_xlabel = 'R_{50} seeing distance (\mum)';
                else
                    if strcmp(seeing_distance_func, 'gaussian')
%                         dist_xlabel = '\sigma of Gaussian seeing distance function (in \mum)';
                        dist_xlabel = '\sigma (\mum)';
                    elseif strcmp(seeing_distance_func, 'exponential')
                        dist_xlabel = '\mu of Exponential seeing distance function (in \mum)';
                    end            
                end
                %%
                figure(median_diff_fig);
                xlabel(dist_xlabel);
                ylabel('Median ori difference');                
                xlim(lims(dist_vals, .05));

                if saveMedOriPctOutlierVsSigmaFigs 
                    hgsave(median_diff_fig, sprintf('%sFig13C_simOri_median_vs_sigma', figureFolder))
                end
                
%                 drawHorizontalLine(med_targets, 'color', 'k', 'linestyle', ':');
                
                figure(pct_outliers_fig);
                xlabel(dist_xlabel);
                ylabel('Percentage of outliers');
                xlim(lims(dist_vals, .05));           
                if saveMedOriPctOutlierVsSigmaFigs
                    hgsave(pct_outliers_fig, sprintf('%sFig13D_simOri_pct_outliers_vs_sigma', figureFolder))
                end


                %%
%                 figure(60+ii); clf;
%                 idx_best_s = indmin( abs(median_values_M(end,:) - med_targets) );
%                 samp_meds = squeeze( median_values_use(end, idx_best_s, :) );
%                 hist(samp_meds);

%                 ci = get95ConfInt(samp_meds, 0.05);
%                 drawVerticalLine(ci, 'color', 'r')
%                 keyboard;
            end
            3;
            if strcmp(runmode, 'plot_v_sigma')
                %%
                interpMedianAtSigmaEst = interp1(dist_vals_X, median_values_use, curEstimateOfGaussSigma);
                interpPoutAtSigmaEst = interp1(dist_vals_X,   pct_outliers, curEstimateOfGaussSigma);
                relDiffSigma = abs(interpMedianAtSigmaEst-median_dOriMU_observed)/median_dOriMU_observed;
%                 relDiffPout  = abs(interpPoutAtSigmaEst-pctOutliers_observed)/pctOutliers_observed;
%                 relDiffPout  = abs(interpPoutAtSigmaEst-pctOutliers_observed)/pctOutliers_observed;
                
%                 [~, idx_bestMaps] = sort( relDiffSigma .* relDiffPout );
                [~, idx_bestMaps] = sort( relDiffSigma );
                idx_bestMap = idx_bestMaps(1);
                idx_secondBestMap = idx_bestMaps(2);
                fprintf('Best map: idx = %d. Median = %.2f (close to %.2f). Pout = %.2f (close to %.2f)\n', idx_bestMap, ...
                    interpMedianAtSigmaEst(idx_bestMap), median_dOriMU_observed,   interpPoutAtSigmaEst(idx_bestMap), pctOutliers_observed );
                fprintf('Second best map: idx = %d. Median = %.2f (close to %.2f). Pout = %.2f (close to %.2f)\n', idx_secondBestMap, ...
                    interpMedianAtSigmaEst(idx_secondBestMap), median_dOriMU_observed,   interpPoutAtSigmaEst(idx_secondBestMap), pctOutliers_observed );
                3;
                %%
                keyboard
            end
        
        end
        
        if abs(a_targets-70.62) < .1
            3;
        end
        median_glob_estimates(rep_i) = a_targets;
        pct_outlier_estimates(rep_i) = pct_outlier_targets;
        3;
    end
%     drawVerticalLine(a_targets(:), 'color', 'k', 'linestyle', ':');
    %%
%     figure(42);
%     plot(a_targets(1:end-1, :), 'o-');
    3;
    fprintf('\n\n');
    nBins = 12;
    alpha = 0.05;
%     med_est_use = median_glob_estimates(1:100);
    med_est_use = median_glob_estimates;
    r50_estimates = med_est_use * gauss_sig_to_R50_3D_factor;
    
    fig_id1 = [70];
    xdata = {med_est_use, r50_estimates, pct_outlier_estimates};
    xlabels = {'\sigma (\mum)',  'r50', 'Percentage of outliers'};
    
    saveFileNames = {'Fig13E_simOri_distSigma', '', 'Fig13F_simOri_distPctOutliers'};
    labl_short = {'\sigma', 'r_{50}', 'p_{out}'};
    suffix = {'\mum', '\mum', '%'};
%%
    doTitle = 0;
    barCol = [.7 .7 .7]; 
    for fig_i = 1:length(xdata)
        fig_id = fig_id1 + fig_i;
        figure(fig_id);
        X = xdata{fig_i};
    
        [n,x] = hist(X, nBins);
        hh(fig_i) = bar(x,n, 1);
        set(hh(fig_i), 'facecolor', barCol);
        mn = mean(X);
        ci   = get95ConfInt(X, alpha);
        xlabel(xlabels{fig_i});
        s = sprintf('Mean %s: %.2f %s [95%% ci: %.2f - %.2f %s]', labl_short{fig_i}, mn, suffix{fig_i}, ci, suffix{fig_i});
        fprintf('%s\n', s)
        if doTitle
            title(s);
        end
        ylabel('Count');
        
        drawVerticalLine(mn, 'color', 'k', 'linewidth', 2);
        drawVerticalLine(ci, 'color', 'k', 'linestyle', ':', 'linewidth', 2);       
        
        set(fig_id, 'color', 'w', 'windowstyle', 'normal', 'position', [figureLB + figureLB_offset*(5+fig_i), figureSize]);
        if ~isempty(saveFileNames{fig_i}) && saveDistSigmaOutlierFigs
            hgsave(fig_id, sprintf('%s%s.fig', figureFolder, saveFileNames{fig_i}))
            
        end
    end
    
    3;
%     title(s);
%     xlabel('');
    
%     title('Distribution of 
%     fprintf('%s\n', s);
    %%
    sigma_mean = mean(median_glob_estimates);
    idx_closest = indmin(abs(median_glob_estimates - sigma_mean));       
    fprintf('Index of repetition which had the median that was closest (to the mean of %.2f) was: %d\n', sigma_mean, idx_closest);
%     figure(71); clf;
        
    %%
    
%     fprintf('Mean pct outlier = %.1f (ci = %.1f - %.1f)\n', pct_outlier_mean, pct_outlier_ci);

    if strcmp(runmode, 'confint')
        keyboard;
    end
    3;

end


function ci = get95ConfInt(x, alpha) 
    %%
    uX = unique(x);
    n = histc(x, uX);
    
    cdf = cumsum(n) / sum(n);
    i1 = find(cdf > alpha/2, 1, 'first');
    i2 = find(cdf < 1-alpha/2, 1, 'last');
    
    ci = uX([i1, i2]);


end

% drawLinesFromAxesToPoint(ax, pt);

function [med_conv, med_unconv, pct_outliers, n_conv] = getMedian_convolvedWithGaussian(all_dOris, all_wgts, gauss_std, binC)
        %%
    nBins = length(binC);
    binE = binCent2edge(binC);
    
    n = histcnt(all_dOris, binE, all_wgts);    
    
    med_unconv = getMedianFromBins(binC, n);
        
    binE_ext = unique([-binE, binE]);    
    binC_ext = binEdge2cent(binE_ext);
    n_ext = [fliplr(n), n]/2;
    
    idx_orig = nBins+1 : nBins*2;
%     idx_rev  = nBins :-1 :1;    
    
    GaussFunc = gaussian(binC_ext, 0, gauss_std);
    GaussFunc = GaussFunc/sum(GaussFunc);
    n_ext_conv = conv(n_ext, GaussFunc, 'same');
    
    n_conv = n_ext_conv(idx_orig);
    med_conv = getMedianFromBins(binC, n_conv);
    
    
    
    idx_outlier = binC > 45;    
    pct_outliers = sum( n_conv(idx_outlier)  ) / sum( n_conv) * 100;    

    show = 0;
    if show        
        %%
        figure(59); clf; hold on;
        plot(binC, n, 'b.-');
        plot(binC, n_conv, 'r.-');
        xlim(binE([1, end]));
        drawVerticalLine(med_unconv, 'color', 'b', 'linestyle', ':');
        drawVerticalLine(med_conv, 'color', 'r', 'linestyle', ':');
        xlabel('Difference in orientation'); ylabel('(weighted) Count');
        legend('original', 'convolved');
        title(sprintf('Median:: unconvolved: %.1f, convolved: %.1f', med_unconv, med_conv));
        3;
    end
        
end
    

function med = getMedianFromBins(binC, dist)
    if length(binC) ~= length(dist)
        error('bins and distribution must be the same length')
    end
    binE = binCent2edge(binC);
    
    dist_norm = dist(:)/sum(dist);
    cumDist = [0; cumsum(dist_norm)];
    
    med = getItpMedian(binE, cumDist);
        
end

function med = getItpMedian(binEdges, cumDist)
    %%
    [uCumDist, idx_first] = unique(cumDist, 'first');

    med = interp1(uCumDist, binEdges(idx_first), 0.5);
    
end


function binV_xy_rewgt = reweightWith3Dweighting( binV_xy, dist_bins_C )
    binV_xy_rewgt = zeros(size(binV_xy));
%     dist_bins_E = binCent2edge(dist_bins_C);
    nDistBins = length(dist_bins_C);
    assert(size(binV_xy, 2) == nDistBins);
    nOriBins= size(binV_xy, 1);
    
    dh = diff(dist_bins_C(1:2));
        
    for h_i = 1:nDistBins
        %%
        h = dist_bins_C(h_i);
        
        nBinsFurther = nDistBins - h_i + 1; % include current bin
        zs = zeros(1,nBinsFurther);
        zi = 1:nBinsFurther;
        zs(zi) = sqrt( (h + (2*zi - 1)*(dh/2)) .^2 - h^2  );
        
        wgts = [zs(1), diff(zs)];
        
        for j = 1:nOriBins
            binV_xy_rewgt(j, h_i + [0:nBinsFurther-1]) = binV_xy_rewgt(j, h_i + [0:nBinsFurther-1]) + wgts * binV_xy(j, h_i);
        end
        
    end
    
end

function binV_conv = convolveBinsWithGaussian(binC, binV, gauss_std_dev, dim)
                                      %%    
    nBins = length(binC);
    if nargin < 4
        dim = find(size(binV)>1,1);
    end

    if ismatrix(binV)
        assert(size(binV, dim) == nBins)        
        if dim == 2
            binV = binV'; % convolve columns
        end
    else        
        binV = binV(:);
    end
    %%
    
%     binE = linspace(0, 90, nBins+1);
%     binC = binEdge2cent(binE);
    
    binC_ext = sort([-binC(:); binC(:)]);
    binV_ext = [flipud(binV); binV];

    idx_orig = nBins+1 : nBins*2;
    assert(isequal(binC(:), binC_ext(idx_orig)));

    GaussFunc = gaussian(binC_ext, 0, gauss_std_dev);
    GaussFunc = GaussFunc/sum(GaussFunc);
    nCols = size(binV_ext, 2);
    if nCols == 1
        binV_ext_conv = conv(binV_ext, GaussFunc, 'same');
    else
        binV_ext_conv = zeros(size(binV_ext));
        for i = 1:nCols
            binV_ext_conv(:,i) = conv(binV_ext(:,i), GaussFunc, 'same');
        end
        
    end
    
    binV_conv = binV_ext_conv(idx_orig, :);
    
    if dim == 2
        binV_conv = binV_conv';
    end
    
end

function x = normalized(x)
    x = x/sum(x);    
end