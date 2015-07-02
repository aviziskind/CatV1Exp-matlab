function simulateV1simpleCells

    Xlim = [0, 8];
    Ylim = [0, 8];
    
    redoParamFile = 1;

    gaborParamDistsFile = [CatV1Path 'GaborParamDists.mat'];
    if ~exist(gaborParamDistsFile, 'file') || redoParamFile;
        createGaborParamDistsFile(gaborParamDistsFile);
    end
    S_dists = load(gaborParamDistsFile);
    dists = S_dists.dists;
        
        
    V1_dens = 500;
    V1_oriStd_deg = 0;
    V1_sig_x = 1;
    V1_sig_y = 2;
    V1_k = (0.5)*(2*pi); % cyc/deg  0.5 * (2*pi) = spat freq such that 1 == width of subregion
    V1_phi = 0;
    LGN_radius = 0;
%     
    X_range = diff(Xlim);
    Y_range = diff(Ylim);
    
    nLGNperV1 = 50;
    
    
    NV1Cells = diff(Xlim)*diff(Ylim)*V1_dens;
%     NV1Cells = 3;
    allV1_mu_x = Xlim(1) + X_range*rand(1, NV1Cells);
    allV1_mu_y = Ylim(1) + Y_range*rand(1, NV1Cells);
    
    allV1_sig_x = randFromDist(dists.sig_x, NV1Cells);
    sig_yx_ratios = randFromDist(dists.sig_yx_ratio, NV1Cells);        
    allV1_sig_y = allV1_sig_x .* sig_yx_ratios;
    
%     dists.Dori_fromMU.mu = 5; %dists.Dori_fromMU.mu/5;
    allV1_theta = deg2rad( randFromDist(dists.Dori_fromMU, NV1Cells) );
    allV1_phi = randFromDist(dists.phi, NV1Cells); 
    allV1_k   = randFromDist(dists.k, NV1Cells); 
    
    
%     log_k = log(V1_k); allV1_k   = exp( randn(1,NV1Cells) );
    
    allLGN_x = zeros(nLGNperV1, NV1Cells);
    allLGN_y = zeros(nLGNperV1, NV1Cells);
    allLGN_sgn = zeros(nLGNperV1, NV1Cells);
    
    
    show = 0;
    for i = 1:NV1Cells
        LGN_x = allV1_mu_x(i) + randn(1, nLGNperV1)*allV1_sig_x(1);
        LGN_y = allV1_mu_y(i) + randn(1, nLGNperV1)*allV1_sig_y(1);
        LGN_xy = [LGN_x(:) LGN_y(:)];
        G = gabor(1, allV1_mu_x(i), allV1_mu_y(i), allV1_sig_x(i), allV1_sig_y(i), allV1_theta(i), allV1_k(i), allV1_phi(i), 0, LGN_xy);
        LGN_sgn = sign(G);        
        if show
            %%
            figure(1); clf; hold on; box on;
            i_ON = LGN_sgn == 1;
            i_OFF = LGN_sgn == -1;
            mksize = 6;

            plot(LGN_x(i_ON), LGN_y(i_ON), 'ro', 'markersize', mksize);
            plot(LGN_x(i_OFF), LGN_y(i_OFF), 'bo', 'markersize', mksize);
            plot(allV1_mu_x(i), allV1_mu_y(i), 'ks', 'markersize', 8);
            3;
        end
        
        allLGN_x(:,i) = LGN_x;
        allLGN_y(:,i) = LGN_y;
        allLGN_sgn(:,i) = LGN_sgn;
        3;
    end
    
    show = 0;
    i_ON = allLGN_sgn == 1;
    i_OFF = allLGN_sgn == -1;
    if show
        %%
        figure(2); clf; hold on; box on;
        mksize = 10;
        plot(allLGN_x(i_ON), allLGN_y(i_ON), 'ro', 'markersize', mksize);
        plot(allLGN_x(i_OFF), allLGN_y(i_OFF), 'bo', 'markersize', mksize);
%         plot(allV1_mu_x, allV1_mu_y, 'k*', 'markersize', 12, 'linewidth', 2);                
        
%         w = 1;
%         plot(10, 10, 'ko', 'markersize', 10);
%         h_rect = rectangle('position', [9.5, 9.5, 1, 1]);
        h_rect = rectangle('position', [Xlim(1), Ylim(1), X_range, Y_range]);
        title(sprintf('LGN ON- and OFF- afferents from %d cells',NV1Cells))
        axis square;
    end

           3;
    
    %%
    nBins = 100;
    xedges = linspace(Xlim(1), Xlim(2),nBins+1);
    yedges = linspace(Ylim(1), Ylim(2),nBins+1);

    x_cents = binEdge2cent(xedges);
    y_cents = binEdge2cent(yedges);
    
%     if LGN_radius == 0
    
%         mHist2d_ON2 = hist2d([allLGN_x(i_ON), allLGN_y(i_ON)],yedges,xedges)';
%         mHist2d_OFF2 = hist2d([allLGN_x(i_OFF), allLGN_y(i_OFF)],yedges,xedges)';
%     else
        LGN_radius = diff(x_cents(1:2))/2;
        LGN_radius = 0.5;
        
        [mHist2d_ON, mHist2d_OFF] = deal(zeros(length(x_cents), length(y_cents)));
        LGN_ON_xy = [allLGN_x(i_ON), allLGN_y(i_ON)];
        LGN_OFF_xy = [allLGN_x(i_OFF), allLGN_y(i_OFF)];

        %%
%         tic;
        xedges_hist = xedges; 
        xedges_hist(1) = -LGN_radius;
        xedges_hist(end) = Xlim(2) + LGN_radius;
        
        [~, LGN_ON_x_bin_idxs] = histcnt(LGN_ON_xy(:,1), xedges_hist);
        [~, LGN_ON_y_bin_idxs] = histcnt(LGN_ON_xy(:,2), xedges_hist);

        [~, LGN_OFF_x_bin_idxs] = histcnt(LGN_OFF_xy(:,1), xedges_hist);
        [~, LGN_OFF_y_bin_idxs] = histcnt(LGN_OFF_xy(:,2), xedges_hist);
        
        y_cand_idx_ON_C = cell(1, nBins);
        y_cand_idx_OFF_C = cell(1, nBins);
        
%         toc;
        %%
        R = ceil(LGN_radius/diff(xedges(1:2)));
        
        progressBar('init-', nBins, 40);
        for i = 1:length(x_cents)
            x_cand_idx_ON = ibetween(LGN_ON_x_bin_idxs, [i-R, i+R]);
            x_cand_idx_OFF = ibetween(LGN_OFF_x_bin_idxs, [i-R, i+R]);
                        
            for j = 1:length(y_cents)
                
                if isempty(y_cand_idx_ON_C{j})
                    y_cand_idx_ON  =  ibetween(LGN_ON_y_bin_idxs, [j-R, j+R]) ;
                    y_cand_idx_OFF =  ibetween(LGN_OFF_y_bin_idxs, [j-R, j+R]) ;
                    y_cand_idx_ON_C{j} = y_cand_idx_ON;
                    y_cand_idx_OFF_C{j} = y_cand_idx_OFF;
                else
                    y_cand_idx_ON = y_cand_idx_ON_C{j};
                    y_cand_idx_OFF = y_cand_idx_OFF_C{j};                    
                end
                    
                bin_xy = [x_cents(i), y_cents(j)];
                
                 
                
                idx_cand_ON = x_cand_idx_ON & y_cand_idx_ON;
                idx_ij_ON = normVsqr( bsxfun(@minus, LGN_ON_xy(idx_cand_ON, :), bin_xy), 2) < (LGN_radius^2);                               
                n_ON = nnz(idx_ij_ON);                
                
                idx_cand_OFF = x_cand_idx_OFF & y_cand_idx_OFF;
                idx_ij_OFF = normVsqr( bsxfun(@minus, LGN_OFF_xy(idx_cand_OFF, :), bin_xy), 2) < (LGN_radius^2);
                n_OFF = nnz(idx_ij_OFF);                

                mHist2d_ON(j,i) = n_ON;
                mHist2d_OFF(j,i) = n_OFF;
                
                doChk = 1;
                if doChk
                    idx_ij_ON =  normVsqr( bsxfun(@minus, LGN_ON_xy, bin_xy), 2) < (LGN_radius^2) ;                
                    idx_ij_OFF = normVsqr( bsxfun(@minus, LGN_OFF_xy, bin_xy), 2) < (LGN_radius^2);
                
                    assert(isequal(n_ON, nnz(idx_ij_ON)));                
                    assert(isequal(n_OFF, nnz(idx_ij_OFF)));    
                end
%                 idx_missed = find( idx_ij_ON & ~(x_cand_idx_ON | y_cand_idx_ON) )
                
%                 [[LGN_ON_x_bin_idxs ( idx_ij_ON)], [LGN_ON_y_bin_idxs( idx_ij_ON)]]
%                 [LGN_ON_x_bin_idxs ( idx_ij_ON), LGN_ON_y_bin_idxs( idx_ij_ON)]
                

                
            end
            progressBar;
        end        
        progressBar('done');
        3;
        
        
%     end
%     mHist2d = hist2d(mYX,yedges,xedges);
%
%   nXBins = length(xedges);
%   nYBins = length(yedges);
%   vXLabel = 0.5*(xedges(1:(nXBins-1))+xedges(2:nXBins));
%   vYLabel = 0.5*(yedges(1:(nYBins-1))+yedges(2:nYBins));
    %%

    mHist2d_ON_OFF = mHist2d_ON-mHist2d_OFF;
%     mHist2d_ON_plus_OFF = mHist2d_ON + mHist2d_OFF;
%     mHist2d_ON_OFF2 = mHist2d_ON2-mHist2d_OFF2;
    
    figure(10); imagesc(x_cents, y_cents, mHist2d_ON); axis xy; axis square; colorbar; title(sprintf('LGN ON density (N = %d V1 cells)', NV1Cells))
    figure(11); imagesc(x_cents, y_cents, mHist2d_OFF); axis xy; axis square; colorbar; title(sprintf('LGN OFF density (N = %d V1 cells)', NV1Cells))
    figure(12); imagesc(x_cents, y_cents, mHist2d_ON_OFF);  axis xy; axis square; colorbar; title(sprintf('LGN ON - OFF density (N = %d V1 cells)', NV1Cells))
    clims = caxis; caxis(max(abs(clims))*[-1, 1])


%     figure(20); imagesc(x_cents, y_cents, mHist2d_ON2); axis xy; axis square; colorbar; title(sprintf('LGN ON density (N = %d V1 cells)', NV1Cells))
%     figure(21); imagesc(x_cents, y_cents, mHist2d_OFF2); axis xy; axis square; colorbar; title(sprintf('LGN OFF density (N = %d V1 cells)', NV1Cells))
%     figure(22); imagesc(x_cents, y_cents, mHist2d_ON_OFF2);  axis xy; axis square; colorbar; title(sprintf('LGN ON - OFF density (N = %d V1 cells)', NV1Cells))
    clims = caxis; caxis(max(abs(clims))*[-1, 1])
    %%
    dBin = diff(x_cents(1:2));
    spatPeriod_pix = (2*pi / V1_k) / dBin;    
        
    
    
    all_oris = [0:5:175];  
    all_ph = 0:20:360;     
    
    [ori_tc, R_ori_ph] = getOriTc(mHist2d_ON_OFF, all_oris, spatPeriod_pix, all_ph);
    
    %%
%     frameDims = [nBins_inLims, nBins_inLims];    
    
    LGN_eval_w = 4;
    LGN_eval_lims_start = [0:2:16];
    nBlocks = length(LGN_eval_lims_start);
    OTCs = cell(nBlocks, nBlocks);
    for st_i = 1:nBlocks
        lims_x = LGN_eval_lims_start(st_i) + [0,LGN_eval_w];
        idx_inLims_x = find( ibetween(x_cents, lims_x) );
        
        for st_j = 1:nBlocks
            lims_y = LGN_eval_lims_start(st_j) + [0,LGN_eval_w];
            idx_inLims_y = find( ibetween(y_cents, lims_y) );
            
            LGN_ON_OFF_inLims = mHist2d_ON_OFF(idx_inLims_x, idx_inLims_y); %#ok<FNDSB>
            [ori_tc, R_ori_ph] = getOriTc(LGN_ON_OFF_inLims, all_oris', spatPeriod_pix, all_ph);
            OTCs{st_i, st_j} = ori_tc;
        end
        
    end
        
    %%
    all_OTCs = [OTCs{:}];
    figure(45);  clf;
    subplot(1,2,1); plot(all_oris, [OTCs{:}]);
    subplot(1,2,2); plot(all_oris, mean([OTCs{:}], 2)); hold on; plot(0,0);
    
    3;
    
%     idx_inLims = find( ibetween(x_cents, LGN_eval_lims) );    
%     LGN_ON_OFF_inLims = mHist2d_ON_OFF(idx_inLims, idx_inLims);
    
    
    
%     getGratingStimulusFrame(
    
%%
%     frame = generateGratingFrame(frameDims, 90, spatPeriod_pix, 0);
%     imagesc(1:nBins_inLims, 1:nBins_inLims, LGN_ON_OFF_inLims);

%%    

%     xx = 1:25;
%     yy = 1:25;
%     [XX, YY] = meshgrid(xx, yy);
%     k = 2*pi/spatPeriod_pix;
%     ZZ = gabor(1, 12.5, 12.5, 6, 6, pi, k, 0, 0, [XX(:), YY(:)]);
%     ZZ = reshape(ZZ, [25, 25]);
%     imagesc(ZZ);

    %%
    3;
    
    
%     gabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, C, XY) 
    
%     gabor(





end




function s = normVsqr(X, dim)
    if nargin == 1
        s = (sum(X(:).^2));  % faster just to use norm for this case, though.
    else
%         if length(size(X)) > 2
        s = (  sum(X.^2, dim)  );
    end
end

function [ori_tc, R_ori_ph] = getOriTc(RF, oris, spfs, phases)
    persistent allFrames allOris allPhs dims
    if isempty(allOris) || ~isequal(allOris, oris) || isempty(allPhs) || ~isequal(allPhs, phases) || ~isequal(size(RF), dims)
        allOris = oris;
        allPhs = phases;
        dims = size(RF);
        allFrames = cell(length(oris), length(phases));
    end
    
    frameDims = size(RF);
    nOris = length(oris);
%     nSpfs = length(spfs);
    nPhs = length(phases);

    R_ori_ph = zeros(nOris, nPhs);    
    for ori_i = 1:nOris
        for ph_i = 1:nPhs
            if ~isempty(allFrames{ori_i, ph_i}) 
                frame = allFrames{ori_i, ph_i};
            else
                frame = generateGratingFrame(frameDims, oris(ori_i), spfs, phases(ph_i) );
                allFrames{ori_i, ph_i} = frame;                
            end
            dt = sum((frame(:) .* RF(:)));
            R_ori_ph(ori_i, ph_i) = rectified(dt);            
%             figure(88); imagesc(frame); colormap gray; 
%             drawnow;
        end
    end    
    
    ori_tc = mean(R_ori_ph, 2);
        
end
    


function dists = createGaborParamDistsFile(filename)

    % need distributions for : sigma x, (sigma y/sigma x), k, diff in ori from MU
    % assume distributions for (theta (uniform), phi (uniform), 

%     diffFromMU = -b_val*log(1-rand(1,10000) );
%     Z = gabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, C, XY);

    dists = struct('sig_x', [], 'sig_yx_ratio', [], 'theta', [], 'k', [], 'k_cpd', [], 'phi', [], 'Dori_fromMU', []);
    dists.Dori_fromMU = struct('type', 'exponential', 'mu', 11.7);    
    dists.theta = struct('type', 'uniform', 'range', [0, pi]);
    dists.phi = struct('type', 'uniform', 'range', [0, 2*pi]);
    

    S = load('allMIDs.mat');
    show = 0;
    Scell = load('flashedGratingCells_GLFcuw8_phase.mat');
    file_Gids = [Scell.allCells.Gid];
    file_cellIds = [Scell.allCells.cellId];
    
    %%
    r_fit_th = .6;
    r_oe_th = 0.6;
    %%
    nMIDs = length(S.allGids);    
    all_r_oe = nan(1, nMIDs);
    for i = 1:nMIDs
        if ~isempty(S.allMIDs_even(i).MID) && ~isempty(S.allMIDs_odd(i).MID)
            all_r_oe(i) = pearsonR(S.allMIDs_even(i).MID, S.allMIDs_odd(i).MID);
        end        
    end
    all_r_fit = sqrt([S.allMIDs.rsqr]);
    %%
    idx_use = find(all_r_fit > r_fit_th & all_r_oe > r_oe_th);
     
    nUse = length(idx_use);
    
    all_params = cat(1, S.allMIDs(idx_use).gparams);    
%     all_params2 = cat(1, S.allMIDs.gparams);

%     diff in orientation: draw from 
        
%     all_theta = all_params(:,6)

%%
    if show
        i = 0;
        %%
%     i = i+1; idx = idx_use(i);
%     i = i+1; idx = small_k_idx(i);
        i = i+1; idx = big_sig_x(i);

        gid = S.allGids(idx);
        s = S.allMIDs(idx);
        figure(58); imagesc(s.MID);
        [mid_fit, xs, ys] = getMIDfitToGabor(gid, s.gparams);
        figure(59); imagesc(xs, ys, mid_fit);    
        p = s.gparams;
        pS = struct('sig_x', p(4), 'sig_y', p(5), 'k', p(7), 'lambda', 2*pi/p(7))
    end
    
    
    %% 1. sigma x
    sig_x = all_params(:,4);
    idx_sig_x_ok = sig_x < 5;    
    p_sig_x = lognfit(sig_x(idx_sig_x_ok));
    dists.sig_x = struct('type', 'lognormal', 'mu', p_sig_x(1), 'sigma', p_sig_x(2));
    
     if show        
        figure(41); clf; hold on;     
        [n,x] = hist(sig_x(idx_sig_x_ok), 50); 
        n = n/(sum(n) * diff(x(1:2)));
        bar(x, n, 1);

        [M,V] = lognstat(p_sig_x(1), p_sig_x(2));
        y = lognpdf( x, p_sig_x(1), p_sig_x(2));
        plot(x,y, 'r:');
%         drawVerticalLine(exp(p_sig_x(1)-p_sig_x(2)^2), 'color', 'g')

        % spatial frequency - log normal distribution.
%         mu = -.87; sigma = 0.51;
    end
    
    3;
    %% 2. ratio of sigma_y / sigma_x
    sig_y = all_params(:,5);
    yx_ratios = sig_x(idx_sig_x_ok) ./ sig_y(idx_sig_x_ok);    
    p_sig_yx_ratio = lognfit(yx_ratios);
    dists.sig_yx_ratio = struct('type', 'lognormal', 'mu', p_sig_yx_ratio(1), 'sigma', p_sig_yx_ratio(2));
    
    if show        
        figure(50); clf; hold on;     
        [n,x] = hist(yx_ratios, 50); 
        n = n/(sum(n) * diff(x(1:2)));
        bar(x, n, 1);
        
        y = lognpdf( x, p_sig_yx_ratio(1), p_sig_yx_ratio(2));
        plot(x,y, 'r:');

        % spatial frequency - log normal distribution.
%         mu = -.87; sigma = 0.51;
    end
   
     %% k : spatial frequency distribution
    k_gabor = nan(nUse, 1);
    k_profile = nan(nUse, 1);
    
    for i = 1:nUse
        idx = idx_use(i);
        Gid = S.allGids(idx);
        cellId = S.allCellIds(idx);
                        
        k_gabor(i) = S.allMIDs(idx).gparams(7);
        
        f_idx = find(file_Gids == Gid & file_cellIds == cellId, 1);
        if ~isempty(f_idx)
            k_profile(i) = Scell.allCells(f_idx).stats.tuningStats.spfStats_si.f_opt;
        end
    end
    
   %%

    useK_profile = 0;
    if useK_profile
        ks_cpd = k_profile;
        ks = ks_cpd * (2*pi);
    else        
        k_gabor_lim = 0.6;
        ks = k_gabor(k_gabor > k_gabor_lim);
        ks_cpd = ks / (2*pi);
    end
    ks = nonnans(ks);
    ks_cpd = nonnans(ks_cpd);
    p_k = lognfit(ks);
    p_k_cpd = lognfit(ks_cpd);
    dists.k = struct('type', 'lognormal', 'mu', p_k(1), 'sigma', p_k(2));
    dists.k_cpd = struct('type', 'lognormal', 'mu', p_k_cpd(1), 'sigma', p_k_cpd(2));
    
    if show
        figure(70);        
        plot(k_gabor, k_profile*2*pi, '.');
        for i = 1:2
            if i == 1
                ks_plot = ks; p = p_k;  txt = 'K (for gabor)';
            else
                ks_plot = ks_cpd; p = p_k_cpd; txt = 'K (cycles per degree)';
            end
        
            figure(71+i); clf; hold on;
            [n,x] = hist(ks_plot, 30); 
            n = n/(sum(n) * diff(x(1:2)));
            bar(x, n, 1);
            xlim([0, x(end)+diff(x(1:2))])
            
            [M,V] = lognstat(p(1), p(2));
            y = lognpdf( x, p(1), p(2));
            plot(x,y, 'r:');
            title({txt, sprintf('mean = %.2f, std = %.2f', M, sqrt(V))});
        end

        % spatial frequency - log normal distribution.
%         mu = -.87; sigma = 0.51;
    end
    
    
  
    
    %%
    save(filename, 'dists');

end

function X = randFromDist(dist, nVals)
   switch dist.type
       case 'uniform', 
           rng = dist.range;               
           X = rand(nVals, 1)*diff(rng) + rng(1);
       case 'exponential',
           X = exprnd(dist.mu, nVals, 1);           
       case 'lognormal',
           X = lognrnd(dist.mu, dist.sigma, nVals, 1);           
   end

end