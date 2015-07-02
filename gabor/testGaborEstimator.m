function testGaborEstimator

%     pLabels = {'A', 'mu_x', 'mu_y', 'sig_x', 'sig_y', '\theta', 'k', '\phi', 'c'};

    S = load('allMIDs.mat');
    [allMIDs, allGids, allCellIds, allRsqr] = deal(S.allMIDs, S.allGids, S.allCellIds, S.rsqr);    
    nRealCells = length(allCellIds);    
    idx_resort = ord(allRsqr, 'descend');
    allMIDs = allMIDs(idx_resort);
    allGids = allGids(idx_resort);
    allCellIds = allCellIds(idx_resort);
    
    [p_curFit, p_savedFit, p_savedEst] = deal([]);
    p_est_all = [];
    [curGid, curCellId, p_curFit, MID_curFit, cur_nFitMax, p_curEst] = deal([]);    
    [rsqr_savedEst, rsqr_savedFit, rsqr_curEst, rsqr_curFit] = deal(nan);
    [resn_savedEst, resn_savedFit, resn_curEst, resn_curFit] = deal(nan);
    [xs, ys, XY_C] = deal([]);
    [prevCellIdx, prevRealCell_tf, params_prev, params_cur] = deal(nan);
    curCellIdx = 32;
    curTransferMode = 'savedEst';
    [Xrange, Yrange] = deal([0 20]);
    Clims = [-1 1];
    [h_cont1a, h_cont1b, h_cont2a, h_cont2b] = deal(0);
    
    displayFunctions = true;    
    if displayFunctions
                
        Mspc = [.01 .01 .01];
        Nspc = [.15 .02 .02];
        figure(502); clf;
        [Z_orig, Z_orig_sm, Z_curEst, Z_curFit, Z_savedEst, Z_savedFit] = deal(zeros(2));
        
        h_im_ax(1) = subplotGap(2,3,1,1, Mspc, Nspc); h_im_orig = imagesc(Z_orig); axis equal tight ij; hold on;  h_tit1 = title('Orig');
        h_dummy_ax = axes('position', [0, .5, .001, .001]); h_col_im = imagesc(Z_orig);
        h_col = colorbar;
        set(h_col, 'position', [.02, .55, .05 .4])        
        h_im_ax(2) = subplotGap(2,3,2,1, Mspc, Nspc); h_im_orig_sm = imagesc(Z_orig_sm); axis equal tight ij; hold on; title('Smoothed'); 
        set([h_cont1a, h_cont1b, h_cont2a, h_cont2b], 'visible', 'on')
        h_im_ax(3) = subplotGap(2,3,1,2, Mspc, Nspc); h_im_est_saved = imagesc(Z_savedEst); axis equal tight ij; title('Saved Estimate');
        h_xlab_savedEst = xlabel(' ');
        h_im_ax(4) = subplotGap(2,3,2,2, Mspc, Nspc); h_im_fit_saved = imagesc(Z_savedFit); axis equal tight ij; title('Saved Fit');
        h_xlab_savedFit = xlabel(' ');
        h_im_ax(5) = subplotGap(2,3,1,3, Mspc, Nspc); h_im_est_cur = imagesc(Z_curEst); axis equal tight ij; title('Current Estimate');
        h_xlab_curEst = xlabel(' ');
        h_im_ax(6) = subplotGap(2,3,2,3, Mspc, Nspc); h_im_fit_cur = imagesc(Z_curFit); axis equal tight ij; title('Current Fit');
        h_xlab_curFit = xlabel(' ');
        

        whichContour_options = {'Saved estimate', 'Saved fit', 'Current estimate', 'Current fit', '[none]'};
        pos_func = @(i) [5 5+25*(i-1) 50, 20];
        uicontrol('style', 'pushbutton', 'string', 'Fit!', 'units', 'pixels', 'position', pos_func(1), 'callback', @doCurrentFit);
        uicontrol('style', 'pushbutton', 'string', 'Estimate!', 'units', 'pixels', 'position', pos_func(2), 'callback', @doCurrentEstimates);
        uicontrol('style', 'pushbutton', 'string', 'SAVE', 'units', 'pixels', 'position', pos_func(3), 'callback', @saveCurrentFit);
        uicontrol('style', 'pushbutton', 'string', 'Xfer', 'units', 'pixels', 'position', pos_func(4), 'callback', @transferFitParams);
        hContSel = uicontrol('style', 'popupmenu', 'string', whichContour_options, 'value', 2, 'units', 'pixels', 'position', pos_func(5), 'callback', @updateGaborContours);
        
        3;
%         set(h_im_ax, 'xtick', [], 'ytick', []);
        
    end
    
%     gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), fliplr(size(XY{1})) )';
    gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), size(XY{1}) );
    
    function updateParams(A, mu_x_rel, mu_y_rel, sig_x, sig_y, theta_deg, k, phi_deg, const, noiseLevel, smoothW, realCell_tf, cell_idx, showSaved, nFitMax)
            gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), size(XY{1}) );
                        
        curCellIdx = cell_idx;
        cell_idx_changed = (curCellIdx ~= prevCellIdx) || (prevRealCell_tf ~= realCell_tf);
        cur_nFitMax = nFitMax;
       
        params_cur = [A, nan, nan, sig_x, sig_y, deg2rad(theta_deg), k, deg2rad(phi_deg), const];
        params_cur(param_isLog) = 10.^(params_cur(param_isLog));
        
        
        paramsChanged = any(params_prev ~= params_cur);
        
        if ~realCell_tf
            theta = deg2rad(theta_deg);
            phi = deg2rad(phi_deg);
            
            nStds = 4;
%             [Xrange, Yrange] = getEffectiveSpaceOfGabor(p_orig, nStds);                
            
            Xrange = [0, 20];
            Yrange = [0, 20];
            
            Clims = [-A, A]+const;            
            
            
            mu_x = Xrange(1) + diff(Xrange)*mu_x_rel;
            mu_y = Yrange(1) + diff(Yrange)*mu_y_rel;
            params_cur(2:3) = [mu_x, mu_y];
            paramsChanged = any(params_prev ~= params_cur);
            p_orig = params_cur;
            
            %%%%    DRAW ORIGINAL GABOR FUNCTION.
            xs = linspace(Xrange(1), Xrange(2), 50);
            ys = linspace(Yrange(1), Yrange(2), 60);
            [xs_grid, ys_grid] = meshgrid(xs,ys);       
            XY_C = {xs_grid, ys_grid};

            p_orig = standardizeGaborP( p_orig );
            Z_orig = gaborFunc_XY(p_orig, XY_C);            
            
%             Z_orig = gaborFunc_XY(p_orig2, XY_C);
%             assert(isequal(Z_orig, Z_orig2));
%             assert( max( abs(Z_orig(:)-Z_orig2(:)) ) < 1e-10);

            set(h_im_est_saved, 'xdata', xs, 'ydata', ys, 'cdata', zeros(size(Z_savedEst)));
            set(h_im_fit_saved, 'xdata', xs, 'ydata', ys, 'cdata', zeros(size(Z_savedFit)));
            set([h_xlab_savedEst, h_xlab_savedFit], 'string', '');
            
            if noiseLevel > 0
                randn('state', 0);
                Z_orig = Z_orig + randn(size(Z_orig)) * max(abs(Z_orig(:))* noiseLevel );
            end
            title_str = 'Simulation';
                        
        elseif realCell_tf 
            
            title_str = sprintf('Gid = %d; cellId = %d', curGid, curCellId);            
            
            if cell_idx_changed
                fprintf('%s\n', title_str)
                p_curFit = [];
                p_est_all = [];

                curGid = allGids(cell_idx);
                curCellId = allCellIds(cell_idx);
                Z_orig = allMIDs{cell_idx};
                sd = siteDataFor(curGid);
                dsamp = 1;

                [xs, ys] = getStimulusXY(sd.stimulusInfo, dsamp);
                Xrange = [xs(1), xs(end)];
                Yrange = [ys(1), ys(end)];
    %             nx = size(Z_orig, 1);
    %             xs = [0:nx-1]*sd.stimulusInfo.degreesPerBlock*dsamp;            
    %             ys = xs;
                [xs_grid, ys_grid] = meshgrid(xs,ys); 
                XY_C = {xs_grid, ys_grid};
            end
            mu_x = Xrange(1) + diff(Xrange)*mu_x_rel;
            mu_y = Yrange(1) + diff(Yrange)*mu_y_rel;
            params_cur(2:3) = [mu_x, mu_y];
            paramsChanged = any(params_prev ~= params_cur);

            if cell_idx_changed
                const_est = median(Z_orig(:));
                L = max( abs(lims(Z_orig) - const_est));
                Clims = [-L, L]+const_est;
                
                Z_curEst = zeros(size(xs_grid));
                Z_curFit = zeros(size(xs_grid));
                
                if showSaved
                    [p_savedFit, rsqr_savedFit2, MID_curFit, s] = mid_getCellGaborParams(curGid, curCellId);
                    p_est_all = s.all_p_est(1,:);
                    p_savedEst = p_est_all(1,:);
                    Z_savedEst = gaborFunc_XY(p_savedEst, XY_C);            
                    Z_savedFit = MID_curFit;
                                        
                    [rsqr_savedEst, resn_savedEst] = getRsqrNorm(Z_orig, Z_savedEst, p_savedEst, XY_C);
                    [rsqr_savedFit, resn_savedFit] = getRsqrNorm(Z_orig, Z_savedFit, p_savedFit, XY_C);                    
                    assert(abs(rsqr_savedFit2 - rsqr_savedFit) < 1e-3);
                    set(h_im_est_saved, 'xdata', xs, 'ydata', ys, 'cdata', Z_savedEst);
                    set(h_im_fit_saved, 'xdata', xs, 'ydata', ys, 'cdata', Z_savedFit);
                    
                    set(h_xlab_savedEst, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_savedEst/resn_savedFit, rsqr_savedEst));
                    set(h_xlab_savedFit, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_savedFit/resn_savedFit, rsqr_savedFit));
                              
                    fprintf('Saved    : %s\n', p2str(p_savedFit(1,:)));
                    [A_saved, C_saved] = dealV( p_savedFit([1, 9]) );
                    Clims = [-A_saved, A_saved]+C_saved;
                    
%                     const_est = median(Z_savedEst(:));
%                     L = max( abs(lims(Z_savedEst) - const_est));
%                     Clims = 1.1*[-L, L]+const_est;                    
                end
                curTransferMode = 'savedFit';
                set(h_im_est_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curEst); 
                set(h_im_fit_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curFit); 

                set(h_xlab_curEst, 'string', ' ');
                set(h_xlab_curFit, 'string', ' ');
                
            end
        end

        
                
        if realCell_tf && paramsChanged && ~cell_idx_changed;
            p_curEst = params_cur;
            Z_curEst = gaborFunc_XY(p_curEst, XY_C);
            set(h_im_est_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curEst); 
           
            [rsqr_curEst, resn_curEst] = getRsqrNorm(Z_orig, Z_curEst);
            set(h_xlab_curEst, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_curEst/resn_savedFit, rsqr_curEst));                            
        end
        
        updateGaborContours;
        set(h_im_orig, 'xdata', xs, 'ydata', ys', 'cdata', Z_orig);
        
        if  (smoothW > 0) % (noiseLevel > 0) || realCell_tf
            circ_flag = 0;
            Z_orig_sm = gaussSmooth(gaussSmooth(Z_orig, smoothW , 1, circ_flag), smoothW, 2, circ_flag);
        else
            Z_orig_sm = Z_orig;
        end
                
        set(h_im_orig_sm, 'xdata', xs, 'ydata', ys', 'cdata', Z_orig_sm)
                
        set(h_tit1, 'string', title_str);
        if showSaved
%             h_im_est_cur
            

            
        end        
                



        %   L = max(abs([Z_orig(:); Z_curEst(:)]));                        
        set(h_im_ax, 'xlim', lims(xs), 'ylim', lims(ys));
        for i = 1:length(h_im_ax)
%             caxis(h_im_ax(i), Clims);
            caxis(h_im_ax(i), 'auto');
        end
        set(h_col_im, 'cdata', [Clims(:), Clims(:)]);

%             set(h_im_ax, 'clim', [-L, L]+const);
        set(h_im_orig, 'xdata', xs, 'ydata', ys, 'cdata', Z_orig);
        set(h_im_orig_sm, 'xdata', xs, 'ydata', ys, 'cdata', Z_orig_sm);

        3;
            
        prevRealCell_tf = realCell_tf;
        prevCellIdx = curCellIdx;
        params_prev = params_cur;
    end

    function doCurrentEstimates(~,~)

            alreadySmoothed_flag = 1;
            %%%%    ESTIMATE GABOR FUNCTION.        
            p_est_all =  estimateGaborParameters(xs, ys, Z_orig_sm, alreadySmoothed_flag);
            fprintf('Obtained %d estimates\n', size(p_est_all, 1));
            p_curEst = p_est_all;
%             if ~realCell_tf
%                 [A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const];   % for comparison
%                 [A0, mu_x0, mu_y0, sig_x0, sig_y0, theta0, k0, phi0, C0] = dealV(p_est);
%             end                
        
            Z_curEst = gaborFunc_XY(p_curEst(1,:), XY_C);
%             Rsqr = corr(Z_curFit(:), Z_orig(:))^2;        
    
            set(h_im_est_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curEst); 
            
            [rsqr_curEst, resn_curEst] = getRsqrNorm(Z_orig, Z_curEst);
            set(h_xlab_curEst, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_curEst/resn_savedFit, rsqr_curEst));                            

            curTransferMode = 'curEst';
            fprintf('Estimate : %s\n', p2str(p_curEst(1,:)));
        updateGaborContours;
    end

    function doCurrentFit(~,~)
        if isempty(p_curEst)
            return; 
        end
        %%%%    CALCULATE FITTED GABOR FUNCTION.

        [p_curFit, rsqr_curFit, Z_curFit, best_est_idx, ord_bestFits] = fitGaborParams(xs, ys, Z_orig, p_curEst, cur_nFitMax);
        p_curEst = p_curEst(best_est_idx, :);

        Z_curFit = gaborFunc_XY(p_curFit, XY_C);
        set(h_im_fit_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curFit); 
        [rsqr_curFit, resn_curFit] = getRsqrNorm(Z_orig, Z_curFit);
        set(h_xlab_curFit, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_curFit/resn_savedFit, rsqr_curFit));                            

        
        Z_curEst = gaborFunc_XY(p_curEst, XY_C);
        set(h_im_est_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curEst); 
        [rsqr_curEst, resn_curEst] = getRsqrNorm(Z_orig, Z_curEst);
        set(h_xlab_curEst, 'string', sprintf('rn = %3.3g (R^2 = %.2f)', resn_curEst/resn_savedFit, rsqr_curEst));                            
                
%         set(h_xlab_curFit, 'string', sprintf('R^2 = %.2f', rsqr_curFit));
        
        fprintf('Fit :      %s\n', p2str(p_curFit));
        
        updateGaborContours;
        curTransferMode = 'curFit';
    end


    function saveCurrentFit(~,~)
        mid_getCellGaborParams(curGid, curCellId, p_curFit, MID_curFit, rsqr_curFit);
        
    end

    function transferFitParams(~,~)
        switch curTransferMode
            case 'savedEst', p = p_savedEst; 
            case 'savedFit', p = p_savedFit;
            case 'curEst',   p = p_curEst(1,:);
            case 'curFit',   p = p_curFit;
        end
        p_curEst = p;
        
        p_set = p;
        p_set(param_isLog) = log10(p_set(param_isLog));
        [mu_x, mu_y] = dealV(p_set(2:3));        
        mu_x_rel = (mu_x-Xrange(1))/diff(Xrange);
        mu_y_rel = (mu_y-Yrange(1))/diff(Yrange);
        p_set(2:3) = [mu_x_rel, mu_y_rel];
        p_set(6) = rad2deg(p_set(6));
        p_set(8) = rad2deg(p_set(8));
        
        if ~isempty(p)
            manipulateSet(vHandles, paramNames, num2cell(p_set));
        end
        
%         Z_curEst = gaborFunc_XY(p_curEst, XY_C);
%         set(h_im_est_cur, 'xdata', xs, 'ydata', ys, 'cdata', Z_curEst);
        
        
    end

    function updateGaborContours(~,~)
        
        contour_option = whichContour_options{ get(hContSel, 'value') };
        
%         whichContour_options = {'Saved estimate', 'Saved fit', 'Current estimate', 'Current fit', '[none]'};
                
        
        switch contour_option
            case 'Saved estimate',   Z_cont = Z_savedEst;
            case 'Saved fit',        Z_cont = Z_savedFit;
            case 'Current estimate', Z_cont = Z_curEst;
            case 'Current fit',      Z_cont = Z_curFit;
            case '[none]',
                set([h_cont1a, h_cont2a, h_cont1b, h_cont2b], 'visible', 'off')
                return;
        end
                
        n = 3;
%         A_fit = p_curFit(1); C_fit = p_curFit(9);
        mn = min(Z_orig(:)); 
        mx = max(Z_orig(:)); 
        C0 = median(Z_orig(:));
        C_fit = C0;

        level_list_a = linspace(C_fit, mx, n+2);% [-A_fit/2, A_fit/2];
        level_list_a = level_list_a(2:end-1);
        level_list_b = linspace(mn, C_fit, n+2);% [-A_fit/2, A_fit/2];
        level_list_b = level_list_b(2:end-1);
                
%                 level_list_a = A_fit*[1/2]+C_fit;% [-A_fit/2, A_fit/2];
%                 level_list_b = A_fit*[-1/2]+C_fit;% [-A_fit/2, A_fit/2];        
        if diff(lims(Z_cont)) > 0

            if h_cont1a == 0;
                [~,h_cont1a] = contour(h_im_ax(1), Z_cont, level_list_a);
                [~,h_cont1b] = contour(h_im_ax(1), Z_cont, level_list_b);                
                [~,h_cont2a] = contour(h_im_ax(2), Z_cont, level_list_a);
                [~,h_cont2b] = contour(h_im_ax(2), Z_cont, level_list_b);
            end

            set([h_cont1a, h_cont2a], 'xdata', xs, 'ydata', ys, 'zdata', Z_cont, 'levelList', level_list_a, 'linecolor', 'k', 'visible', 'on');
            set([h_cont1b, h_cont2b], 'xdata', xs, 'ydata', ys, 'zdata', Z_cont, 'levelList', level_list_b, 'linecolor', 'w', 'visible', 'on');        
       

        end
        
        
%         figure(99); 
%         contour(xs, ys, Z_cont, [level_list_a, level_list_b]);
        
                                
        
    end

    paramNames = {'A', 'mu_x_rel', 'mu_y_rel', 'log_sigma_x', 'log_sigma_y', 'theta', 'log_k', 'phi', 'const'};
    param_isLog = logical([0     0        0       1           1         0        1     0       0]);
    A_range = [.001:.005:.8]; % 1
    rel_mu_range = [-.5:.01:1.5];
    log_sigma_range = [-1:.02:3.5];
    theta_range = [-360:1:360]; %  6
    log_k_range = [-2:.01:2];    
    phi_range = [-360:1:360];
    const_range = [-.01:.0001:.01];

    args = {{'A', A_range, A_range(10)}, {'mu_x_rel', rel_mu_range, .5}, {'mu_y_rel', rel_mu_range, .5}, ...
            {'log_sigma_x', log_sigma_range, 0.4}, {'log_sigma_y', log_sigma_range, 0.3}, ...
           {'theta', theta_range, 45}, {'log_k', log_k_range, 0}, {'phi', phi_range, 0}, {'const', const_range, 0},...
           {'noiseLevel', [0:.05:1]}, ... 
           {'smoothWidth', [0:.05:3]}, {'actualCells', [false, true]}, {'cell_idx', [1:nRealCells], curCellIdx}, ...
           {'showSaved', [true, false], true}, {'nEstimates', [1:100], 5} };            
              
    vHandles = manipulate(@updateParams, args, 'FigId', 600);


end

function [rsqr, resnorm, resnorm_wgt] = getRsqrNormWgt(Z1, Z2, p, XY)
    rsqr = pearsonR(Z1(:), Z2(:))^2;
    n = numel(Z1);
    dZsqr = (Z1(:)-Z2(:)) .^2;
    resnorm = sum( dZsqr )/n;
    
    [A, mu_x, mu_y, sig_x, sig_y, theta] = dealV(p(1:6));
    Z_wgt = gabor(A, mu_x, mu_y, sig_x, sig_y, theta, -100, 0, 0, XY);
    
    resnorm_wgt = wgt_sum( dZsqr, Z_wgt(:) )/n;
    
    
end

function s = p2str(p)
    [A, mu_x, mu_y, sig_x, sig_y, theta_rad, k, phi_rad, C] = dealV(p);
    theta_deg = rad2deg(theta_rad);
    phi_deg = rad2deg(phi_rad);
    s = sprintf('A = %.2f. mu = [%.2f, %.2f], sig = [%.2f, %.2f], theta = %.1f, k = %.2f, phi = %.2f, C = %.2g', ...
        A, mu_x, mu_y, sig_x, sig_y, theta_deg, k, phi_deg, C);

end

%{

            showSurfPlots = 0;
if showSurfPlots
            figure(501); clf;
            h_surf(1) = mesh(zeros(1,3),zeros(1,3),zeros(3)); hold on;
            h_surf(2) = mesh(zeros(1,3),zeros(1,3),zeros(3));
            h_surf(3) = mesh(zeros(1,3),zeros(1,3),zeros(3)); 
            uicontrol(501, 'style', 'checkbox', 'units', 'pixels', 'position', [0 0 20 20], 'value', 1, 'backgroundColor', 'b', 'userdata', 1, 'callback', @updateSurfVisibilities);
            uicontrol(501, 'style', 'checkbox', 'units', 'pixels', 'position', [20 0 20 20], 'value', 1, 'backgroundColor', 'g', 'userdata', 2, 'callback', @updateSurfVisibilities);
            uicontrol(501, 'style', 'checkbox', 'units', 'pixels', 'position', [40 0 20 20], 'value', 1, 'backgroundColor', 'r', 'userdata', 3, 'callback', @updateSurfVisibilities);
        end



            showSurfPlots = 0;
            if showSurfPlots
                set(h_surf(1), 'xdata', xs, 'ydata', ys, 'zdata', Z_orig, 'EdgeColor', 'b')
                set(h_surf(2), 'xdata', xs, 'ydata', ys, 'zdata', Z_curEst, 'EdgeColor', 'g', 'visible', 'on')
                set(h_surf(3), 'xdata', xs, 'ydata', ys, 'zdata', Z_curFit, 'EdgeColor', 'r', 'visible', 'off')
                view(h_ax1, rad2deg(theta0), 30);
            end


%         pCell = num2cell(p_orig); b0Cell = num2cell(p_est); bCell = num2cell(p_curFit);
%         A = {' ', pLabels{1:8} ;
%              'actual: ', pCell{1:8} ;   
%              'comptd: ', bCell{1:8} ;
%              'estimd: ', b0Cell{1:8} }
%              cycleThroughFigures([h_surf(1), h_surf(2), h_surf(3)], 1, .6, {'actual', 'computed', 'estimated'});
%             input('Press Any Key to continue');


    function updateSurfVisibilities(src,~)
        val = get(src, 'value');
        id = get(src, 'userdata');
        
        set(h_surf(id), 'visible', iff(val, 'on', 'off') );        
    end


%}