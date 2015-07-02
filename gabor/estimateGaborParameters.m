function all_p0 = estimateGaborParameters(xs, ys, zs, smoothedAlready_flag)  % xs and ys are vectors, zs is a matrix;

    smooth_w = 1;
    k_nearest = 2;
%     min_peak_ratio = .33;

    nContourLevels = 15;
    contour_RelHeight = .75;
    extrema_nStd_consider = 0.5;
    extrema_nStd_minForPair = 1;    
    fit1Dgabor = 1;

    show3D = 0;
    showIM = 0;
        showOnSameFigure = 1;
        show1Dgabor = 0;
        showEstimateIMs = 1;
        
    zs = double(zs);        

    doSmoothing = ~exist('smoothedAlready_flag', 'var') || isempty(smoothedAlready_flag);        
    if doSmoothing
        zs = gaussSmooth(gaussSmooth(zs, smooth_w, 1), smooth_w, 2);
    end
    Nx = length(xs);
    Ny = length(ys);
    
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    XY_C = {xs_grid, ys_grid};
    gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), size(XY{1}) );
    
    % estimate C: median of all pixels;
    C0 = median(zs(:));
    zs = zs - C0;
    z_std = std(zs(:));

    
    % find Contours
    C_mtx = contourc(xs, ys, zs, nContourLevels);
    [allContours, allLevels] = contourMtx2Cell(C_mtx);
    
    % Find local maxima & minima. These will the the starting points for the parameter estimates.
    [idx_locExtr_y, idx_locExtr_x, locExtr_polarity, locExtr_vals] = findLocalExtrema2D_Matlab('minmax', zs, 1);        
    if isempty(idx_locExtr_y)
        error('No minima/maxima found: try using higher smoothing');
    end
    
    % only use extrema that are at least 1 std deviation from baseline, and are not at the very edge.
    idx_use = find(   abs(locExtr_vals(:)) > z_std*extrema_nStd_consider ...
                & between(idx_locExtr_x, 1, Nx) & between(idx_locExtr_y, 1, Ny) );
    idx_locExtr_x = idx_locExtr_x(idx_use);
    idx_locExtr_y = idx_locExtr_y(idx_use);
    locExtr_polarity = locExtr_polarity(idx_use);
    locExtr_vals = locExtr_vals(idx_use);
    
    % remove duplicates (two extrema within the same contour) use the +1-std and -1std contours for
    % this purpose
    plus_contour_lvl_idx = find(allLevels > z_std*extrema_nStd_minForPair, 1, 'first');
    if isempty(plus_contour_lvl_idx)
        plus_contour_lvl_idx = length(allLevels);
    end        
    idx_peaks = find(locExtr_polarity == 1);
    peak_XY = [xs(idx_locExtr_x(idx_peaks)); ys(idx_locExtr_y(idx_peaks))];    
    peak_contour_idxs = idx_contourClosestToPoints(peak_XY, allContours{plus_contour_lvl_idx} ); % find which contour contains each pea        peak_contour_level_idxs        
    [uPeakContourIdxs, peakContourIdxLists] = uniqueList(peak_contour_idxs);        % find which peaks share the same contour
    highestInList = cellfun(@(cont_idxs)  indmax(locExtr_vals(idx_peaks(cont_idxs)) ), peakContourIdxLists, 'un', 0); 
    idx_peaks_use = cellfun(@(list, idx_choose) list(idx_choose), peakContourIdxLists, highestInList); % pick the highest peak in each contour
    idx_non_nans = ~isnan(uPeakContourIdxs);    
    idx_peaks = idx_peaks(idx_peaks_use(idx_non_nans));

    minus_contour_lvl_idx = find(allLevels < -z_std*extrema_nStd_minForPair, 1, 'last');
    if isempty(minus_contour_lvl_idx)
        minus_contour_lvl_idx = 1;
    end
    idx_troughs = find(locExtr_polarity == -1);
    trough_XY = [xs(idx_locExtr_x(idx_troughs)); ys(idx_locExtr_y(idx_troughs))];    
    trough_contour_idxs = idx_contourClosestToPoints(trough_XY, allContours{minus_contour_lvl_idx} ); % find which contour contains each trough
    [uTroughContourIdxs, troughContourIdxLists] = uniqueList(trough_contour_idxs);        % find which troughs share the same contour
    lowestInList = cellfun(@(cont_idxs)  indmin(locExtr_vals(idx_troughs(cont_idxs)) ), troughContourIdxLists, 'un', 0); 
    idx_troughs_use = cellfun(@(list, idx_choose) list(idx_choose), troughContourIdxLists, lowestInList); % pick the highest trough in each contour
    idx_non_nans = ~isnan(uTroughContourIdxs);    
    idx_troughs = idx_troughs(idx_troughs_use(idx_non_nans));    
    
    idx_unique = [idx_peaks; idx_troughs];
    
    idx_locExtr_x = idx_locExtr_x(idx_unique);
    idx_locExtr_y = idx_locExtr_y(idx_unique);
    locExtr_polarity = locExtr_polarity(idx_unique);
    locExtr_vals = locExtr_vals(idx_unique);
    3;
%     figure(14);
%     plot(idx_locExtr_y(ip), idx_locExtr_x(ip), 'k*')
%     plot(idx_locExtr_y(it), idx_locExtr_x(it), 'ko')
    3;

    if showIM
        figure(545); clf; 
        imagesc(xs, ys, zs); axis equal tight ij; hold on;
        for i = 1:length(idx_locExtr_x)
            b = iff(locExtr_polarity(i) == 1, '()', '[]');            
            text( xs(idx_locExtr_x(i)), ys(idx_locExtr_y(i)), sprintf('%s%d%s', b(1), i, b(2)), 'color', 'w', 'horiz', 'cent', 'vert', 'mid');
        end            
    end
    
    % select pairs of extrema to use.
    
        
    
    nExtrema = length(idx_unique);    
    
    allExtrPairs_mtx = false(nExtrema, nExtrema);
    extrema_pos = [idx_locExtr_x(:), idx_locExtr_y(:)]';
    idx_closest = kNearestNeighbors(extrema_pos, k_nearest);
    for i = 1:nExtrema
        for j = 1:nExtrema
            allExtrPairs_mtx(i,j) = (i == j) || ...  % ie. use this extremum by itself.
                ((locExtr_polarity(i) ~= locExtr_polarity(j)) && ... % or use it with an extremum
                any(idx_closest(:,i) == j)) &&  ...                      % close by of opposite polarity.      
                abs(locExtr_vals(i)) > abs(locExtr_vals(j)) && ... ;
                any(abs(locExtr_vals([i,j])) > z_std*extrema_nStd_minForPair);
%                 locExtr_vals 
%                 r = abs( (peakVal2)/(peakVal1) );
            
        end
    end
    
    nEstimates = nnz(allExtrPairs_mtx);  
%     nEstimates = min(nEstimates, 20);
    [idx_pair_i, idx_pair_j] = find(allExtrPairs_mtx);
    nGaborParams = 9;
    all_p0 = zeros(nEstimates, nGaborParams);
    resnorms = nan(1, nEstimates);
    est_M = floor(sqrt(nEstimates));
    est_N = ceil(nEstimates/est_M);
    zs_est = cell(1,nEstimates);
    est_ok = true(1, nEstimates);
    
%     printProgress = 1;

    progressBar('init-', nEstimates, 30);
    for est_i = 1:nEstimates
        i = idx_pair_i(est_i);
        j = idx_pair_j(est_i);
            
        useTwoPeaks = i~=j;        
        
        xi = idx_locExtr_x(i);
        yi = idx_locExtr_y(i);                
        peakVal1 = zs(yi, xi);        
        assert(peakVal1 == locExtr_vals(i));


        peak1_pos = [xs(xi); ys(yi)];
%         isPeakAtEdge = any(x1i == [1, Nx]) || any(y1i == [1, Ny]);
        peak1_contour_orig = getContourWithExtremaAtRelHeight(allLevels, allContours, peak1_pos, peakVal1, contour_RelHeight);
        peak1_contour = removeIdenticalEnds(peak1_contour_orig);
        peak1_contour_COM = mean(peak1_contour, 2);
        
        if useTwoPeaks 
            % try use a pair of extrema. The pivot point is just the COM of the first one, and
            % theta0 is perpendicular to the line drawn from the first COM to the second.
            xj = idx_locExtr_x(j);
            yj = idx_locExtr_y(j);            
            peak2_pos = [xs(xj); ys(yj)];
            peakVal2 = zs(yj, xj);
            assert(peakVal2 == locExtr_vals(j));
            
%             r = abs( (peakVal2)/(peakVal1) );
%             if (r < min_peak_ratio) % second peak is not significant (is at least 1/3 as big as first peak)
%                 all_p0(est_i,:) = nan;
%                 continue;
%             end

            % otherwise, second peak is also significant -- use it with the first peak to estimate
            % the pivot point.            
            peak2_contour_orig = getContourWithExtremaAtRelHeight(allLevels, allContours, peak2_pos, peakVal2, contour_RelHeight );
            peak2_contour = removeIdenticalEnds( peak2_contour_orig );

            peak2_contour_COM = mean(peak2_contour, 2);    
            
            dx = peak1_contour_COM(1) - peak2_contour_COM(1);
            dy = peak1_contour_COM(2) - peak2_contour_COM(2);

            theta0 = atan2(dy, dx);
            pivot_X = peak1_contour_COM;
               
            point_in_y_dir = pivot_X + [cos(theta0+pi/2); sin(theta0+pi/2)];

        else   % just use the first peak.
            % Use one extremum. The pivot point is just the center of the first one, and
            % theta0 is estimated by looking at the principal components of the contour.

            % get the principal components of the contour of peak1.                         
            peak_cov = cov(peak1_contour');
            [peakEigVec,peakEigVal] = eig(peak_cov);
            v_pc1 = peakEigVec(:, indmax(diag(peakEigVal)));
%             v_pc2 = peakEigVec(:, indmin(diag(peakEigVal)));        

            
            theta0 = - atan2(v_pc1(1), v_pc1(2)) ;        
            pivot_X = peak1_contour_COM; 
            
            point_in_y_dir = pivot_X + v_pc1;
        end

        % we estimate sig_y0 by using the contour        
        [Q1, Q2] = itpContourPointsClosestToLine(peak1_contour_orig, peak1_contour_COM, point_in_y_dir);
        
        % interpolate to find the points exactly along the PC1 line.
        xy_i = [Q1, Q2, peak1_contour_COM];
        Zs_q = interp2(xs_grid,ys_grid,zs, xy_i(1,:), xy_i(2,:) );
        [Z_q1, Z_q2, Z_com] = dealV(Zs_q);

        dx1 = norm(Q1-peak1_contour_COM); 
        Z_peak1 = iff( abs(Z_com) > abs(Z_q1), Z_com, peakVal1);
        sig_est1 = dx1 / sqrt( 2*log(Z_peak1/Z_q1) );

        dx2 = norm(Q2-peak1_contour_COM);        
        Z_peak2 = iff( abs(Z_com) > abs(Z_q2), Z_com, peakVal1);
        sig_est2 = dx2 / sqrt( 2*log(Z_peak2/Z_q2) );
        
        sig_ests = [sig_est1 sig_est2];
        sig_y0 = mean(sig_ests(~isnan(sig_ests) & isreal(sig_ests) )); % one of these could be nan if one of the edges is outside the plot
        if isnan(sig_y0)
            dists = [dx1, dx2];
            sig_y0 = mean(dists(~isinf(dists)));
        end
        
        if isnan(sig_y0)            
            sig_y0 = 1;            
        end
               
        
        % once we have a pivot_X and an estimate for theta0, we shift coordinate axes to what we
        % estimate is the original gabor frame. 
        % (ie. x and y shifted by pivot_X and rotated by theta0)

        dx = diff(xs(1:2));  
        xs_ext = [xs - dx*(Nx), xs]';
        X0 = [xs(1); ys(1)];
        X_xy = [xs_ext(:), ys(1)*ones(length(xs_ext), 1)]';
        rotatedToXframe = shiftRotateShift(X_xy, -X0, theta0, pivot_X);
        
        X_xi = rotatedToXframe(1,:); % x values of xs (in new coordinates)
        X_yi = rotatedToXframe(2,:); % y values of xs (in new coordinates)
        X_zi = interp2(xs_grid,ys_grid, zs, X_xi, X_yi); % interp Z values of y coordinates
        X_nonnans = find(~isnan(X_zi));                  % only keep the ones inside the valid xy range
        X_xi = X_xi(X_nonnans); X_yi = X_yi(X_nonnans); X_zi = X_zi(X_nonnans); 
        x_ok = xs_ext(X_nonnans)';

        
        
        if showIM
            dy = diff(ys(1:2));  
            ys_ext = [ys - dy*(Ny), ys]';
            Y_xy = [xs(1)*ones(length(ys_ext), 1), ys_ext(:)]';
            rotatedToYframe = shiftRotateShift(Y_xy, -X0, theta0, pivot_X);

            % I used to use the Y_zi values to estimate the sigma_y0, but this is sensitive to noise, so
            % i now use contour2 for this purpose
            Y_xi = rotatedToYframe(1,:);  Y_yi = rotatedToYframe(2,:); % ys rotated to new coordinates
            Y_zi = interp2(xs_grid,ys_grid, zs, Y_xi, Y_yi); % interpolated Z values of y coordinates
            Y_nonnans = find(~isnan(Y_zi));                  % only keep the ones inside the valid xy range
            Y_xi = Y_xi(Y_nonnans); Y_yi = Y_yi(Y_nonnans); % Y_zi = Y_zi(Y_nonnans); 
%             y_ok = ys_ext(Y_nonnans)';
        end            
                  
        
        % Estimate parameters of 1D gabor function (in rotated x' frame)
        [A0, mu_est, sig_x0, k0, phi0] = estimate1DgaborParams(x_ok, X_zi, fit1Dgabor);        
        
        
        idx_Xi_mu = indmin(abs(mu_est - x_ok));  % rotate 'mu' back to original coordinates    
        MU = [rotatedToXframe(:, X_nonnans(idx_Xi_mu))];
        mu_x0 = MU(1);
        mu_y0 = MU(2);
        
        
        p0 = [A0, mu_x0, mu_y0, sig_x0, sig_y0, theta0, k0, phi0, C0];
        p0 = standardizeGaborP( p0 );
        
        % calculate the norm of the residuals from this estimate;
        zs_est{est_i} = gaborFunc_XY(p0, XY_C);        
        resnorms(est_i) = norm(zs_est{est_i}(:)-zs(:));
        gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), size(XY{1}) );
                
        
        if showIM            
            if showOnSameFigure
                figure(546); 
                subplot(est_M, est_N, est_i);
            else                
                figure(545 + est_i); clf; 
                if show1Dgabor
                    subplot(5,1,1:4);                
                end
            end
            imagesc(xs, ys, zs); axis equal tight ij; hold on;
            set(gca, 'xtick', [], 'ytick', []);


            plot(X_xi, X_yi, 'b-')
            plot(Y_xi, Y_yi, 'r-')
            plot(pivot_X(1), pivot_X(2), 'go', 'linewidth', 3, 'markersize', 15);

            plot(peak1_contour_orig(1,:), peak1_contour_orig(2,:), 'w-')
            plot(peak1_contour_COM(1), peak1_contour_COM(2), 'w*')            
            
            plot([Q1(1) Q2(1)], [Q1(2) Q2(2)], 'k*');
            plot(mu_x0, mu_y0, 'ko', 'markersize', 10, 'linewidth', 2);            
            
            if useTwoPeaks
                title(sprintf('(%d) %d, %d', est_i, i, j))
                plot(peak2_contour_orig(1,:), peak2_contour_orig(2,:), 'w:')
                plot(peak2_contour_COM(1), peak2_contour_COM(2), 'w+')                
                
            else
                title(sprintf('(%d) %d', est_i, i));
                com = peak1_contour_COM;
                v1 = peakEigVec(:,1)*peakEigVal(1,1);
                v2 = peakEigVec(:,2)*peakEigVal(2,2);

                plot([0, v1(1)]+com(1), [0, v1(2)]+com(2), 'kd-');
                plot([0, v2(1)]+com(1), [0, v2(2)]+com(2), 'ko-');            

    %             plot(peak1_contour(1, idx_closest), peak1_contour(2, idx_closest), 'wp');

                3;
                
            end
            3;
            if ~showOnSameFigure && show1Dgabor
                subplot(5,1,5); plot(x_ok, X_zi); hold on;
                drawVerticalLine(mu_est);
                plot(x_ok, gabor1D(A0, mu_est, sig_x0, k0, phi0, x_ok), 'r:');
                xlim(x_ok([1, end]))
            end


            
        end


        
        est_ok(est_i) = is_param_ok(p0);

        all_p0(est_i, :) = p0;
        
%         if printProgress
%             fprintf('%d/%d',  nEstimates)
%         end
        progressBar(est_i);
    end
    
    if showIM && showEstimateIMs
        figure(746); clf(746); 
        ord_ests = ord(resnorms);
        for est_i = 1:nEstimates
            i = ord_ests(est_i);
            if showOnSameFigure
                figure(746); 
                subplot(est_M, est_N, est_i);
            else                                    
                figure(745 + est_i); clf; 
            end
            imagesc(xs, ys, zs_est{i}); 
            set(gca, 'xtick', [], 'ytick', []);
            title(sprintf('(%d) %.5g', i, resnorms(i)))                
            axis equal tight ij; 
        end
    end

    % remove estimates that were not ok.
    all_p0 = all_p0(est_ok,:);
    resnorms = resnorms(est_ok);
        
    % sort the p0 estimates by their residual-norms (smallest residuals first):
    idx_bestEstimates = ord(resnorms);  % note: nans are placed at the end, which is ok.
    all_p0 = all_p0(idx_bestEstimates, :);    
    
end


function tf = is_param_ok(p)
    tf = 1;
    [A0, mu_x0, mu_y0, sig_x0, sig_y0, theta0, k0, phi0, C0] = dealV(p);
    
    if any( imag(p) ~= 0) || any(isnan(p)) || any(isinf(p)) 
        tf = 0;
    end
    if any([A0 sig_x0 sig_y0 k0] <= 0)
        tf = 0;
    end

end

%             if show3D
%                 stem3(mu_x0, mu_y0, peakVal1*1.5, 'ro', 'linewidth', 3);
%                 stem3(mu_x0, mu_y0, -peakVal1*1.5, 'ro', 'linewidth', 3);
%                 stem3(mu_x_opp, mu_y_opp, -peakVal1*1.5, 'go', 'linewidth', 3);
%                 stem3(mu_x_opp, mu_y_opp, peakVal1*1.5, 'go', 'linewidth', 3);
% 
%     %             plot(mu_x0, mu_y0, 'wo')
%             end


%         if show3D
%             stem3(rotatedToXframe(1,:), rotatedToXframe(2,:), ones(length(rotatedToXframe),1) )
%             stem3(rotatedToYframe(1,:), rotatedToYframe(2,:), ones(length(rotatedToYframe),1), 'g' )
%         end

%{
        if show
%             figure(main_fig_id);
%             if dbug_theta, set(thetaStems, 'Visible', 'off'); end
            hh = stem3(X_xi, X_yi, X_zi*(1.1), 'ro');
            stem3(X_xi(1), X_yi(1), X_zi(1)*(1.5), 'b', 'fill', 'MarkerSize', 10);

            stem3(Y_xi, Y_yi, Y_zi*(1.01), 'g');
            stem3(Y_xi(1), Y_yi(1), Y_zi(1)*(1.5), 'b', 'fill', 'MarkerSize', 10);
        end
%}

function X = shiftRotateShift(X, Xc1, theta, Xc2)  
    X = bsxfun(@plus, X, Xc1);
    X = rotationMatrix(theta) * (X);
    X = bsxfun(@plus, X, Xc2);    
end



function peak_contour = getContourWithExtremaAtRelHeight(allLevels, allContours, peak_pos, peak_val, relHeight)
        
    desiredHeight = peak_val * relHeight;
    level_idx = indmin( abs(allLevels - desiredHeight) );
    
    idx_withExtremum = idx_contourContainingPoint(peak_pos, allContours{level_idx});
    
    if isempty(idx_withExtremum)    
        idx_withExtremum = idx_contourClosestToPoints(peak_pos, allContours{level_idx});
    end
        
    peak_contour = allContours{level_idx}{idx_withExtremum};        
    
end



%{
function peak_contour = getContourWithExtremaAtRelHeight(allLevels, allContours, peak_pos, peak_val, extremumType, relHeight)
    
    if (extremumType == 1) || strcmp(extremumType, 'max')
        level_idx0 = length(allContours);        
        lvl_inc = -1;
        cmpFun = @le;
    elseif (extremumType == -1) || strcmp(extremumType, 'min')
        level_idx0 = 1;
        lvl_inc = +1;
        cmpFun = @ge;
    end    

    % we start from the highest (lowest) contour, and continue down (up) until we get to the desired
    % relative height (or until we can't find the peak in the contours anymore)
    
    % contour with most extreme (absolute) peak:
    desiredHeight = peak_val * relHeight;
    level_idx = level_idx0;
    curLevel = allLevels(level_idx);    
%     idx_withExtremum = point_contourIdx(peak_pos, allContours{level_idx}, isPeakAtEdge);
    idx_withExtremum = idx_contourClosestToPoints(peak_pos, allContours{level_idx});
    
    while ~cmpFun(curLevel, desiredHeight)
        level_idx = level_idx + lvl_inc;                
        curLevel = allLevels(level_idx);
        
%         idx_withExtremum = point_contourIdx(peak_pos, allContours{level_idx}, isPeakAtEdge);
        idx_withExtremum = idx_contourClosestToPoints(peak_pos, allContours{level_idx});
        if isempty(idx_withExtremum)
            assert(level_idx ~= level_idx0);
            level_idx = level_idx - lvl_inc;
            break;
        end        
    end

    if isempty(idx_withExtremum)
        figure(54); clf;
        plot(peak_pos(1), peak_pos(2), 'r*'); hold on;
        for i = 1:length(allContours{level_idx})
            xy = allContours{level_idx}{i};
            plot(xy(1,:), xy(2,:), [color_s(i) '.-']);
        end
        keyboard;
    end
    assert(length(idx_withExtremum) == 1);    
    
    peak_contour = allContours{level_idx}{idx_withExtremum};        
    
%     idx_withExtremum = point_contourIdx(peak_pos, allContours{level_idx});
%     assert(length(idx_withExtremum) == 1);
%     peak_contour = allContours{level_idx}{idx_withExtremum};
    

end
%}
%{
function idx = idx_contourContainingPoint(P, contours, isPeakAtEdge, f)
    nContours = length(contours);
    inAnyContours = false(1, nContours);
    if nargin < 4
        f = 1e-5;
    end
    for ci = 1:nContours        
        P_ci = P;
        if isPeakAtEdge  % adju
            com_ci = mean(contours{ci}, 2);
            v = com_ci - P;
            P_ci = P + v*(f);                    
        end
        inAnyContours(ci) = isPointInPolygon(P_ci, contours{ci});        
    end

%     inAnyContours = cellfun(@(poly) isPointInPolygon(P, poly), contours);
    idx = find(inAnyContours);
    
    
end
%}

function idx = idx_contourContainingPoint(P, contours)
    
    inAnyContours = cellfun(@(poly) isPointInPolygon(P, poly), contours);
    idx = find(inAnyContours);    
    
end


function idxs = idx_contourClosestToPoints(P, contours, nanIfEdge_flag)
    contourCOMs = cellfun(@(c) mean(c, 2), contours, 'un', 0);    
    idxs = zeros(1, size(P,2));
    returnNanIfAtEdge = (nargin > 2) && ~isempty(nanIfEdge_flag);
    for p_i = 1:size(P,2)            
        distFromCOMsToPoint_i = cellfun(@(com) norm(com-P(:,p_i)), contourCOMs);
        cidx = indmin(distFromCOMsToPoint_i);
        if returnNanIfAtEdge && ( ~all( contours{cidx}(:,1) == contours{cidx}(:,end) ) )
            idxs(p_i) = nan;
        else
            idxs(p_i) = cidx;
        end                
    end
    
end


function [idxs, lvl_idxs] = idx_closedContourClosestToPoints(P, contours, contour_lvl_idx_start, lvl_inc)
    nPoints = size(P,2);
    idxs = zeros(1, nPoints);
        
    
    returnNanIfAtEdge = (nargin > 2) && ~isempty(nanIfEdge_flag);
    for p_i = 1:nPoints
        lvl_i = contour_lvl_idx_start;
        contourCOMs_lvl_i = cellfun(@(c) mean(c, 2), contours, 'un', 0);    
                
        distFromCOMsToPoint_i = cellfun(@(com) norm(com-P(:,p_i)), contourCOMs);
        cidx = indmin(distFromCOMsToPoint_i);
        if returnNanIfAtEdge && ( ~all( contours{cidx}(:,1) == contours{cidx}(:,end) ) )
            idxs(p_i) = nan;
        else
            idxs(p_i) = cidx;
        end
        
        
    end
    

end


function contour1 = removeIdenticalEnds( contour1 )
    if all( contour1(:,1) == contour1(:,end) )
        contour1 = contour1(:,1:end-1);
    end
end

function [X_int1, X_int2] = itpContourPointsClosestToLine(C, p1, p2)
    
    % find which points in the contour are closest to the line defined by the two points [p1, p2];
%     if ~all( C(:,1) == C(:,end) ) % end point is the same as the first point
%         C = [C, C(:,end)];         
%     end
    L0 = {p1, p2};
    nPoints = size(C,2);
%     Q = zeros(size(C));
    ccws = zeros(1, nPoints);
    for p_i = 1:nPoints        
        ccws(p_i) = ccw(p1, p2, C(:,p_i));
    end
    idx_crossings = find(diff(ccws) ~= 0);
    
    if isempty(idx_crossings)
        X_int1 = inf(2,1);
        X_int2 = inf(2,1);
        return;
        
    elseif nnz(idx_crossings) >= 1
        idx1 = idx_crossings(1);
        idx1_opp = idx1+1;
        
        line1 = {C(:,idx1), C(:,idx1_opp)};
        [~, X_int1] = doLinesIntersect(line1, L0);
    end
    
    if nnz(idx_crossings) >= 2
        idx2 = idx_crossings(2);
        idx2_opp = idx2+1;
        
        line2 = {C(:,idx2), C(:,idx2_opp)};
        [~, X_int2] = doLinesIntersect(line2, L0);
    else
        X_int2 = inf(2,1);
        return;         
    end
                
    %{ 
        %visualize: 
        figure(93); clf; hold on; 
        plot(C(1,:), C(2,:), 'b.-');  axis(axis);
        dx = p2(1)-p1(1); dy = p2(2)-p1(2)
        plot(p1(1)+dx*[-10, 10], p1(2)+dy*[-10, 10], 'rs-')
        plot(C(1,[idx1 idx1_opp]), C(2,[idx1 idx1_opp]), 'gs-');  plot(X_int1(1), X_int1(2), 'g*')
        plot(C(1,[idx2 idx2_opp]), C(2,[idx2 idx2_opp]), 'rs-');  plot(X_int2(1), X_int2(2), 'r*')    
    %}

end

%{
    line1 = {C(:,idx1), C(:,idx1_opp)};
    [~, X_int1] = doLinesIntersect(line1, L0);
        
    idx2 = idx_sort( find( circDist( idx_sort, idx1, nPoints) > nPoints/4, 1)); % the point on the other side of the contour
    ccw2      = ccw(p1, p2, C(:,idx2), 1);
    idx2_left  = mod((idx2-1)-1, nPoints)+1; ccw2_left  = ccw(p1, p2, C(:,idx2_left), 1);
    idx2_right = mod((idx2-1)+1, nPoints)+1; ccw2_right = ccw(p1, p2, C(:,idx2_right), 1);
    if ccw2 == -ccw2_left
        idx2_opp = idx2_left;  % the point on the other side of the line.
    elseif ccw2 == -ccw2_right
        idx2_opp = idx2_right;
    else
        X_int2 = inf(2,1);
        return;
    end
    
    line2 = {C(:,idx2), C(:,idx2_opp)};
    [~, X_int2] = doLinesIntersect(line2, L0);

    
    %}

%{
%     peak_COM = getCenterOfMass(xs, ys, zs, idx_peak1(2),idx_peak1(1), 0);
%     peak_opp_COM = getCenterOfMass(xs, ys, zs, idx_peak2(2),idx_peak2(1), 0);

function COM = getCenterOfMass(X, Y, Z, xi, yi, n)
    
    [Ny, Nx] = size(Z);
    
    x1_idx = max(xi-n, 0);  x2_idx = min(xi+n, Nx);
    y1_idx = max(yi-n, 0);  y2_idx = min(yi+n, Ny);
    
    idx_around_peak_x = x1_idx:x2_idx;
    idx_around_peak_y = y1_idx:y2_idx;
    
    [idx_peak_x_grid, idx_peak_y_grid] = meshgrid(idx_around_peak_x, idx_around_peak_y);
    z_aroundPeak = Z(idx_around_peak_x, idx_around_peak_y);
    
    COM_x = sum(X(idx_peak_x_grid(:)) .* z_aroundPeak(:)')/sum( z_aroundPeak(:));
    COM_y = sum(Y(idx_peak_y_grid(:)) .* z_aroundPeak(:)')/sum( z_aroundPeak(:));

    COM = [COM_x; COM_y];

end
%}


%{

    %find X_center:     %  integrate half of square of function,
%     pow = 2;
%     abs_smoothed_zs = abs(  smooth2(zs, round(length(xs)/20) ) ) .^ pow;
%     [xi] = quadCenter(sum(abs_smoothed_zs,1));
%     [yi] = quadCenter(sum(abs_smoothed_zs,2));
%     xi = 33;
%     yi = 19;
%     [tmp, X_center] = maxElement( abs_smoothed_zs );
%     X_center = fliplr(X_center)';
%     [xi, yi] = elements(X_center);
%     X_center = [xs(xi); ys(yi)] - Xc;               
%     X_center = [40;25];
%}

    %{
%     x1_idx = max(idx_peak1(2)-1, 0);  x2_idx = min(idx_peak1(2)+1, Nx);
%     y1_idx = max(idx_peak1(1)-1, 0);  y2_idx = min(idx_peak1(1)+1, Ny);
%     
%     idx_around_peak_x = x1_idx:x2_idx;
%     idx_around_peak_y = y1_idx:y2_idx;
%     
%     [idx_peak_x_grid, idx_peak_y_grid] = meshgrid(idx_around_peak_x, idx_around_peak_y);
%     z_aroundPeak = zs(idx_around_peak_x, idx_around_peak_y);
%     
%     peak_COM_x = sum(xs(idx_peak_x_grid(:)) .* z_aroundPeak(:)')/sum( z_aroundPeak(:));
%     peak_COM_y = sum(ys(idx_peak_y_grid(:)) .* z_aroundPeak(:)')/sum( z_aroundPeak(:));
%}

%{
    %     shiftRotateShift(X_xy, X0, theta0, X_peak+X0);
%     rotatedToXframe = rotateAndShift(X_xy, pi/30, -(X_peak-X0));
%     rotatedToYframe = rotateAndShift(Y_xy, pi/30, -(X_peak-X0));

%     rotatedToXframe = shiftAndRotate(X_xy, -(X_peak-X0), 0);
%     rotatedToYframe = shiftAndRotate(Y_xy, -(X_peak-X0), 0);


function X = shiftAndRotate(X, Xc, theta)  
    X = bsxfun(@minus, X, Xc);
    X = rotationMatrix(-theta) * (X);
end

function X = rotateAndShift(X, theta, Xc)  
    X = rotationMatrix(-theta) * (X);
    X = bsxfun(@minus, X, Xc);

end

%}


%{ 
        % I used to use the Y_zi values to estimate the sigma_y0, but this is sensitive to noise, so
        i now use contour2 for this purpose
        Y_xi = rotatedToYframe(1,:);  Y_yi = rotatedToYframe(2,:); % ys rotated to new coordinates
        Y_zi = interp2(xs_grid,ys_grid, zs, Y_xi, Y_yi); % interpolated Z values of y coordinates
        Y_nonnans = find(~isnan(Y_zi));                  % only keep the ones inside the valid xy range
        Y_xi = Y_xi(Y_nonnans); Y_yi = Y_yi(Y_nonnans); Y_zi = Y_zi(Y_nonnans); 
        y_ok = ys_ext(Y_nonnans)';
%}


 %{     
    OLD (LESS ROBUST) ESTIMATION ALGORITHM
    %%% 1. ESTIMATE FOR center of gabor
    X_center = estimateGaborCenter(xs, ys, zs);

    %%% 2. ESTIMATE FOR const
    const = mean(zs(:));                            
    zs = zs - const;

    %%% 3. ESTIMATE FOR A
    A = max(abs(zs(:)));
    
            if dbug
                main_fig_id = 99;
                freq_fig_id = 100;
                zmax = 1.5* maxElements(zs);
                [xs_grid, ys_grid] = meshgrid(xs, ys);
                figure(main_fig_id); clf;
                surf(xs_grid, ys_grid, zs); hold on
                xlabel('x'); ylabel('y'); zlabel('z');
                colormap('gray');
                stem3(X_center(1) + Xc(1), X_center(2) + Xc(2), 1.2*zmax, 'r')
            end
    

    %%% 5. GET THETA (& LAMBDA FOR K)
    [theta, lambda] = findThetaLambdaWithStrongestOscillations(xs, ys, zs, X_center);

    %%% 5. ESTIMATE FOR k.
    k = 2*pi/lambda;
    
    %%% 6, 7.  ESTIMATE FOR sig_x and sig_y
    rotatedToXframe = rotateAndShift(Xs, theta, -(X_center+Xc));
    rotatedToYframe = rotateAndShift(Xs, theta + pi/2, -(X_center+Xc));
        
    X_xi = rotatedToXframe(1,:);  X_yi = rotatedToXframe(2,:);
    Y_xi = rotatedToYframe(1,:);  Y_yi = rotatedToYframe(2,:);

    X_zi = interp2(xs_grid,ys_grid,zs, X_xi, X_yi);
    Y_zi = interp2(xs_grid,ys_grid,zs, Y_xi, Y_yi);
    
    X_nonnans = ~isnan(X_zi);    Y_nonnans = ~isnan(Y_zi);
    X_xi = X_xi(X_nonnans); X_yi = X_yi(X_nonnans); X_zi = X_zi(X_nonnans); 
    Y_xi = Y_xi(Y_nonnans); Y_yi = Y_yi(Y_nonnans); Y_zi = Y_zi(Y_nonnans); 
    
            if dbug
                figure(main_fig_id);
                if dbug_theta, set(thetaStems, 'Visible', 'off'); end
                stem3(X_xi, X_yi, X_zi*(1.01), 'r');
                stem3(X_xi(1), X_yi(1), X_zi(1)*(1.5), 'b', 'fill', 'MarkerSize', 10);

                stem3(Y_xi, Y_yi, Y_zi*(1.01), 'g');
                stem3(Y_xi(1), Y_yi(1), Y_zi(1)*(1.5), 'b', 'fill', 'MarkerSize', 10);
            end

    sig_x = distribStd(dx, abs(X_zi) );    
	sig_y = distribStd(dx, abs(Y_zi) );    
    
    %%% 8. ESTIMATE FOR phi
    spacedXs = (0:length(X_zi)-1) * dx;
    x0 = norm(X_center+Xc - [X_xi(1); X_yi(1)]);
    phi = getPhaseOfCosine(spacedXs, X_zi, k, x0);       
%     if k < dx
%         phi = 0;
%     end
        
    if A < 0  % convention: A > 0;
        A = -A;
        phi = mod(phi + pi, 2*pi);
    end

    beta0 = [A;  X_center(1); sig_x; X_center(2); sig_y; k; phi; theta; const];
    %}

%         % find 2 most extreme values in Z
%         [maxVal_abs, idx_peak1] = maxElement(abs(zs)); 
%         peakVal1 = zs(idx_peak1(1), idx_peak1(2));    
            
%             [maxVal_opp, idx_peak2] = maxElement( -sign(peakVal1) * zs);
%             peakVal2 = zs(idx_peak2(1), idx_peak2(2));            
%             peak2 = [xs(idx_peak2(2)); ys(idx_peak2(1))];

%             idxFurthestFromCOM1 = indmax( normV( bsxfun(@minus, peak1_contour, peak1_contour_COM), 1) );
%             FurthestPoint1 = peak1_contour(:,idxFurthestFromCOM1);
%             dx_2 = FurthestPoint1(1) - peak1_contour_COM(1);
%             dy_2 = FurthestPoint1(2) - peak1_contour_COM(2);
%             theta0 = atan2(dy_2, dx_2)+pi/2;
