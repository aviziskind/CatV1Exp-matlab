function test1DGaborParamsEstimator

%     k = 1.2;
%     phi = deg2rad(20);
%     A = 2;
%     mu = 1;
%     sigma = 2;
    
    
    figure(11); clf;
    h1 = plot(0,0, 'b.'); hold on;
    h_ax = gca;
    h_sig = plot(0,0, 'g-');
    h_sig_est = plot(0,0, 'r-');
    
    h_k = plot(0,0, 'k');
    h_k_est = plot(0,0, 'r-');
    h_peaks = plot(0,0, 'rs');
    
    h_line_mu = line([0 0], [1 1]);
    h_line_mu_est = line([0 0], [1 1]);
    
    h_A = plot(0,0, 'g-');
    h_A_est = plot(0,0, 'r-');
    
    h_phi = plot(0,0, 'k');
    h_phi_est = plot(0,0, 'r-');
    
    tf2onOff = @(tf) iff(tf, 'on', 'off');
    
    h_all_est = plot(0,0, 'r:');
    h_all_fit = plot(0,0, 'g-');
    
    figure(12); clf; hold on;
    h_p_orig = plot(0,0, 'bs');
    h_p_est  = plot(0,0, 'ro');
    h_p_fit  = plot(0,0, 'k*');
    
    
    gabor1D = @(A, mu, sigma, k, phi, X) A*exp( -(X-mu).^2./(2*sigma.^2)) .* cos(abs(k)*(X-mu) + phi);    
    gaborFunc = @(p, X) gabor1D(p(1), p(2), p(3), p(4), p(5), X);
    
    function updatePlot(A, mu, sigma, k, phi_deg, showMu, showSigma, showA, showK, showPhi, showAll, showFit)
                
        randn('state', 0);
        addNoise = 1;
        
        phi = deg2rad(phi_deg);
        x = [-8:.2:8]+8;
        y = A*exp( -(x-mu).^2./(2*sigma.^2)) .* cos(k*(x-mu) + phi);
        if addNoise
            y = y+randn(size(y))*A/3;
            y = gaussSmooth(y, 2);
            y = y(:)';
        end       
        
        y_exp = A*exp( -(x-mu).^2./(2*sigma.^2));
        set(h1, 'xdata', x, 'ydata', y);

        [A, mu, sigma, k, phi, extras] = estimate1DgaborParams(x,y);
        p_est = extras.p_est;
        p_fit = [A, mu, sigma, k, phi];
        [A_est, mu_est, sigma_est, k_est, phi_est] = dealV(p_est);
        
        % show MU
        set(h_line_mu, 'xdata', [mu mu], 'ydata', lims(y, .05), 'color', 'k', 'visible', tf2onOff(showMu) );
        set(h_line_mu_est, 'xdata', [mu_est mu_est], 'ydata', lims(y, .05), 'color', 'r', 'visible', tf2onOff(showMu) );
                        
        % show sigma
        y_exp_est = A*exp( -(x-mu).^2./(2*sigma_est.^2));        
        set(h_sig, 'xdata', x, 'ydata', y_exp, 'color', 'k', 'visible', tf2onOff(showSigma));        
        set(h_sig_est, 'xdata', x, 'ydata', y_exp_est, 'visible', tf2onOff(showSigma) );
        
        % show K
        y_k = A * cos(k*(x-mu) + phi);
        y_k_est = A*cos(k_est*(x-mu)+phi);
%         y_k_est = extras.y_fft;
        
                        
        set(h_k, 'xdata', x, 'ydata', y_k, 'visible', tf2onOff(showK) );
        set(h_k_est, 'xdata', x, 'ydata', y_k_est, 'visible', tf2onOff(showK) );
        set(h_peaks, 'xdata', extras.k_peaks_x, 'ydata', extras.k_peaks_y, 'visible', tf2onOff(showK) );
        
        % show A
        y_A     = A*exp( -(x-mu).^2./(2*sigma.^2));                
        y_A_est = A_est*exp( -(x-mu_est).^2./(2*sigma_est.^2));                
        
        set(h_A, 'xdata', x, 'ydata', y_A, 'color', 'k', 'visible', tf2onOff(showA) );
        set(h_A_est, 'xdata', x, 'ydata', y_A_est, 'color', 'r', 'visible', tf2onOff(showA));
        
        
        % show PHI 
        y_phi = A*cos(k*(x-mu)+phi);
        y_phi_est = A*cos(k*(x-mu)+phi_est);
        set(h_phi, 'xdata', x, 'ydata', y_phi,  'visible', tf2onOff(showPhi) );
        set(h_phi_est, 'xdata', x, 'ydata', y_phi_est,  'visible', tf2onOff(showPhi) );

        % show ALL estimates
        y_all_est =  A_est*exp( -(x-mu_est).^2./(2*sigma_est.^2)) .* cos(k_est*(x-mu_est) + phi_est);
        set(h_all_est, 'xdata', x, 'ydata', y_all_est,  'visible', tf2onOff(showAll), 'color', 'r', 'linestyle', '-' );
        
        % show FIT
        if showFit                        
            y_fit = gaborFunc(p_fit, x);
            set(h_all_fit, 'xdata', x, 'ydata', y_fit, 'visible', 'on' )
        else
            set(h_all_fit, 'visible', 'off' )
        end
        
        p_labels = {'A', 'mu', 'sigma', 'k', 'phi'};        
        p_orig = [A, mu, sigma, k, phi];
        p_est  = [A_est, mu_est, sigma_est, k_est, phi_est];

        if showFit
            rel_diff_est = abs(p_est-p_orig) ./ p_orig; rel_diff_est(5) = circDist(p_est(5), p_orig(5), 2*pi)/2*pi;
            rel_diff_fit = abs(p_fit-p_orig) ./ p_orig; rel_diff_fit(5) = circDist(p_fit(5), p_orig(5), 2*pi)/2*pi;
            set(h_p_orig, 'xdata', 1:5, 'ydata', zeros(1,5))
            set(h_p_est, 'xdata', 1:5, 'ydata', rel_diff_est);                
            set(h_p_fit, 'xdata', 1:5, 'ydata', rel_diff_fit);                

            A = [p_labels; num2cell(p_orig); num2cell(p_est); num2cell(p_fit)]
        end
    end

    tfOp = [true, false];
    args = {{'A', [1:.1:10], 3}, {'mu', [-5:.1:15], 1}, {'sigma', [.1:.1:20], 2}, {'k', [.1:.1:10], 1.2}, {'phi', [0:10:360]}, ...
            {'showMU', tfOp, 0}, {'showSigma', tfOp, 0}, {'showA', tfOp, 0}, {'showK', tfOp, 0}, {'showPhi', tfOp, 0}, {'showAll', tfOp}, {'showFit', tfOp} }; 
    manipulate(@updatePlot, args, 'FigId', 4);
    
    

    % want to estimate 



end



%{



x = -10:.01:10;
sig = 1;

y = 1/sqrt(2*pi) * exp(- (x.^2)/(2) ) ;

y = exp(-(x.^2)) * cos(x);

y = 1/sqrt(2*pi*sig^2).* exp(- (x.^2)/(2*sig^2) ) .* cos(x);
Y = su
%}


%{          
  extrema_vals = y(idx_extrema);
            pair_SignExtremaProd = sign(extrema_vals(1:nExtrema-1) .* extrema_vals(2:nExtrema) );
            pair_MinDistFromCent = min( [distFromCenter(1:nExtrema-1); distFromCenter(2:nExtrema)], [], 1);
            pair_ValDiffs = abs( extrema_vals(1:nExtrema-1) - extrema_vals(2:nExtrema) );
                        
            
            idx_sgnChange = find(pair_SignExtremaProd == -1);
            minDist_withSgnChange = min( pair_MinDistFromCent( idx_sgnChange ) ); 
            idx_minDistPairs = find(pair_MinDistFromCent( idx_sgnChange ) == minDist_withSgnChange);
            
            idx_bestPair = indmax ( pair_ValDiffs ( idx_sgnChange(idx_minDistPairs)) ) ;
            
            idx_pair = idx_sgnChange(idx_minDistPairs(idx_bestPair));

            %}