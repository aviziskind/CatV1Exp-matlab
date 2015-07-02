function [A, mu, sigma, k, phi, extras] = estimate1DgaborParams(x,y, useFitFunction)

    x = x(:)';
    y = y(:)';

%     fitAlgorithm = 'nlinfit';
    fitAlgorithm = 'lsqcurvefit';
%     fitAlgorithm = 'fit';
    if nargin < 3
        useFitFunction = 1;
    end

    show = 0;
    
    kMethod = 'fft';
%     kMethod = 'peaks';
    
    gaborFunc = @(p, X) gabor1D(p(1), p(2), p(3), p(4), p(5), X);

    % MU
    mu_pow = 2;
    mu_est = sum( x .* abs(y).^(mu_pow) ) / sum(abs(y).^(mu_pow));

    % SIGMA
    sigma_est = distribStd(x, abs(y), mu_est);


    % find maxima & minima (for A, and for k (if use 'peaks' method);    
    [idx_extrema, ~] = findLocalMaxima(abs(y), 1);
    idx_extrema = [1, idx_extrema, length(y)];
    nExtrema = length(idx_extrema);        
    idx_closestToMax = idx_extrema( indmin(abs(x(idx_extrema) - mu_est)));
    
    
    if strcmp(kMethod, 'peaks')
        if nExtrema > 1

            idx_mu = indmin(abs(x - mu_est));
            distFromCenter = abs(idx_extrema - idx_mu);

            extrema_vals = y(idx_extrema);
            %             pair_SignExtremaProd = sign(extrema_vals(1:nExtrema-1) .* extrema_vals(2:nExtrema) );
            pair_MinDistFromCent = min( [distFromCenter(1:nExtrema-1); distFromCenter(2:nExtrema)], [], 1);
            pair_ValDiffs = abs( extrema_vals(1:nExtrema-1) - extrema_vals(2:nExtrema) );

            minDist_withSgnChange = min( pair_MinDistFromCent );
            idx_minDistPairs = find(pair_MinDistFromCent == minDist_withSgnChange);

            idx_bestPair = indmax ( pair_ValDiffs ( idx_minDistPairs) ) ;

            idx_pair = (idx_minDistPairs(idx_bestPair));
            idx1 = idx_extrema( idx_pair );
            idx2 = idx_extrema( idx_pair+1 );                        
            
            dx = abs(  x(idx1) - x(idx2)  );

        else
            idx1 = idx_extrema;            
            x1 = x(idx1);
            dx = min( x1-x(1), x(end)-x1 );
        end
        k_est = pi/dx;
    
    elseif strcmp(kMethod, 'fft');

        [frequencies, powers, phases] = powerSpectrum(x, y);
        idx_highestPower = indmax(powers);
        k_est = 2*pi*frequencies(idx_highestPower);
        
    end

            
    
    % estimate A
    y_peak = y(idx_closestToMax);
    x_peak = x(idx_closestToMax);
    A_est = abs(y_peak)/(exp ((x_peak-mu_est)^2 /(2*sigma_est^2)) );
    A_est = max(A_est, abs(y_peak));

    % estimate PHI
    phi_est = getPhaseOfCosine(x, y, k_est, mu_est);
    phi_est = mod(phi_est, 2*pi);

        j = idx_highestPower;
        phi_est_fft = mod( phases(j) - k_est*x(1), 2*pi); % for some reason this isn't always such a good estimate for phi
        y_fft = sqrt(powers(j))*cos(k_est*x  + phi_est_fft);


    % curve-fitting using all estimates
    p0 = [A_est, mu_est, sigma_est, k_est, phi_est];

    if ~useFitFunction
        fitAlgorithm = '';
    end
    
    switch fitAlgorithm
        case 'nlinfit'

        p_fit = nlinfit(x, y, gaborFunc, p0);

        A     = p_fit(1);
        mu    = p_fit(2);
        sigma = abs(p_fit(3));
        k     = abs(p_fit(4));
        phi   = mod(p_fit(5), 2*pi);
        
        %{
         % to plot:
            figure(56); clf; plot(x,y, '.'); hold on;
            plot(x, gaborFunc(p0, x), 'r');
            plot(x, gaborFunc(p_fit, x), 'g');
        %}
        
        case 'lsqcurvefit', 
            lower_bnds = [0, -inf, 0, 0, -inf];
            try
                p_fit = lsqcurvefit(gaborFunc,p0,x,y,lower_bnds,[], optimset('display', 'off'));
            catch
                p_fit = p0;
            end
            
            A     = p_fit(1);
            mu    = p_fit(2);
            sigma = p_fit(3);
            k     = p_fit(4);
            phi   = mod(p_fit(5), 2*pi);
            
        case 'fit'
    
            ft = fittype('gabor1D(A, mu, sigma, k, phi, x)');
            startPoint = [A_est, k_est, mu_est, phi_est, sigma_est];
            lower_bnds = [0, 0, -inf, -inf, 0];
            fo = fitoptions('Method', 'NonlinearLeastSquares', 'startpoint', startPoint, 'lower', lower_bnds);

            cfit = fit(x(:), y(:), ft, fo);

            A     = cfit.A;
            mu    = cfit.mu;
            sigma = cfit.sigma;
            k     = cfit.k;
            phi   = mod(cfit.phi, 2*pi);
            
        case ''
            % skip fit
            A = A_est;
            mu = mu_est;
            sigma = sigma_est;
            k = k_est;
            phi = phi_est;
             
            
        otherwise,
            eror('unknown algorithm');
            
    end
       
    
    extras.p_est = p0;
    extras.idx_closestToMax = idx_closestToMax;
%     extras.k_peaks_x = x([idx1, idx2]);
%     extras.k_peaks_y = y([idx1, idx2]);
    extras.y_fft = y_fft;
    

    if show
        
%         y_exp = A*exp( -(x-mu).^2./(2*sigma.^2));
        figure(101); clf; hold on;
        plot(x, y, '.');
        
        % show MU
%         set(h_line_mu, 'xdata', [mu mu], 'ydata', lims(y, .05), 'color', 'k', 'visible', tf2onOff(showMu) );
%         set(h_line_mu_est, 'xdata', [mu_est mu_est], 'ydata', lims(y, .05), 'color', 'r', 'visible', tf2onOff(showMu) );
%                         
%         % show sigma
%         y_exp_est = A*exp( -(x-mu).^2./(2*sigma_est.^2));        
%         set(h_sig, 'xdata', x, 'ydata', y_exp, 'color', 'k', 'visible', tf2onOff(showSigma));        
%         set(h_sig_est, 'xdata', x, 'ydata', y_exp_est, 'visible', tf2onOff(showSigma) );
%         
%         % show K
%         y_k = A * cos(k*(x-mu) + phi);
%         y_k_est = A*cos(k_est*(x-mu)+phi);
% %         y_k_est = extras.y_fft;
%         
%                         
%         set(h_k, 'xdata', x, 'ydata', y_k, 'visible', tf2onOff(showK) );
%         set(h_k_est, 'xdata', x, 'ydata', y_k_est, 'visible', tf2onOff(showK) );
%         set(h_peaks, 'xdata', extras.k_peaks_x, 'ydata', extras.k_peaks_y, 'visible', tf2onOff(showK) );
%         
%         % show A
%         y_A     = A*exp( -(x-mu).^2./(2*sigma.^2));                
%         y_A_est = A_est*exp( -(x-mu_est).^2./(2*sigma_est.^2));                
%         
%         set(h_A, 'xdata', x, 'ydata', y_A, 'color', 'k', 'visible', tf2onOff(showA) );
%         set(h_A_est, 'xdata', x, 'ydata', y_A_est, 'color', 'r', 'visible', tf2onOff(showA));
%         
%         
%         % show PHI 
%         y_phi = A*cos(k*(x-mu)+phi);
%         y_phi_est = A*cos(k*(x-mu)+phi_est);
%         set(h_phi, 'xdata', x, 'ydata', y_phi,  'visible', tf2onOff(showPhi) );
%         set(h_phi_est, 'xdata', x, 'ydata', y_phi_est,  'visible', tf2onOff(showPhi) );

        % show ALL estimates
        y_all_est =  A_est*exp( -(x-mu_est).^2./(2*sigma_est.^2)) .* cos(k_est*(x-mu_est) + phi_est);
        plot(x, y_all_est, 'r-');
        
        % show FIT
        if any(strcmp(fitAlgorithm, {'nlinfit', 'lsqcurvefit'}))       
            y_fit = gaborFunc(p_fit, x);
            plot(x, y_fit, 'g-');
        elseif strcmp(fitAlgorithm, 'fit')       
            plot(cfit, 'g-');            
        end
        
        
        p_labels = {'A', 'mu', 'sigma', 'k', 'phi'};                
        p_est  = [A_est, mu_est, sigma_est, k_est, phi_est];

        
        dispCell = [p_labels; num2cell(p_orig); num2cell(p_est); num2cell(p_fit)] %#ok<NOPRT,NASGU>
        
    
    end 
    
end