function [p_fit, rsqr, Z_fit, idx_bestFit, order_bestFits, all_p0] = fitGaborParams(xs, ys, zs, all_p0, nMaxFits)

    fitAlgorithm = 'lsqcurvefit';
    if nargin < 5
        nMaxFits = 200;
    end
%     gaborFunc = @(P, X) gabor(P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), X);

    restrictToSigPoints = 0;
    if restrictToSigPoints
        zmax = max(abs(zs(:)));
        [y_sel_ind, x_sel_ind, z_sel] = find( (abs(zs) >= .1*zmax) .* zs);
        XX = [xs(x_sel_ind)', ys(y_sel_ind)'];
        YY = z_sel;
    else
        [x_grid, y_grid] = meshgrid(xs, ys);
        XX = [x_grid(:), y_grid(:)];
        YY = zs(:);    
    end    
    
    nEsts_orig = size(all_p0, 1);
    if nEsts_orig > nMaxFits
        resnorms_est = zeros(1, nEsts_orig);
        zs_col = zs(:);
        for i = 1:nEsts_orig
            zs_est = gaborP(all_p0(i,:), XX);
            resnorms_est(i) = sum( (zs_col-zs_est).^2);            
        end
        [~, idx_resnorm] = sort(resnorms_est, 'ascend');
                
        all_p0 = all_p0( idx_resnorm (1:nMaxFits), :);        
    end
    
    [nEsts, nVar] = size(all_p0);
    assert(nVar == 8 || nVar == 9);
    resnorms = zeros(1, nEsts);
    p_fits = cell(1,nEsts);
    
    
    
    
    if strcmp(fitAlgorithm, 'lsqcurvefit')
        
        Dz = lims(zs, 1); % increase by 100%
        A_bnds = [0, diff(Dz)];
        Dx = xs(end)-xs(1); dx = diff(xs(1:2));
        Dy = ys(end)-ys(1); dy = diff(ys(1:2));
        mu_x_bnds = lims(xs, 5); % increase range by 500%
        mu_y_bnds = lims(ys, 5); 
        sig_x_bnds = [0, inf];
        sig_y_bnds = [0, inf];
        theta_bnds = [-inf, inf];

        spPeriod_pix_min = (min(dx,dy)*2) * (.95);  % 1 cycle every 2 pixels
        spPeriod_pix_max = (max(Dx,Dy)*2)*2;        % 1/2 of a cycle over the entire image
        k_bnds = [2*pi/spPeriod_pix_max, 2*pi/spPeriod_pix_min];
        phi_bnds = [-inf, inf];
        const_bnds = [-inf, inf];

        lower_bnds = [A_bnds(1), mu_x_bnds(1), mu_y_bnds(1), sig_x_bnds(1), sig_y_bnds(1), theta_bnds(1), k_bnds(1), phi_bnds(1), const_bnds(1)];
        upper_bnds = [A_bnds(2), mu_x_bnds(2), mu_y_bnds(2), sig_x_bnds(2), sig_y_bnds(2), theta_bnds(2), k_bnds(2), phi_bnds(2), const_bnds(2)];
    end
        
    doProgressBar = nEsts > 5;
    if doProgressBar
        fprintf('Trying %d (out of %d) fits ... ', nEsts, nEsts_orig);
        progressBar('init', nEsts, 10)
    end
        
    for est_i = 1:nEsts
        p0 = all_p0(est_i,:);
                
        switch fitAlgorithm
            case 'nlinfit',
                [p_fits{est_i}, resid] = nlinfit(XX, YY, @gaborP, p0, statset('maxiter', 200));
                resnorms(est_i) = sum(resid.^2);

            case 'lsqcurvefit',    
                [p_fits{est_i}, resnorms(est_i)] = lsqcurvefit(@gaborP, p0, XX, YY, lower_bnds, upper_bnds, optimset('Jacobian', 'on', 'display', 'off'));
        end
        if doProgressBar
            progressBar(est_i);
        end
    end
    progressBar('done'); fprintf('\n');
            
    idx_bestFit_tmp = indmin(resnorms);
    p_fit = p_fits{idx_bestFit_tmp};
    p_fit = standardizeGaborP( p_fit );
    
    Z_fit = reshape(gaborP(p_fit, XX), size(zs));
    rsqr = (corr(zs(:), Z_fit(:))).^2;

    
    if nEsts > 1
        normDists = zeros(1,nEsts);
%         std_p0 = std(all_p0, [], 1);
        tolfun = 1e-6;
        L = diff(lims(resnorms));
        if L > eps
            L = L*eps;
        end
    
        for i = 1:nEsts            
            z_est = reshape(gaborP(all_p0(i,:), XX), size(zs));
%             normDists(i) = norm( (p_fit-all_p0(i,:)) ./std_p0 );
            normDists(i) = norm( Z_fit-z_est );
            
            diffs = abs(resnorms - resnorms(i));
            idx_tinyDiff = find(between(diffs, 0, tolfun*2));
            if ~isempty(idx_tinyDiff)
                resnorms(idx_tinyDiff) = resnorms(i);
            end
        end        
        
        
        [A, order_bestFits] = sortrows([resnorms(:), normDists(:)]);  % if multiple estimates gave the same final fit, also sort by how close they were to the original.    
%         C = num2cell(A);
        idx_bestFit = order_bestFits(1);

        all_p0 = all_p0(order_bestFits,:);
    else
        idx_bestFit = 1;
        order_bestFits = 1;
    end
        
    
    
        
end