function testGabor
    
%         A = - (rand*10 + .1);
%         Xc = [0;0];
%         mu_x = rand*5 - 2.5;
%         mu_y = rand*5 - 2.5;
%         sig_x =  .1 + rand*3;
%         sig_y =  .1 + rand*3;
%         k = 0;%rand*pi;
%         phi = rand*2*pi;
%         theta = rand*pi;
%         const = 0;%2*rand - 1;

    testGradient = true;

    xlims = [-4 4]; Nx = 100;
    ylims = [-4 4]; Ny = 102;
    xs = linspace(xlims(1), xlims(2), Nx);
    ys = linspace(ylims(1), ylims(2), Ny);
    [xs_grid ys_grid] = meshgrid(xs, ys);
    sizeX = size(xs_grid);
    zs = zeros( sizeX );    
    
    A0 = 5;
    mu_x0 = -.5;
    mu_y0 = 2;
    sig_x0 = 1.2;
    sig_y0 = .5;
    k0 = 5;%rand*pi;
    phi0 = pi/6;
    theta0 = pi/5;
    const0 = 0;%2*rand - 1;

%     gparams = {A0, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const};
    pLabels = {'A', 'mu_x', 'mu_y', 'sig_x', 'sig_y', '\theta', 'k', '\phi', 'c'};
    nParams = 9;

    if ~testGradient
        figure(15); clf;
        h_im = imagesc(xs, ys, zs);        

        axis equal tight xy;
        colormap('gray');    
        colorbar;

        hold on;
        h_cent  = plot(0,0, 'gs', 'markersize', 2, 'color', 'g');
        h_theta = plot(0,0, 'ro');
        
    else
        figure(20); clf; set(20, 'name', 'gradient');
        for sub_i = 1:nParams
            subplot(2,5, sub_i);            
            h_im(sub_i) = imagesc(xs, ys, zs);      %#ok<AGROW>
%             h_ax(sub_i) = gca;                      %#ok<AGROW,NASGU>

            axis equal tight xy;
            colormap('gray');    
            title(pLabels{sub_i});
            set(gca, 'xtick', [], 'ytick', []);
        end
        
        figure(21); clf; set(21, 'name', 'finite difference approximation');
        for sub_i = 1:nParams
            subplot(2,5, sub_i);            
            h_im_est(sub_i) = imagesc(xs, ys, zs);      %#ok<AGROW>
%             h_ax(sub_i) = gca;                      %#ok<AGROW,NASGU>
            h_ylab_est(sub_i) = ylabel('');
            axis equal tight xy;
            colormap('gray');    
            title(pLabels{sub_i});
            set(gca, 'xtick', [], 'ytick', []);
        end
    end
    
%     zs = reshape( gabor(gparams{:}, [xs_grid(:), ys_grid(:)]) , size(xs_grid));            
    
%     Xc = [0;0];
%     Xc + [mu_x; mu_y];
%     dXarrow = [-sin(theta); cos(theta)];
%     plot(Xmid(1), Xmid(2), 'ro');
%     quiver(Xmid(1), Xmid(2), dXarrow(1), dXarrow(2));
%     hold off
        
    function updateGabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const)
                
        Z = gabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const, [xs_grid(:), ys_grid(:)]);
        Z = reshape(Z, size(xs_grid));
        set(h_im, 'cdata', Z);        
%         L = max(abs(Z(:)))+.1;
%         set(h_ax, 'clim', [-L, L]);        
        
        set(h_cent, 'xdata', [mu_x], 'ydata', mu_y, 'visible', 'on');
        set(h_theta, 'visible', 'off');
        
    end


    function updateDGabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const)
                
        dZ = dgabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const, [xs_grid(:), ys_grid(:)]);
        pSiz = [size(xs_grid), nParams];
        dZ = reshape(dZ, pSiz);

        dZ_est = dgabor(A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, const, [xs_grid(:), ys_grid(:)], 1);
        dZ_est = reshape(dZ_est, pSiz);
        
        for p_i = 1:nParams
            set(h_im(p_i), 'cdata', dZ(:,:,p_i));                
        end
        for p_i = 1:nParams
            set(h_im_est(p_i), 'cdata', dZ_est(:,:,p_i));                
            
            D = dZ(:,:,p_i) - dZ_est(:,:,p_i);
            D_abs = mean(abs(D(:)));
            if D_abs < 1e-10;
                D_str = 'ok';
                D_col = 'k';
            else
                D_str = sprintf('D = %.3g', D_abs);
                D_col = 'r';
            end
            set(h_ylab_est, 'string', D_str, 'color', D_col);
            
        end        
        
        
    end


    args = {{'A', [-5:.1:5], A0}, {'mu_x', [-3:.1:3], mu_x0}, {'mu_y', [-3:.1:3], mu_y0}, ...
            {'sig_x', [.1:.1:2], sig_x0}, {'sig_y', [.1:.1:2], sig_y0}, {'theta', [0:.1:pi], theta0}, ...
            {'k', [.1:.1:10], k0}, {'phi', [0:.1:2*pi], phi0}, {'const', [-5:.1:5], 0} };

    if ~testGradient
        manipulate(@updateGabor, args, 'FigId', 3);
    else
        manipulate(@updateDGabor, args, 'FigId', 3);
    end
                
        
    
end


function adjustZlimsFor(M)
    mn = min(M(:));
    mx = max(M(:));
    d = mx-mn;
    f = 1e5;
    if f*d > 1, f = 1/d; end
    zlim([mn- f*d, mx + f*d])
end


%{
        if testGradient
            plotGradients = false;
            
            gaborGradientEstimate = zeros(length(ys), length(xs), length(p));
    %         gaborGradients = zeros(length(xs), length(ys), length(p));
            % estimate gradient empirically
            delta = 1e-5;
            tic;
            gaborGradients = gabor(p, {xs_grid, ys_grid}, 'gradient');
            toc;
            tic;
            for p_j = 1:length(p);
                pplus = p;   pplus(p_j)  = pplus(p_j)  + delta;
                pminus = p;  pminus(p_j) = pminus(p_j) - delta;
                gaborGradientEstimate(:,:,p_j) =  (gabor(pplus, {xs_grid,ys_grid}) - gabor(pminus, {xs_grid,ys_grid}) ) / (2*delta);

                discrep = maxElement(gaborGradientEstimate(:,:,p_j) - gaborGradients(:,:,p_j)) ;
                if discrep > delta
                    disp(['Discrepancy for parameter ' num2str(pi) ' of order ' num2str(discrep)])
                end
                if plotGradients
                    figure(p_j); clf(p_j);
                    mesh(xs_grid, ys_grid, gaborGradientEstimate(:,:,p_j), 'EdgeColor', 'b', 'FaceAlpha', 0); hold on
                    mesh(xs_grid, ys_grid, gaborGradients(:,:,p_j),        'EdgeColor', 'g', 'FaceAlpha', 0); hold off
                    title(pLabels{p_j});

                    if p_j == 9
                        adjustZlimsFor(gaborGradientEstimate(:,:,p_j))
                    end
                end
                    
            end
            toc;
    
        end    
%}