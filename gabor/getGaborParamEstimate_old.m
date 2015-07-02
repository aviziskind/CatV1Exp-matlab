function beta0 = estimateGaborParameters(xs, ys, zs)  % xs and ys are vectors, zs is a matrix;

    global Xc;
    if isempty(Xc), Xc = [0;0]; end;
    dbug = false;
    dbug_theta = true;
    smooth_w = 1.5;
    
    zs_orig = double(zs);
    zs = gaussSmooth(gaussSmooth(zs_orig, smooth_w, 1), smooth_w, 2);
    
%     zs = smooth2(zs, 3);

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
end



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
