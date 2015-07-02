function tf = sig_gabor_pixels(gparam, xs, ys)

%     f = 0.01;

    if nargin < 3
        Gid = xs;
        [xs, ys] = getStimulusXY(Gid);
    end

    [xs_grid, ys_grid] = meshgrid(xs, ys);
    
    X = [xs_grid(:), ys_grid(:)];
    
    % p  = [A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, C]
    [A, mu_x, mu_y, sig_x, sig_y, theta, k, phi, C] = dealV(gparam); %#ok<NASGU>
    
    useMethod = 'sigmas';
    minFracOfMID = 0.01;
    
    switch useMethod
        
        case 'sigmas',
            gparam_noCosine = [1, mu_x, mu_y, sig_x, sig_y, theta, 0, 0, 0];    
            G_noCosine = reshape( gaborP(gparam_noCosine, X), size(xs_grid) );    
            nStd = 3;
            tf = zeros(size(xs_grid));
            while nnz(tf(:)) < numel(tf)*.01
                tf = G_noCosine > exp(-nStd); % 2 std deviations.
                nStd = nStd+1;
            end            
            
        case 'rotatedCov'
            mu = [mu_x; mu_y];
            C_orig = diag([sig_x, sig_y].^2);
            C_rot = rotatedCovMtx(C_orig, -theta);
            C_rot_inv = inv(C_rot);

            [twoStd_ellipse_x, twoStd_ellipse_y] = ellipsoidFromCov(mu, C_rot_inv, 3, 50);
            elps_xy = [twoStd_ellipse_x(:)'; twoStd_ellipse_y(:)'];

            tf = zeros(size(xs_grid));
            for i = 1:numel(tf)
                tf(i) = isPointInPolygon(X(i,:)', elps_xy);
            end

            
    
    end
    
%     z = xAx(Cinv, X', mu)
%     gaussianN
    
        
    
    show = 0;
    if show
        figure(655); clf;
        subplot(1,3,1); imagesc(xs, ys, G_orig); axis equal; hold on; plot(twoStd_ellipse_x, twoStd_ellipse_y, 'r.-')
        subplot(1,3,2); imagesc(xs, ys, tf_1); axis equal; hold on; plot(twoStd_ellipse_x, twoStd_ellipse_y, 'r.-')
        subplot(1,3,3); imagesc(xs, ys, tf);  axis equal;hold on; plot(twoStd_ellipse_x, twoStd_ellipse_y, 'r.-')
        hold on;
    end
        
    3;
    
end





function X = shiftAndRotate(X, Xc, theta)  
    X = bsxfun(@minus, X, Xc);
    X = rotationMatrix(-theta) * (X);
end

function X = rotateAndShift(X, theta, Xc)  
    X = rotationMatrix(-theta) * (X);
    X = bsxfun(@minus, X, Xc);

end





