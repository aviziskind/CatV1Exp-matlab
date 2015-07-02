function [theta, lambda] = findThetaLambdaWithStrongestOscillations(xs, ys, zs, X_center)
    dbug = 1;
    dbug_theta = false;
    
    thetas = 0:pi/50:pi;
    maxfreqs = zeros(1, length(thetas));
    amps     = zeros(1, length(thetas));
    rotateAndShift = @(X, theta, Xcenter)  rotationMatrix(theta) * X - repmat(Xcenter, 1,size(X,2));

    dx_orig = xs(2)-xs(1);
    dx = dx_orig/2;
    new_xs = [xs(1):dx:xs(end)] - mean(xs);
    Xs = [ new_xs; zeros(1,length(new_xs))];

    if dbug && dbug_theta
        thetaStems = stem3(Xs(1,:), Xs(2,:), zmax * ones(1, size(Xs,2)));
        figure(freq_fig_id); clf(freq_fig_id); hold on
        figure(main_fig_id);
        view(2);
    end

    [xs_grid, ys_grid] = meshgrid(xs,ys);

    for th_i = 1:length(thetas)
        [X_rotated] = rotateAndShift(Xs, thetas(th_i), -(X_center+Xc));

        xi = X_rotated(1,:);
        yi = X_rotated(2,:);
        zi = interp2(xs_grid,ys_grid,zs, xi, yi);

        zis = nonnans(zi);
        if ~isempty(zis)   % could be all nans if outside range
            [maxfreqs(th_i), amps(th_i), allFreqs, powers] = findStrongestFrequencies(dx, nonnans(zi), 1);
        else
            maxfreqs(th_i) = 0;
        end

        if dbug && dbug_theta
            figure(main_fig_id);
            set(thetaStems, 'xdata', xi, 'ydata', yi, 'zdata', zi)
            figure(freq_fig_id);
            plot3(allFreqs, ones(size(allFreqs))*thetas(th_i), powers, color(th_i));
            stem3(maxfreqs(th_i), thetas(th_i), amps(th_i), color(th_i));
            view(37, 52);
            xlim([0 .1]);
        end

    end

    % Remove funny high frequency blips.
    blip_threshold = 1/50;
    maxfreqs( amps < median(amps) * blip_threshold) = 0;

    ind_theta = indmax(amps .* maxfreqs);   % (ideally, use maxfreqs alone.)

    theta = thetas(ind_theta);
    lambda = 1/maxfreqs(ind_theta);
    
    %  kx <--> wt.   k = 2*pi/lamda.    w = 2*pi*f  ==> k = 2*pi*f
    %  lambda = 2*pi/k = 
    
%     kx k  
    
    
    
    
    
end