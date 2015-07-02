function phi = getPhaseOfCosine(xs, ys, k, x0, flag)
%%
    % given data y = cos(k*(x-x0) + phi)
    % this function returns the value of phi, given y, x, k, and x0

    dbug = exist('flag', 'var') && ~isempty(flag);
    
    dx = xs(end)-xs(round(end/2));
    dys = 1 - (.2/(dx^2)) * (xs).^2;
    
    [mx, ind_peak] = max(abs(ys).*dys);         % find index of highest/lowest point.
    [mn, ind_center] = min(abs(x0 - xs));  % find index of closest place to center
   
%     distCentToPeak = abs(xs(ind_peak) - xs(ind_center));
    % the cosine reaches its max when (k(x-x0) + phi) = 0.
    if ys(ind_peak) > 0   % peak in front / trough behind
        phi = mod(  -k*(xs(ind_peak) - xs(ind_center))       , 2*pi);
    else                  % trough in front / peak behind
        phi = mod(   pi - k*(xs(ind_peak) - xs(ind_center))  , 2*pi) ;
    end

    if dbug
        %%
        figure(534); clf
        plot(xs,ys, 'b.');
        hold on;
        plot(xs(ind_peak),   ys(ind_peak),   'ro', 'MarkerFaceColor', 'r');  
        plot(xs(ind_center), ys(ind_center), 'ks');
        plot(xs, cos(k*(xs-x0)+phi) * max(abs(ys)) );  
        hold off
    
    end
    
end


% test script:
% k = randn;
% phi = randn;
% x0 = randn;
% x = -4*pi:.02:4*pi; y =  cos(k*(x-x0)+phi) .* gaussian(x, 0, 10); plot(x,y, '.');
% findPhi(x, y, k, x0)