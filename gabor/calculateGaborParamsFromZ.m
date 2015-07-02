function p = calculateGaborParamsFromZ(xs, ys, zs)

    if nargin == 1
        zs = xs;
        [nx, ny] = size(zs);
        xs = 1:nx;
        ys = 1:ny;
    end

	dbug = true;

    %%%%    GET XS AND YS
    [xs, ys, xs_grid, ys_grid] = parseXandY(xs, ys);

%     gaborFunc_XY = @(P, XY) reshape(gaborP(P, [XY{1}(:), XY{2}(:)]), size(XY{1}));
%     zs_estm = gaborFunc_XY(p_est, XY_C); 
        
    XY_all = [xs_grid(:), ys_grid(:)];
    
    if dbug
        smooth_w = 1.5;
        zs_sm = gaussSmooth(gaussSmooth(zs, smooth_w, 1), smooth_w, 2);
        
        figure(100); clf(100);
%         h1 = mesh(xs, ys, zs_sm, 'EdgeColor', 'b'); hold on;
        h1 = mesh(xs, ys, zs_sm'); hold on;
        xlims = [xs(1) xs(end)];
        ylims = [ys(1) ys(end)];
        xlabel('X');
        ylabel('Y');
    end
    
    %%%%    ESTIMATE GABOR FUNCTION.
%     gaborFunc = @(P, X) gabor(P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), X);
    all_p0 =  estimateGaborParameters(xs, ys, zs);        
    p_fit = fitGaborParams(xs, ys, zs, all_p0);
    
    if dbug
        z_est = reshape(gaborFunc(p0, XY)', size(xs_grid'));
        
%         h2 = mesh(xs, ys, z_est); % @(X) gaborFunc(p0, X), xlims, ylims, 'spacing', [diff(xlims)/100, diff(ylims)/100]);  
%         set(h2, 'EdgeColor', 'g');
    end
    
    %%%%    CALCULATE ACTUAL GABOR FUNCTION.
%     minimumValue = .05 * max(abs(zs(:)));
%     [y_sel_ind, x_sel_ind, z_sel] = find( (abs(zs) >= minimumValue) .* zs);

    
    XX = [xs_grid(:), ys_grid(:)];
    YY = zs_sm';
    YY = YY(:);
    p = nlinfit(XX, YY, gaborFunc, p0);
%     gaborFunc(p0, XX)
    
    if dbug        
        h3 = fmesh(@(X) gaborFunc(p, X), xlims, ylims);  set(h3, 'EdgeColor', 'm');
        3;
    end
    3;

end

function [xs, ys, xs_grid, ys_grid] = parseXandY(xs, ys)
    if (isvector(xs) && isvector(ys))
        [xs_grid, ys_grid] = meshgrid(xs, ys);
    else 
        xs_grid = xs;
        ys_grid = ys;
        xs = xs_grid(1,:);
        ys = xs_grid(:,1);
    end
end