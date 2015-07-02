% test_mode = true;
function test_selectWithLogAxes
    close all;

    dofig1 = 1;
    dofig2 = 1;

    xlog = 0;
    ylog = 1;

    shortTest = true;

    if shortTest
        x_scl = .5;
        y_scl = 10;
        x = [-x_scl, 0, 0, x_scl];
        y = [0, -y_scl, y_scl, 0];
        N = length(x);
    else

        N = 30;
        rnd_state = rand('state');
        rand('state', 0);
        x = rand(N,1);
        y = rand(N,1);

        rand('state', rnd_state); %#ok<*RAND>
        idx = ord(x);
        x = x(idx); y = y(idx);

    end
    b_x = 1;%exp(rand*2);
    b_y = 1;%exp(rand*2);

    if xlog    
%         x = 10.^((x-.5)*b_x);
        x = 10.^(x*b_x);
    end
    if ylog
%         y = 10.^((y-.5)*b_y);
        y = 10.^(y*b_y);
    end

%     th = linspace(0, 2*pi, 50);


    if dofig1
        figure(1); clf; plot(x,y, '.', 'userdata', 'main'); hold on;
        if xlog
            set(gca, 'xscale', 'log');
        end
        if ylog
            set(gca, 'yscale', 'log');
        end
%         xlim1 = xlim; ylim1 = ylim;
%         xlim1 = xlim1 + diff(xlim1)*[-1, 1]/10;
%         ylim1 = ylim1 + diff(ylim1)*[-1, 1]/10;
%         axis([xlim1 ylim1]);

        for i = 1:N
            h1(i) = plot(0,0, 'g');
        end
        set(h1, 'userdata', 'skipSelect');
        glob_i = 1;
        title(' ');        
        
        set(gcf, 'windowButtonDownFcn', @selectPointsInFigure)
        set(gcf, 'resizefcn', @redrawcircles);        

        % daspect([1 1 1]);


    end

    if dofig2
        figure(2); clf; hold on;
        x2 = x;
        y2 = y;

        if xlog
            x2 = log10(x);
        end
        if ylog
            y2 = log10(y);
        end
        plot(x2, y2, 'r.', 'userdata', 'main');

        xlim1 = xlim; ylim1 = ylim;
        xlim1 = xlim1 + diff(xlim1)*[-1, 1]/10;
        ylim1 = ylim1 + diff(ylim1)*[-1, 1]/10;
        axis([xlim1 ylim1]);        
        
        da2 = daspect;
        d = da2(2)/da2(1)
        ax2 = axis;
        axis(ax2);

        xlims2 = xlim;
        ylims2 = ylim;
        r_x2 = diff(xlims2)/15;
        r_y2 = diff(ylims2)/15;

        for i = 1:N
            h2(i) = plot(0, 0, 'g')    ;
        end
        set(h2, 'userdata', 'skipSelect');
        
        set(2, 'windowButtonDownFcn', @selectPointsInFigure)
        set(2, 'resizefcn', @redrawcircles);        
        
    end



    function redrawcircles(src, evnt);
%         return;
        da = daspect;
        3;
        % d = da(2) /da(1) 
%         d = log(da(2)/da(1))
        % d = 1;

        ax = get(src, 'children');
        xlog_here = strcmp(get(ax, 'xscale'), 'log');
        ylog_here = strcmp(get(ax, 'yscale'), 'log');
        ax_r = axaspect(ax);
        
        hnds = get(ax, 'children');        
        ud = get(hnds, 'userdata');
        
        hnd1 = hnds( strcmp(ud, 'main') );
        xx = get(hnd1, 'xdata');
        yy = get(hnd1, 'ydata');
        hnd_circ = hnds( strcmpi(ud, 'skipSelect') );
        n = length(xx);
        glob_i = glob_i + 1;
        set(get(ax,'Title'), 'string', num2str(glob_i));
        
        xlims1 = get(ax, 'xlim');
        ylims1 = get(ax, 'ylim');
        Q = 1/15;
        if xlog_here
            r_x1 = Q*diff(log10(xlims1)) /ax_r;
        else
            r_x1 = Q*diff(xlims1);
        end

        if ylog_here    
            r_y1 = Q*diff(log10(ylims1));
        else
            r_y1 = Q*diff(ylims1);
        end
        r_ratio = r_y1/r_x1;
        d_ratio = da(2)/da(1);
%         assert(r_ratio == d_ratio);
        
        theta = linspace(0, 2*pi, 50);
        for j = 1:n
            if xlog_here
                x_c = 10.^ ( log10( xx(j) ) + (r_x1*cos(theta)*ax_r ) );
            else
                x_c = xx(j) + r_x1*cos(theta)*ax_r;
            end
            if ylog_here
                y_c = 10.^ ( log10( yy(j) ) + (r_y1*sin(theta) ) );
            else
                y_c = yy(j) + r_y1*sin(theta);
            end
            set(hnd_circ(j), 'xdata', x_c, 'ydata', y_c, 'color', color_s(glob_i))    
        end
        

    end
end



