function slides4



    doSlides = 6;

    
    if any(doSlides == 1) % little minature phase tuning curves for shuffling

        w = 3;
%         shifts = [20, -30];
%         shifts = [-25, 50];
        shifts = [40, 20];
        states = [10, 9];    
        for i = 1:length(shifts)
            figure(40+i); clf; hold on; box on;
            randn('state', states(i));
            n = 100;
            x = 1:n;
            y{i} = gaussSmooth(randn(1,n), n/7, [], 1);
            h(1,i) = plot(x,y{i}, 'b-', 'linewidth', w);
            set(gca, 'xtick', [], 'ytick', [])
            xlim([1 n]);
            ylim(lims(y{i}, .2));
        end
        3;
        for i = 1:length(shifts)
            figure(40+i);        
        
            y2{i} = circshift(y{i}, shifts(i));

            h(2,i) = plot(x,y2{i}, 'r-', 'linewidth', w);        
            set(h(1,i), 'linewidth', 2, 'linestyle', ':')
            3;
        end
        3;

    end



    
    
    

%%
    if any(doSlides == 2)
        th = linspace(0, 2*pi, 100);
        th_deg = rad2deg(th);
        y = sin(th+pi);
        idx = 1:4;
        Y = bsxfun(@plus, y(:)*.25, idx);

        figure(20); plot(th_deg, Y);
        axis ij
        set(gca, 'xtick', [0:90: 360], 'ytick', idx)
        xlim([0, 360]);
        xlabel('phase');

    end


    if any(doSlides == 4) % little minature phase shifts as x-axes
        ph_shift = [0, 90, 180, 270];
        ph = [0:360];
        for i = 1:4           
            figure(40+i); clf; 
            plot(ph, cos( deg2rad( ph+ph_shift(i) )), 'k-', 'linewidth', 5);
            set(gca, 'xtick', [], 'ytick', []);
            box off;
            xlim([-90, 360+90])
            ylim([-1.5, 1.5])
        end
        3;
        
    end

    % cartoon pair of phase tuning curves
    %%
    
    if any(doSlides == 5)
        %%
        x = [0:72];
        y0 = gaussian(x, 36, 8)';
        y1 = circshift(y0, -10);
        y2 = circshift(y0, 25);
        figure(77); box on;
        h = plot(x, y1, 'b', x, y2, 'r');
        set(h, 'linewidth', 3)
        xlim([x(1), x(end)]);
        set(gca, 'xtick', [], 'ytick', []);
        3;


    end
    
    
    if any(doSlides == 6)
        %%
        x = -1:.005:1;
        y = x;
        [xx, yy] = meshgrid(x,y);
        
%         phi = -pi/2;
        phi = 0;
        zz = gabor(1, 0, 0, .4, .38, 0, 4/2*pi, phi, 0, [xx(:), yy(:)]);
        Z = reshape(zz, size(xx));
        R = 0.9;
        th = 0.1;
        Z(Z > th) = 1;
        Z(Z < -th) = -1;
        Z(ibetween(Z, -th, th)) = 0;
        
        % outside the circle
        Z( (xx .^2 + yy.^2 > R.^2) ) = .5;

        figure(66); clf;
        imagesc(x, y, Z);
        colormap gray;        
%         axis([-1.2, 1.2, -1.2, 1.2])
        axis square tight;
        set(gca, 'xtick', [], 'ytick', []);
        caxis([-1, 1]);
        
        3;


    end    
    

end