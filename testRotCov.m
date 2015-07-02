function testRotCov


    th = linspace(0,2*pi, 41); th = th(1:end-1);

    figure(81); clf;
    h = plot(0,0, '.-');
    hold on;
    axis equal;
    h_ax = gca;

    randn('state', 0);
    n = 50000;
    x0 = randn(1,n);
    x0 = x0-mean(x0); x0 = x0/std(x0);
    y0 = randn(1,n);
    y0 = y0-mean(y0); y0 = y0/std(y0);
    
    he1 = plot(0,0, 'g-');
    he2 = plot(0,0, 'r-');
    
    XY = mvnrnd([0;0], [sig_x, .3; .3, sig_y], 5000)';
    function updatecov(sig_x, sig_y, theta_deg)
%         n = 10000;
%         XY = mvnrnd([0;0], [sig_x, .3; .3, sig_y], 5000)';
        
%         x = x0*sig_x;
%         y = y0*sig_y;
        x = XY(1,:);
        y = XY(2,:);
        
%         x = cos(th)*(sig_x^2);
%         y = sin(th)*(sig_y^2);

        XY = [x(:)'; y(:)'];
        
        theta = deg2rad(theta_deg);
        R = rotationMatrix(theta);
        XY_rot = R*XY;
        x_rot = XY_rot(1,:);
        y_rot = XY_rot(2,:);        
        
        s = sin(theta);
        c = cos(theta);

        Cov1 = cov(XY');
        Cov2 = cov(XY_rot');
        
%         Cov1_pred = [sig_x^2, 0; 0, sig_y^2]
%         Cov2_pred = [sig_x^2*c^2 + sig_y^2*s^2,   (sig_x^2-sig_y^2)*s*c;
%                     (sig_x^2-sig_y^2)*s*c,         sig_y^2*c^2 + sig_x^2*s^2]
        Cov2_pred = rotatedCovMtx(Cov1, theta);
        Cov2-Cov2_pred
        3;
        
        [xe1, ye1] = ellipsoidFromCov([0;0], Cov2, 2, 100);
        [xe2, ye2] = ellipsoidFromCov([0;0], Cov2_pred, 2, 100);
                
        set(h, 'xdata', x_rot, 'ydata', y_rot, 'linestyle', 'none');
        
        set(he1, 'xdata', xe1, 'ydata', ye1);
        set(he2, 'xdata', xe2, 'ydata', ye2);
        3;
        
    end

    args = { {'sig_x', [1:.1:5]}, {'sig_y', [1:.1:5]}, {'theta', [0:5:360]} };

    manipulate(@updatecov, args, 'FigId', 82);

end