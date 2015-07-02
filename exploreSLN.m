function exploreSLN

    f = 1e-10:.0001:4;
    y = zeros(size(f));
    figure(11); clf;
    h = plot(f, y); hold on;
    h_ax = gca;
    
    h_line_min = plot(0,0, 'k');
    h_line_max = plot(0,0, 'r');
    
    f_opt0 = 1;
    r_bkg0 = 0;
    r_max0 = 1;
    w0 = 1;
    s0 = 0;
    r_max = 1;

    %     SLN = @(f_opt, r_bkg, w, s, f) (r_max - r_bkg)./(1-exp(-1/s.^2)) * ( exp( - ( log(f./f_opt)./(w + s.*log(f./f_opt) ) ).^2 ) - exp(-1./s.^2) ) + r_bkg;

    function updatePlot(f_opt, r_max, r_bkg, w, s, xscale)
        
       f = .0001:.001:10;
       switch xscale
           case 'linear'
                f_input = f;
                y = skewLogNormal(f, f_opt, r_max, r_bkg, w, s);           
                xscale_input = 'linear';
           case 'log',
                f_input = f;
                y = skewLogNormal_exp(f, f_opt, r_max, r_bkg, w, s);
                xscale_input = 'log';
           case 'linear axis with log values'                
                f_input = log10(f);
                y = skewLogNormal_exp( (f), (f_opt), r_max, r_bkg, w, s);
                xscale_input = 'linear';
       end
       
       
       set(h, 'xdata', f_input, 'ydata', y, 'marker', '.');        
       set(h_line_min, 'xdata', f_input([1, end]), 'ydata', [r_bkg, r_bkg]);
       set(h_line_max, 'xdata', f_input([1, end]), 'ydata', [r_max, r_max]);

%        warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
       set(h_ax, 'xscale', xscale_input);
       drawnow;
       xlim(h_ax, [f_input(1), f_input(end)]);
%        warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');
       
    end

    args = {{'fopt', [.001:.01:3.5], 1}, {'rmax', [.1:.01:2], 1}, {'rbkg', [0:.01:2]}, {'w', [-1:.003:1], .98}, {'s', [-5:.05:5], .2}, {'xscale', {'linear', 'log', 'linear axis with log values'}}};
        
    
    manipulate(@updatePlot, args, 'FigId', 12);
    
%     SLNfunc = @(b, f) SLN(b(1), b(2), b(3), b(4), f);



end