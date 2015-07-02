function testSLN

    f = [0.01:.01:10];
    
    figure(2);
    
    hMain = plot(0,0);  
    
    line([min(f) max(f)], [0 0], 'Color', 'k', 'LineStyle', ':');
%     line([0 0],  [-100 100], 'Color', 'k', 'LineStyle', ':');
    function y = SLN(f, f_opt, w, s, r_max, r_bkg)
        
        A = (r_max - r_bkg)./(1-exp(-1/(s^2)))
        exp1s2 = (1-exp(-1/(s^2)));
        
%         logf_fopt_frac = zeros(size(f));
%         i = (f ~= 0);
%         logf_fopt_frac(~i) = 0;
%         logf_fopt_frac(i) = log(f(i)./f_opt)./(w + s.*log(f(i)./f_opt));
        logf_fopt_frac = log(f./f_opt)./(w + s.*log(f./f_opt));
        
        y =  A .* ( exp( -  ( logf_fopt_frac ).^2)  -exp(-1/s^2) ) + r_bkg;
    end
%     SLN = @(f, f_opt, w, s, r_max, r_bkg) ...
%         (r_max - r_bkg)./(1-exp(-1./s^2)) .* ( exp( -( (log(f./f_opt)./(w + s.*log(f./f_opt)) ).^2))  -exp(-1./s.^2) ) + r_bkg;
    

    function plotMain(f_opt, w, s, r_max, r_bkg)            
        y = SLN(f, f_opt, w, s, r_max, r_bkg);
        set(hMain, 'xdata', f, 'ydata', y);
    end
    
    func_handle = @(f_opt, w, s, r_max, r_bkg) plotMain(f_opt, w, s, r_max, r_bkg);
    args = { {'f_opt', [0:.01:10], 1}, {'w', [0:.01:10], 1}, {'s', [0:.01:10], 1}, {'r_max', [0:.01:10], 1}, {'r_bkg', [0:.01:10], 1} };
    manipulate(func_handle, args, 'FigId', 1);



end