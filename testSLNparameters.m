function testSLNparameters

    r_max = 1;

    r_bkg = 0;
    f_opt = 2;
    log_f_opt = log(f_opt);
    
    ws = [.05:.05:1];
    ss = [-1:.1:1];
   
    Nw = length(ws);
    Ns = length(ss);
    [LOs, HIs, Cnt1, Cnt2] = deal(zeros(Nw, Ns));
    
    f = log(.1:.01:4);
    tic;
    for wi = 1:Nw
        for si = 1:Ns
                         
            [lo, hi, Cnt1(wi, si), Cnt2(wi, si)] = getHalfHiLo(@(x) skewLogNormal_exp(x, f_opt, r_max, r_bkg, ws(wi), ss(si)), f_opt);
            LOs(wi, si) = lo;
            HIs(wi, si) = hi;
            
        end
    end
    toc;
    
%     figure(21); surf(ws, ss, LOs'); xlabel('w'); ylabel('s'); title('lo');
%     figure(22); surf(ws, ss, HIs'); xlabel('w'); ylabel('s'); title('hi');
    R = (HIs-f_opt)./(f_opt-LOs);
    figure(23); surf(ws, ss, R'); xlabel('w'); ylabel('s'); title('ratio');
    S = HIs+LOs-f_opt*2;
    D = HIs-LOs;
    
    W_est = D./(2);
    S_est = S./(1.25.*W_est);

    
    
    figure(24); surf(ws, ss, S'); xlabel('w'); ylabel('s'); title('sum');
    figure(25); surf(ws, ss, D'); xlabel('w'); ylabel('s'); title('diff');
    
    figure(26); surf(ws, ss, W_est'); xlabel('w'); ylabel('s'); title('W (est)'); view([-pi/2, 0])
    figure(27); surf(ws, ss, S_est'); xlabel('w'); ylabel('s'); title('S (est)'); view([92, 0])
    3;
    
    
end


function [lo, hi, count1, count2] = getHalfHiLo(f, x_peak)

    offset = .1;
    hi = x_peak;
    count1 = 0;
    count2 = 0;

    % verify that is peak
%     [x_peak_search] = fminsearch(@(x) -f(x), x_peak);
%     assert( abs(x_peak - x_peak_search) < 1e-3 );

    while (x_peak >= hi)
        [hi, ~, exitflag] = fzero(@(x) f(x)-.5, x_peak + offset);
        if exitflag ~= 1
            hi = x_peak;
        end        
        offset = offset * -1.5;
        count1 = count1+1;
    end
    
    offset = .1;
    lo = x_peak;
    while (x_peak <= lo)
        [lo, ~, exitflag] = fzero(@(x) f(x)-.5, x_peak - offset);
        if exitflag ~= 1
            lo = x_peak;
        end
        offset = offset * (-4);
        count2 = count2+1;
    end
    
    if (count1 > 50) || (count2 > 50)
        3;
    end
%     if lo > hi
%         [lo, hi] = deal(hi, lo);
%     end
end