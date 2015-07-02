function testMedFiltHi
    N = 1e6;
    samplingRateHz = 10000;
        
    t = [1:N]/sr;
    
    X = randn(1, N);
    
    hi_freq1 = 600;
% %     hi_freq2 = 800;
    pow2db = @(p) (10.*log10(p));
    
%     filterings = {{1, hi_freq1, 'high', 'butter'}, {50, hi_freq1, 'high', 'fir1'}};
    
    filterings = {{1, hi_freq1, 'high', 'med'}, {50, hi_freq1, 'high', 'med'}};

    ns = [11:2:25];    
    freqs = (2./ns)*(sr/2);
    ns = (2./freqs)*(sr/2);
    

    for i = 1:length(filterings)
        
        filtering_i = filterings{i};
        offset = floor(filtering_i{2}/2);
        samp_filt{i} = filterSamples(X, filtering_i, samplingRateHz);
        figure(11);
        plot(samp_filt{i}(offset+1:end), ['.-' color_s(i)]);        
     
        [Pow{i}, freq{i}] = pwelch(samp_filt{i}, [], [], [], samplingRateHz);
        figure(12);
        plot(freq{i}, pow2db(Pow{i}), 'color', color_s(i)); 
        
    end    

    
    
        [b,a] = filterCoefficients(1, hi_freq1/(sr/2), 'high', 'butter');
        X_f{i} = filter(b, a, X);

        [b,a] = filterCoefficients(50, hi_freq1/(sr/2), 'high', 'fir1');
        X_f2 = filter(b, a, X);

        Wn = hi_freq1/(sr/2);
        X_f3 = medfilt1_hi(X, Wn);
    
    if odd(floor(n))
        n = floor(n);
    else
        n = ceil(n);
    end
    if ~odd(n)
        n = n+1;
    
    end
    
%     figure(1); clf;
%     pwelch(X, [], [], [], sr);
        
    %%
    nSm = 100;
    figure(1); clf; hold on;
%     g_conv = gaussian(-10:10, 0, 1); g_conv= g_conv/sum(g_conv);
    
    [Pow_raw, freq_raw] = pwelch(X, [], [], [], sr);      Pow_raw = fastmedfilt1d(Pow_raw, nSm);
    plot(freq_raw, pow2db(Pow_raw));     
    
    [Pow_f1, freq_f1] = pwelch(X_f1, [], [], [], sr);     Pow_f1 = fastmedfilt1d(Pow_f1, nSm);
    plot(freq_f1, pow2db(Pow_f1), 'g'); 

    [Pow_f2, freq_f2] = pwelch(X_f2, [], [], [], sr);     Pow_f2 = fastmedfilt1d(Pow_f2, nSm);
    plot(freq_f2, pow2db(Pow_f2), 'r'); 

    [Pow_f3, freq_f3] = pwelch(X_f3, [], [], [], sr);     Pow_f3 = fastmedfilt1d(Pow_f3, nSm);
    plot(freq_f3, pow2db(Pow_f3), 'm'); 
    
    
    drawVerticalLine(hi_freq1);
    
%     xlim([0 1000]); 
    xlabel('Hz');
%     set(gca, 'xscale', 'log')
    legend({'Raw', 'Butter', 'Fir1', 'med'})
    3;
%     linkaxes(h_ax([2, 1]), 'y');
    3;


end


    
function data_filtered = filterSamples(samp_raw, filter_opt, samplingRateHz)

    samp_raw = bsxfun(@minus, samp_raw, samp_raw(:,1));

    % Filter settings :
    normFreq = filter_opt{1}/(samplingRateHz/2);                
    filterOrder = filter_opt{2};    
    funcName = filter_opt{3};
    ftype = 'high';

    if ~strcmp(funcName, 'med')
    
        [b,a] = filterCoefficients(filterOrder, normFreq, ftype, funcName);    
        filter_dim = 2;    
    
        data_filtered = filter(b, a, samp_raw, [], filter_dim); 
    else        
        data_filtered = medfilt1_hi(samp_raw, normFreq);
    end        
    
end