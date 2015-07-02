function barPlotTwoSamples(x1, x2, label1, label2, measure_label)
    if nargin < 5
        measure_label = {};
    else
        measure_label = {measure_label};
    end
    xx = [1, 2]; 
    yy = [mean(x1), mean(x2)]; 
    ee = [stderr(x1), stderr(x2)];
    bar(xx, yy); 
    hold on;
    errorbar(xx, yy, ee, 'r.');
    set(gca, 'xtick', [1 2], 'xticklabel', {label1, label2});
    xlim([.5, 2.5]);
    [~, p_t] = ttest(x1, x2);
    p_w = ranksum(x1, x2);
%     title(sprintf('p_t = %.2g. p_W = %.2g', p_t, p_w))    
    title({measure_label{:}, ...
        sprintf('%s : mean : %.2g. median: %.2g ', label1, mean(x1), median(x1)), ...
        sprintf('%s : mean : %.2g. median: %.2g ', label2, mean(x2), median(x2)), ...
        sprintf('p_t = %.2g. p_W = %.2g', p_t, p_w)});

    
end