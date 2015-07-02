function p = tmp_ttest(m, tail)
    x = [-.1:.01:1];%randn(100,1);
    [h,p] = ttest(x, m, .05, tail);
end