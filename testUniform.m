function [tf, pval] = testUniform(x, y)
    alpha = 0.01;
    [tf, pval, ksstat, cv] = kstest(y, [0 0; 1 1], alpha)

end
