function [p, ms] = getPrctilesMeanStd(x, prctiles_vals)
    p = prctile(x, prctiles_vals);
    ms = [mean(x), std(x)];

end