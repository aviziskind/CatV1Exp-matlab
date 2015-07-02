function m = distribMean(xs, ys)
    dx = xs(2) - xs(1);
    A = sum (ys)*dx;    % normalization;
    m = sum( (ys*dx) .* xs)/A;
end