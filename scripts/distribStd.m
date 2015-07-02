function [s, mu] = distribStd(xs, ys, mu)
    if length(xs) == length(ys)
        dx = xs(2) - xs(1);
    elseif length(xs) == 1
        dx = xs;
        xs = (0:length(ys)-1)*dx;
    end
    A = sum (ys)*dx;    % normalization;
    if nargin < 3
        mu = sum( (ys*dx) .* xs)/A;
    end
    s = sqrt(   sum((ys*dx) .* (xs - mu).^2 )/A    );
end
