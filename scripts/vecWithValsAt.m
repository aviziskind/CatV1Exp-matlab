function V = vecWithValsAt(n, idxs, vals)
    V = zeros(n, 1);
    V(idxs) = vals;
end