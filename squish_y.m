function Z2 = squish_y(Z, a)

    x = ones(5,1);
    y = reallocate(x, 3)
    
%     Z2 = Z;
%     ny = size(Z,2);
%     as = (ny:-1:0)*a/ny;
%     [.1] -> [.1 , .5,  0]
%     
%     
%     1 2 3 4 5
        


end

    
function y = reallocate(x, n)
    x = x(:)';
    N = length(x);
    y = zeros(1,n);
    parts = [0, (1:n)*N/n];

    for j = 1:n
        lower = parts(j);
        upper = parts(j+1);
        inds = [floor(lower) : floor(upper)]+1;        
        wgts = ones(size(inds));
        if length(wgts) == 1
            wgts = upper - lower;
        else
            wgts(1) = ceil(lower) - lower;
            wgts(end) = upper - floor(upper);
        end
        q = find(abs(wgts) < 1e-5);
        wgts(q) = []; inds(q) = [];
        y(j) = sum(x(inds) .* wgts);    
    end
    
end