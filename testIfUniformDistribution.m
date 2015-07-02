function testIfUniformDistribution(px)

    N = 1000;
    n = length(px);
    px_rnd = round( N * px/sum(px) );
    x = linspace(0, 1, n+1);
    x = x(1:n);
    
    x_new = zeros(1, sum(px_rnd));    
    inds = [1 cumsum(px_rnd)];
    for i = 1:length(px_rnd)
        x_new([inds(i):inds(i+1)]) = x(i);
    end
    x_new;
    figure(1); hist(x_new, n);

    [h, p] = kstest_prob( x_new, [0 0; max(x_new) 1] )
    
end