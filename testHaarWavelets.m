function testHaarWavelets

    function y = haar(x, a, b)
        x = (x-b)/a;
        s = 1/sqrt( abs(a) );
        
        if (x < 0) || (x > 2)
            y = 0;
        elseif x < 1
            y = s;
        elseif x <= 2
            y = -s;
        else
            3;
        end
            
    end
       
    j = 1:4;
    k = 1:4;
    [j_grid, k_grid] = meshgrid(j,k);
    as = arrayfun(@(J,K) 2.^(-J), j_grid, k_grid);  
    bs = arrayfun(@(J,K) 2.^(-J)*K, j_grid, k_grid);

    i = 1;
    for ji = 1:length(j)
        for ki = 1:length(k)
            subplot(length(j), length(k), i);
            a = as(ji,ki);
            b = bs(ji,ki);
            x = linspace(0, 3, 100);
            y = zeros(size(x));
            for xi = 1:length(y)
                y(xi) = haar(x(xi), a, b);
            end
            plot(x,y);            
            title(sprintf('j = %d, k = %d, a = %.2f, b = %.2f', j(ji), k(ki), a, b))
            i=i+1;
        end
    end
        
   


end

