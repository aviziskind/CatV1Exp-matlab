function testAllCombos(a,b,c)
    X = {1, a, 1-a, b, 1-b, c, 1-c, sqrt(a), sqrt(1-a), sqrt(b), ...
         sqrt(1-b), sqrt(c), sqrt(1-c), a.^2, (1-a).^2, b.^2, (1-b).^2, c.^2, (1-c).^2 };
    n = length(X);
    for i = 1:n
       for j = 1:n
            Y = X{i}./X{j};
            idx_max = indmax(Y);
            if ibetween(idx_max, 2, 6)
                3;
                figure(55);
                plot(Y, 'o-');
            end        
       end
    end
    
    
end
    