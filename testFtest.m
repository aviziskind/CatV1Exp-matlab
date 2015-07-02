function testFtest


    show = 1;

    as = [0:.01:1];
    ps = zeros(size(as));

    
    randn('state', 0);
    x = [0:.5:10];
    y_sqr = x.^2;
    y_lin = 2.*x;
    y_rnd = randn(1, length(x))*2;
    
    for ai = 1:length(as)
        a = as(ai);

        y = y_lin*a + (1-a)*y_sqr + y_rnd;

        [poly1, S1] = polyfit(x,y,1);
        [poly2, S2] = polyfit(x,y,2);

        y1 = polyval(poly1, x);
        y2 = polyval(poly2, x);
        
        if show
            figure(45); clf;
            plot(x,y,'.');
        
            hold on
            plot(x, y1, 'r:');
            plot(x, y2, 'g:');
        end
        
        p1 = nestedFtest(S1, S2);
        p2 = nestedFtest(x, y, @(x) polyval(poly1, x), @(x) polyval(poly2, x), 1, 2);
        p3 = nestedFtest(x, y, y1, y2, 1, 2);

        assert(abs(p1-p2) < 1e-10)
        assert(abs(p1-p3) < 1e-10)
        ps(ai) = p1;
        

        
%         SS2=S2.normr^2; SS1 = S1.normr^2;
%         DF2 = S2.df;    DF1 = S1.df;
% 
%         fstat = ((SS1-SS2)/SS2) / ((DF1-DF2)/DF2);
%         p = 1-fcdf(fstat, DF1-DF2, DF2);
%         ps(ai) = p;

    end 

    3;

end