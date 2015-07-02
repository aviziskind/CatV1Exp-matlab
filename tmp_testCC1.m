function tmp_testCC1

    Nt = 1000;
    n = 10;
    b = 2;
    
    
    xs1 = zeros(n, Nt);
    xs2 = zeros(n, Nt);
    ccs = zeros(1, Nt);
    for i = 1:Nt        
        x1 = randvec(n,b);
        x2 = randvec(n,b);
        cc = pearsonR(x1, x2);
        if cc < -.1
            figure(4);
            plot(1:n, x1, 'bo-', 1:n, x2, 'g.-')
            
            3;
        end
            
        
        xs1(:,i) = x1;
        xs2(:,i) = x2;   
        ccs(i) = cc;        
    end
%     ccs = pearsonR_v(xs1, xs2);
    
    figure(5); clf;
    hist(ccs, 30);
    


end

function y = randvec(n,b)
    y = zeros(1,n);    
    idxs = randperm(n); 
    idxs = idxs(1:b);
    vals = abs(randn(1,b));
    y(idxs) = vals;
end