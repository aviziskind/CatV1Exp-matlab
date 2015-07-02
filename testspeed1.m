% function testspeed1
%     t = 0:.1:10000;
%     tic;
%     y = sin(t);
%     toc;
%     tic;
%     y = arrayfun(@(T) sin(T), t);
%     toc;
    
% end

    N = 1000000;
    
    nB = 10; nC = 10;
    B = randi(nB,1,N);
    C = randi(nC,1,N);
    
    tic;
    X1 = zeros(nB,nC);
    for b = 1:nB
        for c = 1:nC
            X1(b,c) = nnz( (B == b) & (C == c));            
        end
    end
    toc;
    t1 = toc;    
    
    
    tic;
    X2 = zeros(nB,nC);
    for i = 1:N
        X2(B(i),C(i)) = X2(B(i),C(i))+1;        
    end
    toc;
    t2 = toc;
    
    assert(all(X1(:) == X2(:)));
%     tic;
%     X3 = zeros(nB,nC);
%     for i = 1:N
%         X2(B(i),C(i)) = X2(B(i),C(i))+1;        
%     end
%     toc;
    disp(num2str(t1/t2));
    

    