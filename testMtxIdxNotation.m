% function testMtxIdxNotation
    N = 200;
    m = 15;
    
    doProjTest = 1;
    doDiagTest = 1;

    if doProjTest
        fprintf('Projection test\n');
        P = rand(N,1);
        A = P*P';

        % test 1 (vector x)
        x = rand(N,1);
        a1 = x'*A*x;

        a2 = 0;
        for I = 1:N
            for J = 1:N
                a2 = a2 + x(I)*A(I,J)*x(J);
            end
        end
        assert(abs(a1 - a2)<1e-5)

        n1 = 20;

        % test 2 (matrix x)
        X = rand(m,N);
        tic;
        for aa = 1:n1
            B1 = X*A*X';
        end
        t1 = toc;
        t1 = t1/n1;

        tic;
        B2a = zeros(m,m);
        for i = 1:m
            for j = 1:m            
                B2a(i,j) = X(i,:) * A * X(j,:)';                        
            end
        end
        t2a = toc;
        assert(all(abs(B1(:)-B2a(:))<1e-5));

        tic;
        B2b = zeros(m,m);
        for I = 1:N
            for J = 1:N
                B2b = B2b +  X(:,I) * A(I,J) * X(:,J)';                        
            end
        end
        t2b = toc;
        assert(all(abs(B1(:)-B2b(:))<1e-5));

        tic;
        B3 = zeros(m,m);
        t3 = 0;
    %     for i = 1:m
    %         for j = 1:m            
    %             
    %             for I = 1:N
    %                 for J = 1:N
    %                     B3(i,j) = B3(i,j)  + X(i,I) * A(I,J) * X(j,J);
    %                 end
    %             end
    %                 
    %         end
    %     end
    %     t3 = toc;    
    %     assert(all(abs(B1(:)-B3(:))<1e-5));
    %     return;

        nC = 20;
        tic;
        for i = 1:nC
            B4a = mtxMult_proj(X, P);
        end
        t4a = toc;
        t4a = t4a/nC;
        assert(all(abs(B1(:)-B4a(:))<1e-5));

        tic;
        for i = 1:nC
            B4b = mtxMult_proj2(X, P);
        end
        t4b = toc;
        t4b = t4b/nC;
        assert(all(abs(B1(:)-B4b(:))<1e-5));

        tic;
        for i = 1:nC
            B4c = mtxMult_proj3(X, P);
        end
        t4c = toc;
        t4c = t4c/nC;
        assert(all(abs(B1(:)-B4c(:))<1e-5));


        t_ref = t2a;
    %     t_ref = t1;
        fprintf('%.2f, vs  %.2f, and %.2f. C1 : %.2f, C2 : %.2f, C3 : %.2f\n', t2a/t_ref, t2b/t_ref, t3/t_ref, t4a/t_ref, t4b/t_ref, t4c/t_ref)

    end
    
    
    

    if doDiagTest
        fprintf('Diagonal test\n');
        D = rand(N,1);
        A = diag(D);

        % test 1 (vector x)
        x = rand(N,1);
        a1 = x'*A*x;

        a2 = 0;
        for I = 1:N
            for J = 1:N
                a2 = a2 + x(I)*A(I,J)*x(J);
            end
        end
        assert(abs(a1 - a2)<1e-5)

        n1 = 20;

        % test 2 (matrix x)
        X = rand(m,N);
        tic;
        for aa = 1:n1
            B1 = X*A*X';
        end
        t1 = toc;
        t1 = t1/n1;

        tic;
        B2a = zeros(m,m);
        for i = 1:m
            for j = 1:m            
                B2a(i,j) = X(i,:) * A * X(j,:)';                        
            end
        end
        t2a = toc;
        assert(all(abs(B1(:)-B2a(:))<1e-5));

        tic;
        B2b = zeros(m,m);
        for I = 1:N
            B2b = B2b +  X(:,I) * A(I,I) * X(:,I)';                        
        end
        t2b = toc;
        assert(all(abs(B1(:)-B2b(:))<1e-5));

        tic;
        B3 = zeros(m,m);
        t3 = 0;
    %     for i = 1:m
    %         for j = 1:m            
    %             
    %             for I = 1:N
    %                 for J = 1:N
    %                     B3(i,j) = B3(i,j)  + X(i,I) * A(I,J) * X(j,J);
    %                 end
    %             end
    %                 
    %         end
    %     end
    %     t3 = toc;    
    %     assert(all(abs(B1(:)-B3(:))<1e-5));
    %     return;

        nC = 20;
        tic;
        for i = 1:nC
            B4a = mtxMult_diag(X, D);
        end
        t4a = toc;
        t4a = t4a/nC;
        assert(all(abs(B1(:)-B4a(:))<1e-5));

        tic;
        for i = 1:nC
            B4b = mtxMult_diag(X, D);
        end
        t4b = toc;
        t4b = t4b/nC;
        assert(all(abs(B1(:)-B4b(:))<1e-5));

%         tic;
%         for i = 1:nC
%             B4c = mtxMult_proj3(X, D);
%         end
%         t4c = toc;
%         t4c = t4c/nC;
%         assert(all(abs(B1(:)-B4c(:))<1e-5));

        t_ref = t2b;
    %     t_ref = t1;
        fprintf('%.2f, vs  %.2f, and %.2f. C1 : %.2f, C2 : %.2f,\n', t2a/t_ref, t2b/t_ref, t3/t_ref, t4a/t_ref, t4b/t_ref)
        fprintf('C1/C2 = %.2f\n', t4a/t4b)

    end
