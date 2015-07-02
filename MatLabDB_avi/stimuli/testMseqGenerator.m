function testMseqGenerator

    % test 5 bits;
    nMseqBits = 14;
    ncols = 3;
    tapRegister = 123;
    
    msequenceGenerator('init', tapRegister, nMseqBits);
    n = 2^nMseqBits-1;
    N = ncols*n;
    
    A = msequenceGenerator([nMseqBits, N])';
    B = zeros(N,1);
    for i = 1:N
        row = num2str(A(i,:)')';
        B(i) = bin2dec(row);
    end
    B = reshape(B, n, ncols);
    disp(B);
    C = B(:,1);
    if ( length(C) == length(unique(C)) );
        disp('All rows unique');
    else
        disp('not unique');
    end
3;



end