function B = runLengths(A, inverseFlag)
    if exist('inverseFlag', 'var') && ~isempty(inverseFlag) % inverse operation: from runlengths -> orig matrix
        B = zeros(1, sum(A(:,2)) );
        b_ind = 1;
        for i = 1:size(A,1)
            [val, n] = elements( A(i,:) );
            B(b_ind:b_ind+n-1) = val;
            b_ind = b_ind+n;
        end
        return;        
    end

    n = length(A);
    B = zeros(n,2);
    indA = 1;
    indB = 1;
    while (indA <= n)
        B(indB,1) = A(indA);
        count = 1;
        while (indA+count <= n) &&  (A(indA) == A(indA+count))
            count = count + 1;
        end
        B(indB,2) = count;        
        
        indB = indB + 1;
        indA = indA + count;
    end
    B(indB:end,:) = [];
end