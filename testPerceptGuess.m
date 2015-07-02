function testPerceptGuess

    n1 = 100;
    n2 = 100;
    D = 2;
    
    O1 = randi(10,2,1)-5;
    O2 = randi(10,2,1)-5;
    X1 = bsxfun(@plus, randn(D,n1), O1);
    X2 = bsxfun(@plus, randn(D,n2), O2);
    X = [X1, X2];
    figure(1); clf; hold on;
    plot(X1(1,:), X1(2,:), 'b.'); 
    plot(X2(1,:), X2(2,:), 'g.');
    axis equal;
    v = axis;
    
    M1 = mean(X1,2);
    M2 = mean(X2,2);
    
    W = getInitialGuess(M1,M2);

    x1s = linspace(min(X(1,:))-1, max(X(1,:)+1), 150);
    x2FromWt = @(wt, x1) (-wt(1)-wt(2)*x1)/wt(3);

    x2s = x2FromWt(W, x1s);                    
    plot(x1s, x2s, 'r:');
    axis(v);

end



function W = getInitialGuess(A,B)
    normVec = B-A;    
    midPoint = (A+B)/2; 
    w0 = -sum(normVec .* midPoint);
    W = [w0; normVec(:)];    
end
