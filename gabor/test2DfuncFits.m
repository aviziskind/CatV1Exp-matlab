function test2DfuncFits
    
    xs = -1:.05:1;
    ys = -1:.05:1;
    [X, Y] = meshgrid(xs, ys);
    a = 5;
    b = 10;
    Z = simpFun2d(a,b,X,Y);
    Z = Z + randn(size(Z))*.3;
 
    % try using 'fit' algorithm
    
    ft = fittype('simpFun2d(a,b,x1,x2)', 'independent', {'x1', 'x2'});
    fo = fitoptions('method', 'nonlinearleastsquares', 'startpoint', [3 1]);
    f = fit( [X(:), Y(:)], Z(:), ft, fo );
    
    plot( f, [X(:), Y(:)], Z(:) );
    
    
    % try using nlinfit
    beta0 = [3 5];
    func = @(b, X) simpFun2d(b(1), b(2), X(:,1), X(:,2));    
    Be = nlinfit([X(:), Y(:)], Z(:), func, beta0);
    
    
    3;
    
    
end


% function y = simpFun2d_int(a, b, X)
%     y = simpFun2d(a, b, X(:,
% 
% end


%{
 
%}
 


