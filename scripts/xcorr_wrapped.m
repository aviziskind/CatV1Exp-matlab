function cc = xcorr_wrapped(x, y, nmax)

    wrap = true;
    
    Nx = length(x);
    if (length(y) ~= Nx)
        error('x and y must be the same lngth');
    end
    if nargin < 3
        nmax = Nx;
    end    
    
    Nc = min(Nx, nmax);
    
    cc = zeros(Nc, 1);
    mu_x = mean(x);
    mu_y = mean(y);
    sig_x = sqrt(var(x, 1));
    sig_y = sqrt(var(y, 1));
    
    x = (x-mu_x);
    y = (y-mu_y);
    
    inds_x  = [1:Nc]';
           
    for xi = 1:Nc 
        
        inds_y = inds_x+(xi-1);
        if wrap
            inds_y = mod(inds_y-1, Nx)+1;            
            X = x(inds_x);
            Y = y(inds_y);
        else
            idx = inds_y < Nx;            
            X = x(inds_x(idx));
            Y = y(inds_y(idx));
        end
                                
        cc(xi) = mean(X .* Y)/(sig_x*sig_y);
        
    end
    
    
%     if wrap
%         if odd(Nc)
%             mid = floor(Nc/2);            
%         else
%             mid = round(Nc/2);
%         end
%         cc = cc( [mid+1:end, 1:mid]);        
%     end
    
%     length(cc)
end
