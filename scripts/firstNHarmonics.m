function [y, t_itp, y_itp] = firstNHarmonics(y, n, itpM, itpMethod, t)

    % keeps the first 'n' harmonics of the vector y, and outputs back to y.
    % (if n is empty, all the harmonics are kept)
    % optionally if more than 2 input arguments, also interpolates between
    % the points, using itpM-point interpolation.

    if isvector(y)
        y = y(:);
    end    
    N = size(y,1);
    Y = fft(y);
        
%     idx_half = round(N/2);
    if ~isempty(n)
        idx_to_blank = [n+2: N-n];
        Y(idx_to_blank,:) = 0;
    end
    y = ifft(Y, 'symmetric');
            
    if (nargout > 1)
        if ~exist('itpMethod', 'var') || isempty(itpMethod)
            itpMethod = 'harmonics';
        end
        if ~exist('itpM', 'var') || isempty(itpM)
            itpM = 20;
        end        
        if ~exist('t', 'var') || isempty(t)
            t = [0:N-1]';
        end        

        if itpM == 1
            t_itp = t;
            y_itp = y;
            return;            
        end
        dt = diff(t(1:2));        
        t_ext = [t(:); t(end)+dt];
        
        switch itpMethod
            case 'harmonics'                
                t_itp_ext = [0 : 1/itpM : N]*dt;
                t_itp = t_itp_ext(1:end-1);                                                

%                 y_itp = zeros(size(t_itp));
%                 for j = 0:length(Y)-1
%                     y_itp = y_itp + exp(1i * j * t_itp *2*pi/N) * Y(j+1)/N;
%                 end
%                 if (sum(abs(imag(y_itp))) < 1e-10)
%                     y_itp=real(y_itp);
%                 end
                
                mid = round(N/2);
                if odd(N)
                    Y_itp = [Y(1:mid); zeros(length(Y)*(itpM-1), 1); Y(mid+1:end)]*itpM;                                    
                else
                    Y_itp = [Y(1:mid); Y(mid+1)/2; zeros(length(Y)*(itpM-1)-1, 1); Y(mid+1)/2; Y(mid+2:end)]*itpM;                
                end
                y_itp = ifft(Y_itp, 'symmetric');                
                
                t_itp = t_itp + t(1);
                
%                 y_itp=abs(y_itp);
%                 3;

            case 'linear'
                t_itp_ext = [0 : 1/itpM : N]*dt;
                t_itp = t_itp_ext(1:end-1);                                
                y_wrp = [y(:); y(1)];
                y_itp = interp1q(t_ext, y_wrp, t_itp');                
                
            case 'spline' % this ends up being almost exactly the same as the harmonics method.
                t_itp = [0 : 1/itpM : 3*N]*dt;                
                t3 = [0:3*N-1]'*dt;
                y3 = [y;y;y];  
                idx_mid = itpM*N+1:2*itpM*N;
                
                y_itp = interp1(t3, y3, t_itp', 'spline');                

                t_itp = t_itp(idx_mid); t_itp = t_itp-t_itp(1);
                y_itp = y_itp(idx_mid);
                                
        end
        
    end
    
    
end
    
   %{ 
    
    
    
    
    % 1. determine dt.
    if (length(t) == 1)  % allow for t = dt
        dt = t;        
        t = (0:length(f)-1)*dt;
    elseif (length(t) ~= length(f))
        error('t must either be of length 1 or length f');
    end       

    % 2. determine T.
    if nargin < 3
        T = t(end) + diff(t(1:2)); % even if t was initially scalar, it is a vector now.
    end

    % rescale time period from 0:T to 0:2pi
    w = (2*pi)/T;
    t = t(:)*w;
    f = f(:);
    dt = t(2) - t(1);
    dts = [diff(t); dt];
    nCyc = (t(end)+dt-t(1))/(2*pi);

    % F1 component
    

    

    
% 
    fft_1 = fourierTransform(1, dts, t, f, T);
%     fft1 = fourierTransform(1, dts, t, f, T);
%     assert(abs(fft_1) == abs(fft1));

    fft_2 = fourierTransform(1, dts, t, f, T);
    fft2 = fourierTransform(1, dts, t, f, T);
    assert(abs(fft_2) == abs(fft2));
    
    
    F1 = 2*sqrt(2)* abs(fft_1);
    DC = fourierTransform(0, dts, t, f, T);

    if (F1 < 1e-10*DC)
        phi = NaN;
    else
        phi = angle(fft_1);
        phi = mod(phi, 2*pi);
    end        
    
    % convert negative phases (-pi:0) to positive (pi:2*pi);
           
    
    
    
end

% function F = fourierTransform(k, dts, t, f, T)  
%     y = 1/sqrt(T) * sum( dts .* f .* exp(-1i * 2*pi* k * t/T) );
% end
% function f = invFourierTransform(w, dts, t, F, T)  
%     y = 1/sqrt(T) * sum( dts .* F .* exp(1i * 2*pi* k * t/T) );
% end

% function F = fourierTransform(f, j)  
%     N = length(f);
%     k = [0:N-1];
%     if isempty(j)
%         j 
%     F = 1/sqrt(N) * sum( exp(-2*pi*1i*j*k/N) .* f );
%     
% end
% function f = invFourierTransform(F, j)  
%     N = length(F);
%     k = [0:N-1];
%     f = 1/sqrt(N) * sum( exp(2*pi*1i*j*k/N) .* F );
% end



% 
% 
%     f = f(:);
%     if (length(t) == length(f))
%         dt = t(2) - t(1);
%         t = t(:);
%     elseif (length(t) == 1)
%         dt = t;
%         t = [(0:length(f)-1)*dt]';
%     end    
%     if nargin < 3
%         L = T;%t(end);
%     else
%         L = t(end)-t(1);
%     end
%     % rescale T to 0:2pi
% %     t = t(:)*(2*pi)/T;
% 
%     dts = [diff(t); dt];
%     
%     
%     % DC component
%     DC = abs( sum(f .* dts) ) / L;
% 
%     % F1 component
%     fourierTransform = @(w)   (1/sqrt(2))* sum( dts .* exp(1i * w * t) .* f) /L;
%     f1 = fourierTransform(1);
%     f_1 = fourierTransform(-1);
%     c = 2;
%     F1 = c * sqrt( (abs(f1))^2 + (abs(f_1))^2 );

%}