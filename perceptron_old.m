function perceptron(varargin)

%     showWorking = true;

    if nargin == 0   % test mode
        D = 2;
        K = 2;
        N = 1000;
%     randn('state', 0);
        for k = 1:K
            mu(k,:) = randi(10, 1,D); %#ok<AGROW>
            sig(k,:) = randi(3, 1,D); %#ok<AGROW>
        end
%         mu1 = [-2,3]; sig1 = [1,1];
%         mu2 = [2,1]; sig2 = [2,1];
        Tops = [-1, 1];
        Tid = randi(K,N,1);
        T = Tops(Tid);
        X = zeros(N,D);
        for i = 1:N
            for d = 1:D
                X(i,d) = mu( Tid(i), d ) + randn * sig( Tid(i), d );
            end
        end
        
%         x1a = [sig1(1) * randn(N1,1) + mu1(1)];
%         x2a = [sig1(2) * randn(N1,1) + mu1(2)];
% 
%         x1b = [sig2(1) * randn(N2,1) + mu2(1)];
%         x2b = [sig2(2) * randn(N2,1) + mu2(2)];
        cols = 'brg';
        switch D
            case 2
                figure(1); clf; hold on;
                for k = 1:K
                    x1 = X(Tid == k, 1);
                    x2 = X(Tid == k, 2);
                    plot(x1, x2, [cols(k) '.']);
                end                      
            case 3
                figure(1); clf; hold on
                for k = 1:K
                    x1 = X(Tid == k, 1);
                    x2 = X(Tid == k, 2);
                    x3 = X(Tid == k, 3);
                    plot3(x1, x2, x3, [cols(k) '.']);
                end
        end

        axis equal; 
        ax = axis;
    
    end
%     return;
        
%     w0 = -5;
%     W = [-5; 3];
%     beta0 = [w0; W];
    beta0 = rand(D+1,1);
    
    
    axis(ax);
%     zeroaxes;
    
   
    function E = errorfunc(beta)
        w0 = beta(1);
        W = beta(2:end);
%         W = rand(2,1); 
%         w0 = 0;
        Y = W' * X' + w0;

        E = sum( (Y(:) - T(:)).^2);
    end

    beta1 = fminsearch(@errorfunc,beta0);

    
    switch D
        case 2
            x1s = linspace(min(X(:,1))-1, max(X(:,1)+1), 15);
            x2FromBeta = @(beta, x1) (-beta(1)-beta(2)*x1)/beta(3);    
%             x2s0 = x2FromBeta(beta0, x1s);
            %     plot(x1s, x2s0, 'g.-');    

            x2s = x2FromBeta(beta1, x1s);
            plot(x1s, x2s, 'g-');     
        case 3
            x1s = linspace(min(X(:,1))-1, max(X(:,1)+1), 10);
            x2s = linspace(min(X(:,2))-1, max(X(:,2)+1), 10);
            [x1s_g, x2s_g] = meshgrid(x1s, x2s);
            
%             x3 = -(w0 + w(1)*x(1) + w2.x2)/w3
            x3FromBeta = @(beta, x1, x2) -(beta(1)+beta(2)*x1+beta(3)*x2)/beta(4);
            
%             x3s0 = x3FromBeta(beta0, x1s_g, x2s_g);
            %     plot(x1s, x2s0, 'g.-');    

            x3s = x3FromBeta(beta1, x1s_g, x2s_g);
            surf(x1s, x2s, x3s, 'facealpha', .2, 'cdata', 1, 'facecolor', 'g');

    end
        
    
        
    
end