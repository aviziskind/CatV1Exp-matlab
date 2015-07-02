function [Wt, E_tot, E_specific] = perceptron(X, Tid, loss_v, showWorking)
    
    
    if nargin < 4
        showWorking = false;
    end    

    [D,N] = size(X);
    Tid = Tid(:);
    if (length(Tid) ~= N)
        error('Number of columns of X must equal number of Ids');
    end
    
    [uTid, T_lists] = uniqueList(Tid);
    Tcounts = cellfun(@length, T_lists);
    
    K = length(uTid);    
    K_needed = iff(K > 2, K, 1);
    if K < 2
        error('Must have at least 2 classes');
    end        
    
    T = zeros(K_needed,N);
    for kk = 1:K_needed
        T(kk,Tid == kk) = 1;
        T(kk,Tid ~= kk) = -1;
    end    
    
    
    if ~exist('loss_v', 'var')
        loss_v = ones(K,1);  % treat errors in all categories equally (ie. symetrically)
    end        
    if isvector(loss_v)
        if length(loss_v) ~= K
            error('loss vector must be of length K');
        end
        loss_m = ones(K);
        for kk = 1:K
            loss_m(kk,[1:kk-1, kk+1:K]) = loss_v(kk);
        end

    elseif isamatrix(loss_v)
        loss_m = loss_v;
        if any(size(loss_m) ~= K)
            error('loss matrix must be of size KxK');
        end
    end

    Xt = [ones(1,N); X];

    L = [];
    if K >= 3
        L = zeros(K,N);
        for kk = 1:K
            L(kk,:) = loss_m( (kk-1)*K + Tid(:) );  % L(i) = loss_m(Tid(i), c(i));
        end
    end
    
    
    X_sorted = sort(X,2);
    X_diffs = diff(X_sorted,1,2);
    X_scale = mean( median(X_diffs, 2))*2;    
    
    
    
    if showWorking   % test mode
        cols = 'brg';
        switch D
            case 1
                figure(1); clf; hold on;
                for kk = 1:K
                    x = X(Tid == kk);                    
                    xs{kk} = x;                     %#ok<AGROW>
%                     plot(x1, x2, [cols(kk) '.']);
                end                      
                displayNHist(xs);
                
            case 2
                figure(1); clf; hold on;
                for kk = 1:K
                    x1 = X(1, Tid == kk);
                    x2 = X(2, Tid == kk);
                    plot(x1, x2, [cols(kk) '.']);
                end                      
            case 3
                figure(1); clf; hold on
                for kk = 1:K
                    x1 = X(1, Tid == kk);
                    x2 = X(2, Tid == kk);
                    x3 = X(3, Tid == kk);
                    plot3(x1, x2, x3, [cols(kk) '.']);
                end
        end
        
        if D > 1
            axis equal; 
            axis(axis);
        end
    end
    

    
    function c = classifyY(Y)
        
        if K == 2
            c = 2*ones(N,1);            
            c(Y > 0) = 1;
        end        
        
        if K >= 3
            assert(all(size(Y) == [K,N]));
            [tmp, c] = max(Y,[], 1);
        end        
        
    end
                        

    function E = errorfunc_tanh(Wt)
        
        Y = Wt' * Xt;   
        c = classifyY(Y); 
        
        if K == 2
            L = loss_m( (c(:)-1)*K + Tid(:) );  
        end
        diff = tanh(Y(:)/X_scale) - T(:);
        E = sum( L(:) .* (abs( diff ) ) );
    end

    pow = .3;
%     h = zeros(1,K);
    function [E_tot, E_specific] = errorfunc_lin(Wt)
        
        Y = Wt' * Xt;   
        c = classifyY(Y); 
        
        if K == 2
            L = loss_m( (c(:)-1)*K + Tid(:) );  % L(i) = loss_m(Tid(i), c(i));
        end
        
%         E = sum( L(:) .* (abs( Y(:) - T(:) ).^pow ) );
        diff = Y(:) - 4*T(:);
        E_tot = sum( L(:) .* (abs( diff  ).^pow ) );
            
        
        if showWorking
%             for k = 1:K
%                 idx = (Tid == k);
%                 h(k) = stem3orUpdate(h(k), Xt(2,idx), Xt(3,idx), Y(idx), cols(k));
%             end            
%             fprintf('E = %.3f\n', E)
%             3;
        end
        
        
    end


    function drawplane(wt, xlims, ylims)
        xs = linspace(xlims(1), xlims(2), 100);
        ys = linspace(ylims(1), ylims(2), 100);
        [xs_g, ys_g] = meshgrid(xs, ys);
        n = length(xs);
        X = [ones(1, n*n);  xs_g(:)';  ys_g(:)'];
        z = wt' * X;
        z(z > 2) = 2;
        z(z < -2) = -2;
        mesh(xs, ys, reshape(z, n,n));        
    end

    function [E_tot, E_specific] = actualError(Wt)

        Y = Wt' * Xt;   
        c = classifyY(Y);
        E_specific = zeros(K,1);
        err_idxs = find(Tid ~= c);
        for k = 1:K
            E_specific(k) = nnz(Tid(err_idxs) == k); %#ok<FNDSB>
        end
        E_tot = sum(E_specific);
        
    end


    function drawWt(Wt, col)
            
        switch D

            case 1  % 1-Dimensional case
                3;
                k0 = -Wt(1)/Wt(2);
                drawVerticalLine(k0, 'color', col)
                3;

            case 2  % 2-Dimensional case
                                
                x1s = linspace(min(X(1,:))-1, max(X(1,:)+1), 150);
                x2FromWt = @(wt, x1) (-wt(1)-wt(2)*x1)/wt(3);
    %             x2s0 = x2FromWt(Wt0, x1s);
                %     plot(x1s, x2s0, 'g.-');    
%                 drawplane(Wt, xlim, ylim);
                3;
                for k = 1:K_needed
                    x2s = x2FromWt(Wt(:,k), x1s);                    
                    plot(x1s, x2s, [col '-']);
                end
                
                if (K == 3)
                    for k = 1:K
                        [k1, k2] = dealV(setdiff(1:3, k));
                        Wt12 = Wt(:,k1)-Wt(:,k2);
                        x2s = x2FromWt(Wt12, x1s);
                        Xline = [ones(1,length(x1s)); x1s; x2s];
                        
                        Y = Wt' * Xline;
                        [tmp, inds] = max(Y, [], 1);
                        idx = (inds == k1) | (inds == k2);                        
                        
                        plot(x1s(idx), x2s(idx), ['k-']);                        
                    end
                end


            case 3
                x1s = linspace(min(X(1,:))-1, max(X(1,:)+1), 10);
                x2s = linspace(min(X(2,:))-1, max(X(2,:)+1), 10);
                [x1s_g, x2s_g] = meshgrid(x1s, x2s);

    %             x3 = -(w0 + w(1)*x(1) + w2.x2)/w3
                x3FromWt = @(Wt, x1, x2) -(Wt(1)+Wt(2)*x1+Wt(3)*x2)/Wt(4);

    %             x3s0 = x3FromWt(Wt0, x1s_g, x2s_g);
                %     plot(x1s, x2s0, 'g.-');    

                x3s = x3FromWt(Wt, x1s_g, x2s_g);
                surf(x1s, x2s, x3s, 'facealpha', .2, 'cdata', 1, 'facecolor', col);

        end
        
    end

    D_needed = D+1; %iff(D == 1, D, D+1);
    if K_needed == 1
        M1 = mean(X(:,T_lists{1}), 2);
        M2 = mean(X(:,T_lists{2}), 2);
        Wt0 = -getInitialGuess(M1,M2);
    else
        Wt0 = rand(D_needed,K_needed);
    end

    opt = optimset('maxfunevals', 5000, 'maxiter', 5000);
    
%     [Wt1, fval, exitflag(1)] = fminsearch(@errorfunc_lin,Wt0, opt);
%     Wt20 = fminsearch(@errorfunc_tanh,Wt);
    [Wt2, fval, exitflag(2)] = fminsearch(@errorfunc_tanh,Wt0, opt);
%     Wt1 = fminsearch(@errorfunc_lin,Wt2);
%     [Wt3, fval, exitflag(3)] = fminsearch(@actualError,Wt0, opt);

    
%     [E_tot0, E_specific0] = actualError(Wt0);
%     [E_tot1, E_specific1] = actualError(Wt1);
    [E_tot2, E_specific2] = actualError(Wt2);
%     [E_tot3, E_specific3] = actualError(Wt3);    
    if any(exitflag == 0);
        3;
    end
    
    [Wt, E_tot, E_specific] = deal(Wt2, E_tot2, E_specific2);

    E_tot = E_tot/length(X);
    E_specific = E_specific(:)./Tcounts(:);
    
    if showWorking
        indiv_temp = '[k%d : %d]';
        fprintf('Error for Wt1 : %d (%.3f%%), %s\n', E_tot1, E_tot1*100/N, sprintf(indiv_temp, [1:K; E_specific1(:)']) );

        fprintf('Error for Wt2 : %d (%.3f%%), %s\n', E_tot2, E_tot2*100/N, sprintf(indiv_temp, [1:K; E_specific2(:)']) );

        fprintf('Error for Wt3 : %d (%.3f%%), %s\n', E_tot3, E_tot3*100/N, sprintf(indiv_temp, [1:K; E_specific3(:)']) );
    
%         drawWt(Wt1, 'b');
        drawWt(Wt2, 'g');
%         drawWt(Wt3, 'k');

        %%%
        viewErrFunc = @actualError;
        Wt = Wt2;
        %%%
    
%     W = Wt(2:end);
%     nsamples = 50;
%     x0 = -Wt(1)/mag(W);
%     xs_samples = linspace(x0-X_scale*50, x0+X_scale*50, nsamples);
%     errs = zeros(K+1,nsamples);
%     
%     
%     for i = 1:nsamples
%         wt_sample = Wt;
%         wt_sample(1) = -xs_samples(i)*  norm(W);        
%         [errs(K+1,i), errs(1:K,i)] = actualError(wt_sample);
%         drawWt(wt_sample, 'k');    
%     end
% %     [E_tot, E_spec] = viewErrFunc(
%     figure(2); clf;    
%     plot(xs_samples, errs(1,:), 'b.-', xs_samples, errs(2,:), 'g.-', xs_samples, errs(3,:), 'r.-');  hold on;
%     plot(x0, E_specific(1), 'bo', x0, E_specific(2), 'go', x0, E_tot, 'ro'); 
%     
%     legend('1', '2', 'tot');
    
    end
end



function W = getInitialGuess(A,B)
    normVec = B-A;    
    midPoint = (A+B)/2; 
    w0 = -sum(normVec .* midPoint);
    W = [w0; normVec(:)];    
end

%     function plotW(Wt)
%         switch D
%             case 2
% 
%                 switch K
%                     x1s = linspace(min(X(1,:))-1, max(X(1,:)+1), 15);
%                     case 2,
%                         x2FromWt = @(Wt, x1) (-Wt(1)-Wt(2)*x1)/Wt(3);
%                         x2s = x2FromWt(Wt1(:,k), x1s);
%                         plot(x1s, x2s, [cols(k) '-']);
%                     
%                     case {3,4,5}
%                         for k = 1:K
%                             x2s = x2FromWt(Wt1(:,k), x1s);
%                             plot(x1s, x2s, [cols(k) '-']);
%                         end
% %                         if (K == 3)
% %                             for k = 1:K
% %                                 [id1, id2] = setdiff(1:3, k);
% 
%                     
%                     
%                         
%                 end
%                 
%                 
%                 
%             case 3
%                 
%                 
%         end
%         
%     end


%     function h = stem3orUpdate(h, xdata, ydata, zdata, varargin)
%         if isempty(h) || (h <= 0)
%             h = stem3(xdata, ydata, zdata, varargin{:});
%         else
%             set(h, 'xdata', xdata', 'ydata', ydata, 'zdata', zdata);
%         end        
%     end
