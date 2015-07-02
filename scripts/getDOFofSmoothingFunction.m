function dof = getDOFofSmoothingFunction(smoothFunc, n, nTrials)

    % should be in the syntax: r_sm = smoothFunc(r)    
    if nargin < 3        
        nTrials = 2000;
    end
    std2dof = @(st)  1./(st).^2 + 1;    
    ccs = getCCs(nTrials, n, smoothFunc);
    cc_std = std(ccs);  
    
    dof = std2dof(cc_std);
end



function std2dof = getStd2DofCurve %#ok<DEFNU>
    showWorking = true;
                
    dof2std = @(beta, x) beta(1)./((x-beta(2)).^beta(3));

    Ls = [2:100, 105:5:500];%, 25:5:500]; % parameter values
    nTrials = 10000;
    ccs = getCCs(nTrials, Ls);
    cc_stds_vs_n = std(ccs, [], 1);        
    
%         beta = nlinfit(Ls, cc_stds_vs_n, dof2std, [1 1 1]);
    beta = [1 1 .5];

    std2dof = @(beta, x)  (beta(1)./x).^(1./beta(3)) + beta(2);        

    if showWorking
        figure(45); clf;
        plot(Ls, cc_stds_vs_n, '.'); hold on;
        axis([0 Ls(end)+1, 0 1.1]);
        fplot(@(x) dof2std(beta, x), [1.1, Ls(end)], 'r')
        xlabel('n (length of vectors)')
        ylabel('std of cc distribution')
        legend('simulations', 'y = 1/sqrt(x-1)')        
    end
   
end


function ccs = getCCs(nTrials, Ls, smoothFunc)
    
    nL = length(Ls);
    if (nargin < 3) || isempty(smoothFunc)
        smoothFunc = @(x) x;
    end
    
    ccs = zeros(nTrials, nL);

    for l_i = 1:nL
        L = Ls(l_i);

        X = randn([L,nTrials]);
        Y = randn([L,nTrials]);

        for t = 1:nTrials
            x = smoothFunc( X(:,t) );
            y = smoothFunc( Y(:,t) );            
            ccs(t, l_i) = pearsonR(x, y);            
        end

%         X = smoothFunc( randn([L,nTrials]) );
%         Y = smoothFunc( randn([L,nTrials]) );
% 
%         for t = 1:nTrials
%             x = X(:,t);
%             y = Y(:,t);            
%             ccs(t, l_i) = pearsonR(x, y);            
%         end

    end
    
end