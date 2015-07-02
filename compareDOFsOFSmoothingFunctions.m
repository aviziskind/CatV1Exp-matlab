function compareDOFsOFSmoothingFunctions

    nPh = 8;
    
    doGauss = true;
    doFermi = false;
    doAverage = false;
    
    if doGauss
        figure(51); clf; hold on;
        switch nPh
            case 8, ws = [0:.1:3];
            case 60, ws = [0:.01:4, 4.1:.1:15, 16:1:40];
        end
        progressBar('init-', length(ws)*2, 30)
        cols = 'br';
        gauss_dofs = zeros(2,length(ws));
        for i = 1:2
            circFlag = iff(i==1, 1, []);            
            for wi = 1:length(ws)
                progressBar;                
                smoothFunc = @(x) gaussSmooth(x, ws(wi), [], circFlag);
                gauss_dofs(i,wi) = getDOFofSmoothingFunction(smoothFunc, nPh, 2000);
            end
            plot(ws, gauss_dofs(i,:), [cols(i) '.-']);
        end        
        progressBar('done');                
        
        figure(60); clf; hold on;
        idx = find(ws > 1);
        w1 = ws(idx);
        g1 = gauss_dofs(1,idx);
%         gaussWtoDOF = @(b, x) max(min( [b(1)./((x-b(2)).^((x).^b(3) + b(4) ) ) + b(5)], nPh*ones(size(x)) ), 3*ones(size(x))) ;

        
        switch nPh
            case 8, beta0 = [12.9415   -1.3612    2.0189    0.7425    2.9712];  
                    gaussWtoDOF = @(b, x) max(min( [b(1)./((x-b(2)).^(x.*b(3) + b(4))) + b(5)], nPh*ones(size(x)) ), 3*ones(size(x))) ;
            case 60, beta0 = [87.2120   -0.6554   -1.1062    1.6746    2.0664];
                    gaussWtoDOF = @(b, x) max(min( [b(1)./((x-b(2)).^((x).^b(3) + b(4) ) ) + b(5)], nPh*ones(size(x)) ), 3*ones(size(x))) ;                    
        end
        
        beta = nlinfit(w1, g1, gaussWtoDOF, beta0);
%         beta = nlinfit(w1, g1, gaussWtoDOF, [10 -1.5 2 3]);
        plot(ws, gauss_dofs(1,:), ['b.']); hold on;
        fplot(@(x) gaussWtoDOF(beta, x), [0, ws(end)], 'r')
        
        % nPh = 8:  beta = [12.9415  -1.3612  2.0189  0.7425  2.9712]
        % nPh = 60:  beta = [12.9415  -1.3612  2.0189  0.7425  2.9712]
        
        legend({'circular', 'non-circular'})
        title('DOF''s of tuning curves vs width of gaussian smoothing');
        xlabel('std of gaussian'); ylabel('degrees of freedom')
        3;
    end
    
    if doFermi
        figure(52); clf; hold on;
        ks = [1:60];        
        progressBar('init-', length(ks), 30)

        fermi_dofs = zeros(size(ks));
        for ki = 1:length(ks)
            progressBar;           
            k = ks(ki); b = k/50;
            smoothFunc = @(x) fermiSmooth(x, k, b, 1);
            fermi_dofs(ki) = getDOFofSmoothingFunction(smoothFunc, nPh, 20000);
        end
        progressBar('done');                
        plot(ks, fermi_dofs, 'o-');

        title('DOF''s of tuning curves vs frequency bandpass (k_{hp})');
        xlabel('k_{hp}'); ylabel('degrees of freedom')
    end

    
    if doAverage
        figure(53); clf; hold on;
        ns = [2:1:60];        
        progressBar('init-', length(ns), 30)

        average_dofs = zeros(size(ns));
        for ni = 1:length(ns)
            progressBar;           
            n = ns(ni); 
            smoothFunc = @(x) alias(x, n, 1);
            average_dofs(ni) = getDOFofSmoothingFunction(smoothFunc, nPh, 50000);
        end
        progressBar('done');                
        plot(ns, average_dofs, 'o-');

        title('DOF''s of tuning curves vs length after aliasing');
        xlabel('n'''); ylabel('degrees of freedom')
    end    
    
    
    3;
end