function testGetSpikeCrossCorrelogram

    doBasicTest = 1;
    doPathologicalCasesTest = 0;
    doCompareMatlabAndC = 1;
    doCompareSwitch = 0;

    
    if doBasicTest % Very basic test for c version
        % test cross correlogram        
        Y = [1:5];
        X = [.9, 1.1, 1.9, 2.1, 2.2,  2.8, 4.75, 10] ;
        maxDist = .5;

        cc_s = getSpikeCC_slow(X, Y, maxDist);
        cc_M = getSpikeCrossCorrelogram_Matlab(X, Y, maxDist);
        cc_c = getSpikeCrossCorrelogram(X, Y, maxDist);
        assert(isequal(cc_s, cc_M));
        assert(isequal(cc_M, cc_c));
        
        % test cross correlogram with index output
        [cc_s, idxX_s, idxY_s] = getSpikeCC_slow(X, Y, maxDist);
        [cc_M, idxX_M, idxY_M] = getSpikeCrossCorrelogram_Matlab(X, Y, maxDist);
        [cc_c, idxX_c, idxY_c] = getSpikeCrossCorrelogram(X, Y, maxDist);
        assert(isequal(cc_M, cc_s));
        assert(isequal(idxX_M, idxX_s));
        assert(isequal(idxY_M, idxY_s));

        assert(isequal(cc_M, cc_c));
        assert(isequal(idxX_M, idxX_c));
        assert(isequal(idxY_M, idxY_c));
    

        % test auto correlogram  : basic
        X = [.9, 1.1, 1.9, 2.1, 2.2,  2.8, 4.75, 10] ;
        maxDist = 5;
        cc_M = getSpikeCrossCorrelogram_Matlab(X, 'auto', maxDist);
        cc_c = getSpikeCrossCorrelogram(X, 'auto', maxDist);
        assert(isequal(cc_M, cc_c));

        
        % test auto correlogram  : regular
        X = sort(rand(1, 1000)*100);
        maxDist = 20;
        cc_M = getSpikeCrossCorrelogram_Matlab(X, 'auto', maxDist);
        cc_c = getSpikeCrossCorrelogram(X, 'auto', maxDist);
        assert(isequal(cc_M, cc_c));
        

        
        % test auto correlogram  : basic  : with index output
        X = [.9, 1.1, 1.9, 2.1, 2.2,  2.8, 4.75, 10] ;
        maxDist = 5;
        [cc_M, idxX1_M, idxX2_M] = getSpikeCrossCorrelogram_Matlab(X, 'auto', maxDist);
        [cc_c, idxX1_c, idxX2_c] = getSpikeCrossCorrelogram(X, 'auto', maxDist);
        assert(isequal(cc_M, cc_c));
        assert(isequal(idxX1_M, idxX1_c));
        assert(isequal(idxX2_M, idxX2_c));
        
        
        % test auto correlogram  : regular : with index output
        X = sort(rand(1, 1000)*100);
        maxDist = 20;
        [cc_M, idxX1_M, idxX2_M] = getSpikeCrossCorrelogram_Matlab(X, 'auto', maxDist);
        [cc_c, idxX1_c, idxX2_c] = getSpikeCrossCorrelogram(X, 'auto', maxDist);
        assert(isequal(cc_M, cc_c));
        assert(isequal(idxX1_M, idxX1_c));
        assert(isequal(idxX2_M, idxX2_c));

        
               
        
    end
    
    
    
    
    if doPathologicalCasesTest
        % 1. empty X, or empty Y;
        nxs = [0:2]; nys = [0:2]; maxDist = 1;
        for i = 1:length(nxs)
            for j = 1:length(nys)
                X = sort(randn(1,nxs(i))); Y = sort(randn(1,nys(j)));
                cc_c = getSpikeCrossCorrelogram(X, Y, maxDist);        
                cc_S = getSpikeCC_slow(X, Y, maxDist);        
                assert(isequal(cc_c(:), cc_S(:)));
            end
        end

        %2. in different ranges
        X = 1:5; Y = 100:105;
        cc_c = getSpikeCrossCorrelogram(X, Y, maxDist);        
        cc_S = getSpikeCC_slow(X, Y, maxDist);        
        assert(isequal(cc_c(:), cc_S(:)));

        %3. all out of range of each other
        X = 1:100:500; Y = 110:115; maxDist = 1;
        cc_c = getSpikeCrossCorrelogram(X, Y, maxDist);        
        cc_S = getSpikeCC_slow(X, Y, maxDist);        
        assert(isequal(cc_c(:), cc_S(:)));
                
    end
    
    
    if doCompareMatlabAndC % 1. test on sparse / dense X & Y
        N = 5000;
        maxDist = 5;
        getIdxs = 1;
        nOut = iff(getIdxs, 3, 1);
            
        
        X = [1:N] + randn(1,N); 
        Y = randn(1,1000)*sqrt(N)+N/2;
        X = sort(X); Y = sort(Y);

        tic;
        [vout_M1{1:nOut}] = getSpikeCrossCorrelogram_Matlab(X, Y, maxDist);
        t_M = toc;
        tic;        
        [vout_c1{1:nOut}] = getSpikeCrossCorrelogram(X, Y, maxDist);
        t_c = toc;
        assert(isequal(vout_M1, vout_c1));
        fprintf('Sparse : faster by %.2f\n', t_M/t_c);

        % 2. test on roughly equally spaced X & Y
        X = [1:N] + randn(1,N); 
        Y = [1:N] + randn(1,N); 
        X = sort(X); Y = sort(Y);

        tic;
        [vout_M2{1:nOut}] = getSpikeCrossCorrelogram_Matlab(X, Y, maxDist);
        t_M = toc;
        tic;
        [vout_c2{1:nOut}] = getSpikeCrossCorrelogram(X, Y, maxDist);
        t_c = toc;
        assert(isequal(vout_M2, vout_c2));

        fprintf('Dense : faster by %.2f\n', t_M/t_c);

        3;
    end
    
    if doCompareSwitch % assess speed increase when switch Y & X (should be no increase, because have implemented automatic switching of X & Y if ny > nx)
        N = 1000;
        X = [1:N] + randn(1,N); 
        Y = [1:10:N];% + randn(1,N/3); 
        X = sort(X); Y = sort(Y);
        maxDist = .5;

        B = 1000;
        tic;
        for b = 1:B
            cc_1 = getSpikeCrossCorrelogram(X, Y, maxDist);
        end
        t1 = toc;
        tic;
        for b = 1:B
            cc_2 = getSpikeCrossCorrelogram(Y, X, maxDist);
        end
        assert(isequal(sort(-cc_1(:)), cc_2(:)));
        t2 = toc;

        fprintf('X,Y : %.2f,   Y,X : %.2f,  faster by %.2f\n', t1, t2, t2/t1);
    end
    

    
    
    
    
end


function [cc, idxX, idxY] = getSpikeCC_slow(X,Y, maxDist)

    dists = bsxfun(@minus, Y(:), X(:)');
    idx = abs(dists) <= maxDist;
    if nargout > 1
        [idxY_v, idxX_v] = find(idx);        
        [cc, newOrd] = sort(dists(idx(:)));        
        idxX = idxX_v(newOrd);
        idxY = idxY_v(newOrd);        
    else
        cc = sort(dists(idx(:)));
    end
    
end



%     X = 1:2e4;
%     randn('state', 0);
%     Y = sort(X + randn(size(X)));
%     B = 10;
%     maxDist = 2;
%     ts = zeros(1,B);
%     
%     for i = 1:B
%         tic;
%         getSpikeCrossCorrelogram(X, Y, maxDist);        
%         ts(i) = toc;
%     end
%     fprintf('%.4f +- %.4f\n', mean(ts), stderr(ts));
%     
%     return;
