function testGenerateListOfCellPairs
    
    ALL_PAIRS_TEST = 1;
    SUB_RANGES_TEST = 2;
    RAND_PAIRS_TEST = 3;
    
    testsToDo = [ALL_PAIRS_TEST, SUB_RANGES_TEST]; % RAND_PAIRS_TEST

    nSites = 20;
    nCellsPerSite = randi(5,1,nSites)+1;  
%     nCellsPerSite = [3 5 3 0 4];
    cellIds = arrayfun(@(n) [iff(rand<0, [], 0) , 1:n], nCellsPerSite, 'un', false);
    
    nUnits= length([cellIds{:}]);

    A = zeros(nUnits);
    for i = 1:nUnits
        A(i,i) = 1;
    end    
    A(1,1) = -1;
    DIAG = A;
    
    if any(testsToDo == ALL_PAIRS_TEST)
    
        [WccPairs, WcmPairs, BccPairs, BcmPairs] = generateListOfCellPairs(cellIds, 1, [], [], 'orig');
    
    %     allPairs = sort([WccPairs; WcmPairs; BccPairs; BcmPairs]);
    %     idxMtx = zeros(nUnits, nUnits);
    %     idxMtx(sort(allPairs)) = 1:length(allPairs);    
    %     figure(44);
    %     imagesc(idxMtx);
    %     3;
    
        A( WccPairs ) = .5;
        A( WcmPairs ) = -.5;
        A( BccPairs ) = 1;
        A( BcmPairs ) = -1;    
    
        figure(1);
        imagesc(A); axis equal tight;
        title('All Pairs');
a    end
    
    if any(testsToDo == SUB_RANGES_TEST)
    
        rngs = {[1:4], [5:8], [9:15]}; 
        B = DIAG;
        [WccPairs_sub, WcmPairs_sub, BccPairs_sub, BcmPairs_sub] = generateListOfCellPairs(cellIds, 1, rngs, [], 'orig');
        B(WccPairs_sub ) = .5;
        B(WcmPairs_sub) = -.5;
        B(BccPairs_sub ) = 1;
        B(BcmPairs_sub) = -1;

        figure(2);
        imagesc(B); axis equal tight;
        title('Pairs with ranges');
        

    end
    
    
    if any(testsToDo == RAND_PAIRS_TEST)
        nResamples = 1000;
        [WccPairs_r, WcmPairs_r] = generateListOfCellPairs(cellIds, nResamples);
        C = DIAG;
        D = DIAG;
        [uWcc_r, wcc_count] = uniqueCount([WccPairs_r{:}]);
        [uWcm_r, wcm_count] = uniqueCount([WcmPairs_r{:}]);
        
        C(uWcc_r) = 1;
        C(uWcm_r) = -1;

        D(uWcc_r) = wcc_count;
        D(uWcm_r) = wcm_count;

        figure(3);
        imagesc(C); axis equal tight;
        title('Randomized pairs: coverage');

        figure(4);
        imagesc(D); axis equal tight;
        title('Randomized pairs: multiplicity');
                colorbar;
        
    end
end
    
    
    
    
    
    
    
    