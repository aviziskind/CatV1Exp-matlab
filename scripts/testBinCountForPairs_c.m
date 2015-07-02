% rand('state', 0);
nData = 1000000;

pairIdxs = uint32(unique( randi(nData, 1, nData/2)));
mask = rand(1, nData) <= .7;
% nBins = 5; nMaxCats = 3;
nBins = 10; nMaxCats = 5;

binIds = uint8(randi(nBins,1, nData));
catIds = uint8(randi(nMaxCats,1, nData));
% catIds = [];
vals = rand(1,nData);

% test with 1 output
binN1 = binCountForPairs_Matlab(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
binN2 = binCountForPairs_c(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
assert(isequal(binN1, binN2));

% test with 2 outputs
[binN1,valsOut1] = binCountForPairs_Matlab(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
[binN2,valsOut2] = binCountForPairs_c(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
assert(isequal(binN1, binN2));
assert(isequal(valsOut1, valsOut2));

% test with 3 outputs (and compare how long they each take)
tic;
[binN1,valsOut1, idx1] = binCountForPairs_Matlab(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
t1 = toc;
tic;
[binN2,valsOut2, idx2] = binCountForPairs_c(pairIdxs, mask, binIds, catIds, nBins, nMaxCats, vals);
t2 = toc;
assert(isequal(binN1, binN2));
assert(isequal(valsOut1, valsOut2));
assert(isequal(idx1, idx2));
fprintf('Outputs match: [C version is %.2f times faster]\n', t1/t2);

