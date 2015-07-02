
nData = 2000000;
nBins = 10;
nCats = 4;
tf = double(rand(nData, 1) < .8);
pairIdxs = unique(randi(nData, round(nData/10), 1))';
allVals = randn(nData, 1);

binIds = randi(nBins, nData,1)';
catIds = randi(nCats, nData,1)';
% tf;
tic;
[N1, vals1] = binCountForPairs(  pairIdxs, tf, binIds, [], nBins, nCats, allVals);
t1 = toc;
toc;
3;

tic;
[N2, vals2] = binCountForPairs_c(pairIdxs, tf, binIds, [], nBins, nCats, allVals);
t2 = toc;
toc;

t1/t2
ok1 = isequal(N1, N2);
ok2 = isequal(vals1, vals2);
if ~ok1 || ~ok2
    error('Wrong!');
end
