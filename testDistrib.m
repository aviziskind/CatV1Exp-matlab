% function testDistrib
    fmt = '%.3f';
    str_mn_std = @(x) [num2str(nanmean(x), fmt) ' \pm ' num2str(nanstd(x)/sqrt(nnz(~isnan(x))), fmt)];
    Nt = 1e4;
    
    n= 8;
    cc = zeros(1, Nt);
    progressBar('init', Nt, 20);
    for i = 1:Nt
        progressBar(i);
        x = zeros(1,n);
        y = zeros(1,n);
        x(randi(n-1,1,1)+1) = 1;
        y(randi(n-1,1,1)+1) = 1;
        cc(i) = pearsonR(x, y);
    end
    progressBar('done');
    
    figure(1);
    hist(cc);
    title(str_mn_std(cc));
    
    
% end