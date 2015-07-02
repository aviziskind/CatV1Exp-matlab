% S = load('allWindowOSPs_mean')
tic;
fn = fieldnames(S);
for i = 1:length(fn)
    v = S.(fn{i});
    n = nnz(v.l_bins > 0);
%     n_set = n;
    n_set = max(n, 10);
    v.l_bins(n_set+1:end) = [];
    v.r_bins(n_set+1:end) = [];
    v.osps(n_set+1:end) = [];
    v.osps_odd(n_set+1:end) = [];
    v.osps_even(n_set+1:end) = [];    
    
    S.(fn{i}) = v;
end
toc;