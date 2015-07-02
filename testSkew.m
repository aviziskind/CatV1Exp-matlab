a = exp(randn(1000,1)*5+3);
hist(a);


vr  = @(mu, sig) (exp(sig.^2)-1)*exp(2*mu+sig.^2);
skw = @(sig) (exp(sig.^2)+2)*(sqrt(exp(sig.^2)-1));
kur = @(sig) exp(4*sig.^2) + 2*exp(3*sig.^2) + 3*exp(2*sig.^2) - 6;


%want to find mu, sigma such that var = V, skewness = S.
V = 2;
S = 3;

sigma = fzero(@(s) skw(s)-S, 1);
mu = fzero(@(mu) vr(mu, sigma)-V, 1);
X = exp(randn(B,1)*sigma + mu);
[v,s] = deal(var(X), skewness(X))
s = skewness(X)


