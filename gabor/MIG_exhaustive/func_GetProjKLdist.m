function kld = func_GetProjKLdist(prj, nspk, Nbins);
kld = [];

ma_prj = unique(max(prj));
mi_prj = unique(min(prj));
d_prj = ( ma_prj - mi_prj )/(Nbins - 1);

Nprj = length(prj);

edges_prj = (mi_prj:d_prj:ma_prj);
h = histc(prj,edges_prj);
h1 = zeros(size(h));

inds = 1 + fix((prj - mi_prj)/d_prj);

for k = 1:length(h)
  inds_k = find(inds == k);
  h1(k) = sum(nspk(inds_k));
end

h1 = h1./sum(h1);
h = h./sum(h);

inds = find(h1 > 0 & h > 0);
h1 = h1(inds);
h = h(inds);

kld = sum(h1.*log2(h1./h));
