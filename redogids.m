S_f = load('cellsGroups_movie_fg.mat');
fg_gids = [S_f.movieGroups_fg.Gid];
S_d = load('cellsGroups_grating_dSf.mat');
dg_gids = [S_d.gratingGroups_dSf.Gid];
Gids = [fg_gids(:); dg_gids(:)];

nGids = length(Gids);
n_ok = zeros(1, nGids);
progressBar('init-', nGids)
for i = 1:length(Gids)
    Gid = Gids(i);
    progressBar;
    groupSpikes = getSpikes(Gid); 
    n1 = size(groupSpikes,1);
    pca = getGroupWaveformCoefficients('PCA', Gid, 2, 'separate', 1);
    n2 = size(pca,1);
    n_ok(i) = (n1 == n2);
end
                    
                    
