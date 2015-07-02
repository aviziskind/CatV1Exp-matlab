function testQuadProdGaussians3

%%
D = 1;
N = 1000;

s = 1;


% X2 = X2_orig;

if D == 1
    m1 = 0;       
    m2 = m1;
    s1 = 1;      
    s2 = s1*1.5;

elseif D == 2
end


all_s = [0:.5:8]; nS = length(all_s);


dists_norm = nan(1, nS);


showGaussians = 1;

if D == 1
    %%
    figure(200); clf; 
    subplot(11,1,1:5);
    hold on; box on;    
    h_ax = gca;
    h_gs = plot([0 0], zeros(2,3), '-');
    xlabel('');

    subplot(11,1,7:11);
    hold on; box on;
    h_dist = plot(0,0);
    h_curdist = plot(0,0);
    h_ax_dist = gca;
    
    dists_norm_max = -quadProdGaussians(m1, s1, m1+all_s(end), s2, 'log', [], [], 1);
    set(h_ax_dist, 'xlim', lims(all_s), 'ylim', [0, dists_norm_max]);    
    ylabel('-log_{10} overlap');
    xlabel('Distance between Centers');
    
    
    


elseif D == 2
    %%
    figure(201); clf; hold on; box on;    
    h_ax = gca;
    for i = 1:3
        h_surf(i) = surf([0, 0], [0, 0], zeros(2));                
    end
end


    
%%

    %%
    for si = 17;% 1:nS; %1:nS
        %%
        
        if D == 1
            m1 = 0;       
            m2 = m1 + all_s(si);
            s1 = 1;      
            s2 = s1*1.5;
            
        elseif D == 2
            M1 = mean(x1,2);
            M2 = mean(x2,2);
            C1 = cov(x1');
            C2 = cov(x2');
        end
                        
        
        if showGaussians 
            if D == 1
                %%
            
                L = [m1-4*s1, m1+4*s2 + all_s(end)];
                % L = [-.75, 1.75];
                X = linspace(L(1), L(2), 200);
                                                    
                f = 1;
                gs1 = gaussian(X, m1, s1); gs1 = gs1/f;
                gs2 = gaussian(X, m2, s2); gs2 = gs2/f;
                gs_prod = gs1.*gs2;  
                
                gs_prod = gs_prod/max(gs_prod)*max(gs1)*max(gs2)*1*(1/(all_s(si)+all_s(2)));
                col1 = 'b'; col2 = [0 .5 0]; col3 = 'r';

                set(h_gs(3), 'xdata', X, 'ydata', gs1, 'color', col1, 'linewidth', 3, 'linestyle', '-')
                set(h_gs(2), 'xdata', X, 'ydata', gs2, 'color', col2, 'linewidth', 3, 'linestyle', '-');
                set(h_gs(1), 'xdata', X, 'ydata', gs_prod, 'color', col3, 'linewidth', 3, 'linestyle', ':');

                dists_norm(si) = -quadProdGaussians(m1, s1, m2, s2, 'log', [], [], 1);
                
                set(h_dist, 'xdata', all_s, 'ydata', dists_norm, 'marker', '.');
                set(h_curdist, 'xdata', all_s(si), 'ydata', dists_norm(si), 'marker', 'o', 'color', 'r', 'linewidth', 2);
            3;
            
            elseif D == 2
                
              
                %%
                figure(201); clf; hold on; box on;
                h_ax = gca;
                n = 100;
                for i = 1:3
                    h_surf(i) = surf(1:n, 1:n, zeros(n));
                end
                view(3);
            %%
                L1 = lims([x1(1,:), x2(1,:) + all_s(end)], .1);
                L2 = lims([x1(2,:), x2(2,:)], .1);
                % L = [-.75, 1.75];
                
                X = linspace(L1(1), L1(2), 50);
                Y = linspace(L2(1), L2(2), 40);
                [X_grid, Y_grid] = meshgrid(X, Y);                                

    %         m1 = mean(pca_1_proj); s1 = std(pca_1_proj);
    %         m2 = mean(pca_2_proj); s2 = std(pca_2_proj);
                ii = 1:2;
                m1 = M1(ii); m2 = M2(ii); c1 = C1(ii, ii); c2 = C2(ii,ii);
%                 m_c = m1+m2
    
                gs1_v = gaussianN([X_grid(:), Y_grid(:)], m1, c1);
                gs1 = reshape(gs1_v, length(Y), length(X));

                gs2_v = gaussianN([X_grid(:), Y_grid(:)], m2, c2);
                gs2 = reshape(gs2_v, length(Y), length(X));
                
                gs_prod = gs1.*gs2*3;
                col1 = 'b'; col2 = [0 .5 0]; col3 = 'r';

                c1 = ones(size(gs1));
                rescale = @(x, mx, mn) (x-mn)/(mx-mn);
                mn_val = -.05;
                
                gs1_r = rescale(gs1, max(gs1(:)), mn_val);
                gs2_r = rescale(gs2, max(gs2(:)), mn_val);
                gs_prod_r = rescale(gs_prod, max(gs_prod(:)), mn_val);
                
                
                col1 = cat(3,    mn_val*c1,  mn_val*c1, gs1_r .*  c1);
                col2 = cat(3,    mn_val*c1, gs2_r.*c1, mn_val*c1);                
                col3 = cat(3,    gs_prod_r.*c1, mn_val*c1, mn_val*c1);
                
                y_use = find(Y > -.5);
%                 y_use = true(size(Y));
                
                set(h_surf(1), 'xdata', X, 'ydata', Y(y_use), 'zdata', gs1(y_use,:), 'cdata', col1(y_use,:,:), 'facealpha', .5); 
                set(h_surf(2), 'xdata', X, 'ydata', Y(y_use), 'zdata', gs2(y_use,:), 'cdata', col2(y_use,:,:), 'facealpha', .5); 
                set(h_surf(3), 'xdata', X, 'ydata', Y(y_use), 'zdata', gs_prod(y_use,:), 'cdata', col3(y_use,:,:), 'facealpha', .5); 
                set(h_ax, 'box', 'off', 'xlim', lims(X), 'ylim', lims(Y(y_use)))
%                 axis tight;
%                 xlabel('x');
%                 ylabel('y');
                xlabel(''); ylabel('');
                set(gcf,'color', 'w');
                set(h_ax, 'ztick', []);
%                 ylim auto;
                
                %%
%                 figure(444);
%                 surf(x1, x2, gs1);
%                 set(h_gs(2), 'xdata', binC_fine, 'ydata', gs2, 'color', col2, 'linewidth', 3);
%                 set(h_gs(3), 'xdata', binC_fine, 'ydata', gs_prod, 'color', col3, 'linewidth', 3, 'linestyle', ':');

            3;
               

            end
        end
        
        drawnow;
        pause(.3);
        
    end

end            
                        


