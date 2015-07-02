function [opt_gbparams, s] = fc_DirXYPsearch(xs, ys, phs, Gb_params, Gt_params_mtx, nspk, nfrm, lagnmb, nbin, mode, method)
% method = 'ro_rect', 'muinf', 'kld'

opt_gbparams = [];
s = [];

% nbin = 10;

% %==================================
% ori_deg = Gt_params_mtx(:,1);
% sp_pix = Gt_params_mtx(:,2);
% ph_deg = Gt_params_mtx(:,3);
%
% s = zeros(length(ys),length(xs),length(phs));
%
% for k = 1:( nfrm - lagnmb + 1)
%   r = nspk(k+lagnmb-1);
%   if r == 0
%     continue;
%   end
%   Sin_params = [0, 0, ori_deg(k), sp_pix(k), ph_deg(k)];
%   prj = func_GaborProj(Gb_params, Sin_params, xs, ys, phs);
%   sz = size(prj);
%   prj = prj(:);
%   prj(find(prj < 0)) = 0;
%   prj = reshape(prj,sz);
%   s = s + r*prj;
%   edbDispProgress(k,  nfrm, 100);
% end
%
% [mm, idx] = max(s(:));
% [ind_y, ind_x, ind_ph] = ind2sub(size(s),idx);
%
% x0 = xs(ind_x);
% y0 = ys(ind_y);
% ph0 = phs(ind_ph);
%
% opt_gbparams = [x0, y0, Gb_params(3), Gb_params(4), Gb_params(5), Gb_params(6), ph0];
% %==================================



ori_deg = Gt_params_mtx(:,1);
sp_pix = Gt_params_mtx(:,2);
ph_deg = Gt_params_mtx(:,3);

s = zeros(length(ys),length(xs),length(phs));

prjs = zeros(( nfrm - lagnmb + 1), length(ys),length(xs),length(phs));
rs = zeros(nfrm - lagnmb + 1,1);
for k = 1:( nfrm - lagnmb + 1)
    r = nspk(k+lagnmb-1);
    if r == 0
        continue;
    end
    Sin_params = [0, 0, ori_deg(k), sp_pix(k), ph_deg(k)];
    
    [prj_sin, prj_cos] = fc_GaborProj_SC(Gb_params, Sin_params, xs, ys, phs);
    %   prj = func_GaborProj(Gb_params, Sin_params, xs, ys, phs);
    if mode == 0  % if simple cell
        prjs(k,:,:,:) = prj_sin;
    else
        prjs(k,:,:,:) = prj_sin; % prj_cos; % sqrt(prj_sin.^2 + prj_cos.^2);
    end
    rs(k) = r;
    edbDispProgress(k,  nfrm, 100);
end

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
switch method
    case {'PPR_ro'}  % -------------------------------------------------------
        
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    p_s(find(p_s < 0)) = 0;   % --- rectification ---
                    %           [p_s, rs] = fc_Median3_ProjNspk(p_s, rs);
                    
                    %           ro = sum(p_s.*rs);         % --- crossproduct ---
                    if std(p_s) > 0 & std(rs) > 0
                        ro_mtx = corrcoef(p_s, rs);
                        ro = ro_mtx(1,2);
                    else
                        ro = 0;
                    end
                    s(ky,kx,kp) = ro;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        figure; plot(p_s, rs, '.k');
    case {'PPR_pos_poly'}  % -------------------------------------------------------
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    inds = find(p_s >= 0);
                    ps_p = p_s(inds);
                    rs_p = rs(inds);
                    S_tot = sum(rs_p.^2);
                    [P,S,MU] = polyfit(ps_p,rs_p,3);
                    [rs_hat, Delta] = POLYVAL(P,ps_p,S,MU);
                    e_res = rs_p - rs_hat;
                    S_res = sum(e_res.^2);
                    MF = 1 - S_res/S_tot;
                    s(ky,kx,kp) = MF;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        figure; plot(p_s, rs, '.k');
    case {'PPR_posneg_poly'}  % -------------------------------------------------------
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    inds_p = find(p_s >= 0);
                    inds_n = find(p_s < 0);
                    p_p = p_s(inds_p);
                    p_n = p_s(inds_n);
                    rs_p = rs(inds_p);
                    rs_n = rs(inds_n);
                    
                    S_tot_p = sum(rs_p.^2);
                    S_tot_n = sum(rs_n.^2);
                    
                    [P,S,MU] = polyfit(p_p,rs_p,3);
                    [rs_hat, Delta] = POLYVAL(P,p_p,S,MU);
                    e_res_p = rs_p - rs_hat;
                    S_res_p = sum(e_res_p.^2);
                    
                    [P,S,MU] = polyfit(p_n,rs_n,3);
                    [rs_hat, Delta] = POLYVAL(P,p_n,S,MU);
                    e_res_n = rs_n - rs_hat;
                    S_res_n = sum(e_res_n.^2);
                    
                    S_tot = S_tot_p + S_tot_n;
                    S_res = S_res_p + S_res_n;
                    
                    MF = 1 - S_res/S_tot;
                    s(ky,kx,kp) = MF;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        figure; plot(p_s, rs, '.k');
    case {'PPR_poly'}  % -------------------------------------------------------
        S_tot = sum(rs.^2);
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    [P,S,MU] = polyfit(p_s,rs,4);
                    [rs_hat, Delta] = POLYVAL(P,p_s,S,MU);
                    e_res = rs - rs_hat;
                    S_res = sum(e_res.^2);
                    MF = 1 - S_res/S_tot;
                    s(ky,kx,kp) = MF;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        figure; plot(p_s, rs, '.k');
    case {'PPR_FS'}  % -------------------------------------------------------
        
        for kx = 1:length(xs)
            %       x = xs(kx);
            for ky = 1:length(ys)
                %         y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    merit_value = fc_FriedmanStuedzleMerit(p_s, rs);
                    s(ky,kx,kp) = merit_value;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        
    case {'PPR_muinf'}  % -------------------------------------------------------
        
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    %           [p_s, rs] = fc_Median3_ProjNspk(p_s, rs);
                    
                    h_p12 = jointDistr([p_s, rs], nbin);
                    h_p1 = squeeze(sum(h_p12,2));
                    h_p2 = squeeze(sum(h_p12,1));
                    H1 = edbGetEntropy(h_p1);
                    H2 = edbGetEntropy(h_p2);
                    H12 = edbGetEntropy(h_p12);
                    
                    mui = H1 + H2 - H12;
                    
                    s(ky,kx,kp) = mui;
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
        
    case {'PPR_kld'}  % -------------------------------------------------------
        kld_max = -1;
        for kx = 1:length(xs)
            x = xs(kx);
            for ky = 1:length(ys)
                y = ys(ky);
                for kp = 1:length(phs)
                    p_s = prjs(:,ky,kx,kp);
                    kld = func_GetProjKLdist(p_s, rs, nbin);
                    s(ky,kx,kp) = kld;
                    if kld > kld_max
                        kld_max = kld;
                        ind_y = ky;
                        ind_x = kx;
                        ind_ph = kp;
                    end
                end
            end
            edbDispProgress(kx,  length(xs), 1);
        end
    otherwise  % -------------------------------------------------------
        return;
end
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

[mm, idx] = max(s(:));
[ind_y, ind_x, ind_ph] = ind2sub(size(s),idx);
p_s = prjs(:,ind_y, ind_x, ind_ph);
h_p12 = jointDistr([p_s, rs], nbin);
figure; imagesc(log10(h_p12 + eps));

x0 = xs(ind_x);
y0 = ys(ind_y);
ph0 = phs(ind_ph);

figure; plot(prjs(:,ind_y, ind_x, ind_ph), rs, '.k');

opt_gbparams = [x0, y0, Gb_params(3), Gb_params(4), Gb_params(5), Gb_params(6), ph0];



% if nargin < 9
%   mode = 'kld';
% end
%
% ori_deg = Gt_params_mtx(:,1);
% sp_pix = Gt_params_mtx(:,2);
% ph_deg = Gt_params_mtx(:,3);
%
% s = zeros(length(ys),length(xs),length(phs));
%
% prjs = zeros(( nfrm - lagnmb + 1), length(ys),length(xs),length(phs));
% rs = zeros(nfrm - lagnmb + 1,1);
% for k = 1:( nfrm - lagnmb + 1)
%   r = nspk(k+lagnmb-1);
%   if r == 0
%     continue;
%   end
%   Sin_params = [0, 0, ori_deg(k), sp_pix(k), ph_deg(k)];
%   prj = func_GaborProj(Gb_params, Sin_params, xs, ys, phs);
%   prjs(k,:,:,:) = prj;
%   rs(k) = r;
%   edbDispProgress(k,  nfrm, 100);
% end
%
% nbin = 10;
% % nbin = 20;
%
% switch mode
%
% case {'minfo'}
% 	for kx = 1:length(xs)
%     x = xs(kx);
%     for ky = 1:length(ys)
%       y = ys(ky);
%       for kp = 1:length(phs)
%         p_s = prjs(:,ky,kx,kp);
%         h_p12 = jointDistr([p_s, rs], nbin);
%         h_p1 = squeeze(sum(h_p12,2));
%         h_p2 = squeeze(sum(h_p12,1));
%         H1 = edbGetEntropy(h_p1);
%         H2 = edbGetEntropy(h_p2);
%         H12 = edbGetEntropy(h_p12);
%
%         mi = H1 + H2 - H12;
%
%         s(ky,kx,kp) = mi;
%       end
%     end
%     edbDispProgress(kx,  length(xs), 1);
% 	end
%
% 	[mm, idx] = max(s(:));
% 	[ind_y, ind_x, ind_ph] = ind2sub(size(s),idx);
% 	p_s = prjs(:,ind_y, ind_x, ind_ph);
% 	h_p12 = jointDistr([p_s, rs], nbin);
% 	figure; imagesc(log10(h_p12 + eps));
%
% 	case {'kld'}
% 		kld_max = -1;
% 		for kx = 1:length(xs)
%       x = xs(kx);
%       for ky = 1:length(ys)
%         y = ys(ky);
%         for kp = 1:length(phs)
%           p_s = prjs(:,ky,kx,kp);
%           kld = func_GetProjKLdist(p_s, rs, nbin);
%           s(ky,kx,kp) = kld;
%           if kld > kld_max
%             kld_max = kld;
%             ind_y = ky;
%             ind_x = kx;
%             ind_ph = kp;
%           end
%         end
%       end
%       edbDispProgress(kx,  length(xs), 1);
% 		end
% 	end
%
% x0 = xs(ind_x);
% y0 = ys(ind_y);
% ph0 = phs(ind_ph);
%
% figure; plot(prjs(:,ind_y, ind_x, ind_ph), rs, '.k');
%
% opt_gbparams = [x0, y0, Gb_params(3), Gb_params(4), Gb_params(5), Gb_params(6), ph0];

