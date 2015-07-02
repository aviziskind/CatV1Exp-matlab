function [gbparams, s] = fc_FIT_STRF(hnd, Gid, Cellnmb, params_ini, xs, ys, phs, lagnmb, mode)
gbparams = [];
s = [];

if nargin < 11
  mode = 'lin';
end

nbin = 10;

nspk_mtx = edbGetMovieNspkPerFrame(hnd, Gid, Cellnmb);

hdr = nspk_mtx(1,:);
nspk_mtx = nspk_mtx(2:end,:);

inds_200 = find(hdr == 200);
inds_201 = find(hdr == 201);
inds_213 = find(hdr == 213);
inds_214 = find(hdr == 214);

mtx_all = [];

if ~isempty(inds_200)
	nspk_200 = squeeze(sum(nspk_mtx(:,inds_200),2));
	nspk_201 = squeeze(sum(nspk_mtx(:,inds_201),2));
	
	Nfrm = size(nspk_200,1);
	
	[ori_deg_200, sp_pix_200, ph_deg_200] = func_GetIndexMtx(200);
	[ori_deg_201, sp_pix_201, ph_deg_201] = func_GetIndexMtx(201);
	
	mtx_200 = [ori_deg_200, sp_pix_200, ph_deg_200];
	mtx_201 = [ori_deg_201, sp_pix_201, ph_deg_201];
	
	mtx_all = [mtx_200; mtx_201];
	nspk_all = [nspk_200; nspk_201];
	Nfrm = 2*Nfrm;
else
  if ~isempty(inds_213)
		nspk_213 = squeeze(sum(nspk_mtx(:,inds_213),2));
		nspk_214 = squeeze(sum(nspk_mtx(:,inds_214),2));
		
		Nfrm = size(nspk_213,1);
		
		[ori_deg_213, sp_pix_213, ph_deg_213] = func_GetIndexMtx(213);
		[ori_deg_214, sp_pix_214, ph_deg_214] = func_GetIndexMtx(214);
		
		mtx_213 = [ori_deg_213, sp_pix_213, ph_deg_213];
		mtx_214 = [ori_deg_214, sp_pix_214, ph_deg_214];
		
		mtx_all = [mtx_213; mtx_214];
		nspk_all = [nspk_213; nspk_214];
		Nfrm = 2*Nfrm;
  else
    mid = unique(hdr);
% 		[ori_deg, sp_pix, ph_deg] = func_GetIndexMtx(mid);
        [uori, usp, uph, phase_deg, sp_pix, ori_deg] = GetOSPidxs(mid);
		mtx_all = [ori_deg, sp_pix, ph_deg];
		nspk_all = nspk_mtx;                                                                                                                
		Nfrm = size(nspk_all,1);
  end
end

if isempty(mtx_all)
  disp('Wrong Mid');
  return;
end

switch mode
  case {'Q'}
    [gbparams, s] = fc_DirXYPsearch(xs, ys, phs, params_ini, mtx_all, nspk_all, Nfrm, lagnmb, nbin, 0, 'XXX'); 
  case {'C'}
    [gbparams, s] = fc_DirXYPsearch(xs, ys, phs, params_ini, mtx_all, nspk_all, Nfrm, lagnmb, nbin, 0, 'PPR_ro');
  case {'A'}
    [gbparams, s] = fc_DirXYPsearch(xs, ys, phs, params_ini, mtx_all, nspk_all, Nfrm, lagnmb, nbin, 0, 'PPR_muinf');
end
