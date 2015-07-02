function [ori_deg, sp_pix, ph_deg] = func_GetIndexMtx(Mid);
ori_deg = [];
sp_pix = [];
ph_deg = [];

fpath = '...\Input\'; % change it for a real directory
ftitle = ['mid',num2str(Mid),'.txt'];
fname = [fpath, ftitle];
mtx = load(fname);
ori_deg = mtx(:,1);
sp_pix = mtx(:,2);
ph_deg = mtx(:,3);
