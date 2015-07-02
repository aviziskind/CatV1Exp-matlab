Gid = 4470;
cellId = 2;

hnd = dbOpenExpDb;
Did  = edbFindDID(hnd,Gid);


% 	x0_gb = gb_params(1);
% 	y0_gb = gb_params(2);
% 	Ori_deg_gb = gb_params(3);
%     Sp_pix_gb = gb_params(4);
% 	Sx2 = gb_params(5);
% 	Sy2 = gb_params(6);
% 	ph0_deg_gb = gb_params(7);
    
	x0_gb = 0;
	y0_gb = 1;
	Ori_deg_gb = 30;
    Sp_pix_gb = 10;
	Sx2 = 1;
	Sy2 = 1.2;
	ph0_deg_gb = 90;
    
    Mids = dbLookup('Mid',  'Did', Did);
    Mid1 = Mids(1);
    
    [uori, usp, uph, phase_deg, sp_pix, ori_deg] = GetOSPidxs(Mid1);
    
    
    
    Gb_params = [x0_gb; y0_gb; Ori_deg_gb; Sp_pix_gb;  Sx2;  Sy2;  ph0_deg_gb];
    
    
    
    lagnmb = 1;
    mode = 'A';

    xs = 1:10;
    ys = 1:10;
    phs = 0:30:180-30;
    
[gbparams, s] = fc_FIT_STRF(hnd, Gid, cellId, Gb_params, xs, ys, phs, lagnmb, mode);


