% Use Sergei's code in this one, but don't divide the data into 2
tic;

Gid = 4470;
Cellnmb =2;

% % Gid = 4712;
% % Gid = 4726;
BkgrSkip_ms = 200;


[spk_t, frm_t, spk_i_pres, spk_i_frm, Mids, N_frm_pres, frame_len_ms, N_bkg, bkg_t_ms] = ...
                    GetParsedSpikesForMovie(  Gid, Cellnmb, BkgrSkip_ms );

bckg_rate = N_bkg/bkg_t_ms * 1000; 


rel_t = spk_t - frm_t;

iPres = 1;
% we assume that all movies have the same set of oris and spfs
for Mid = Mids'
  [uori, usp, uph, phase_deg, sp_pix, ori_deg] = GetOSPidxs(Mid);
  if( iPres == 1)
		N_ori = length(uori); 
		N_spf = length(usp);
		osp = zeros(N_ori, N_spf);
  end
  for iOri = 1:N_ori
    for iSpf = 1:N_spf
      ori = uori(iOri);
      spf = usp(iSpf);
      frm_idx = find((sp_pix==spf)&(ori_deg==ori));
      for iFrm = frm_idx'
          osp(iOri, iSpf) = osp(iOri, iSpf) + length(find((spk_i_pres == iPres) & (spk_i_frm == iFrm)));          
      end
    end
  end
  iPres = iPres+1;
end     

N_lag = 1;
blnCalcDiff = 0;
if blnCalcDiff
  N_grad = 511;
else
  N_grad = 256;  % # grey levels (?)
end

hnd = dbOpenExpDb;
% Data file ID
Did  = edbFindDID(hnd,Gid);
% Frame size
[nx, ny] = edbGetMovieDim(hnd, Did);
if isempty(nx)
    return;
end
n_x=nx(1);
n_y=ny(1);

% we assume that all movies had same dimensions
N_pix = n_x * n_y; 

tot_hist = zeros( N_grad, N_pix, N_lag );
tot_ref_hist = zeros( N_grad, N_pix);

movie_hist = zeros( N_grad, N_pix, N_lag);
movie_ref_hist = zeros( N_grad, N_pix);

iPres = 1;
for Mid = Mids'

	Movie_fname = dbGetMoviePathAndFilename(hnd, Mid);
	
	N_frm = N_frm_pres(iPres);
    N_spk_frm = [];
  % find spikes per frame for this presentation
  for iFrm = 1:N_frm;
      N_spk_frm = [N_spk_frm; length(find((spk_i_pres == iPres) & (spk_i_frm == iFrm)))]; %#ok<AGROW>
  end

    [movie_hist, movie_ref_hist] = PixelsHist( Movie_fname, N_pix, N_frm, N_spk_frm, N_lag, blnCalcDiff );
    % add the hist for this presentation to the total one                     
    tot_hist = tot_hist + movie_hist;
    tot_ref_hist = tot_ref_hist + movie_ref_hist;


  iPres = iPres+1;
end     



if blnCalcDiff
  aver_int = 256;
else
  aver_int = 127;
end

% kern = zeros( N_pix, N_lag );

if N_lag == 1
    kern = zeros( 1, N_pix );
else
    kern = zeros( N_pix, N_lag );
end    

for grad = 1: N_grad
    kern = kern + squeeze(tot_hist(grad, :, :))*(grad - aver_int);
end

kern = reshape(kern, n_y, n_x, N_lag);


% length(cell_t)/max(cell_t) * SamplingRateHz
N_all_frm = sum(N_frm_pres);
N_bins = 20;
% 
figure(1);

ha_0 = subplot(1,3,1);
plotPSTH( rel_t, N_bins, frame_len_ms, 1/N_all_frm, bckg_rate );
axis square;

% hi = max(  max(abs(osp))  );

subplot(1,3,2);
imagesc(uori, usp, osp); 
axis square;
hi = max(  max(abs(kern))  );

subplot(1,3,3);
imagesc(kern, [-hi hi]); 
axis square xy;

strSel = 'SELECT TBL_ANIMALS.TXT_LAB_NAME, TBL_PENETRATIONS.PENETRATION_ID, TBL_LOCATIONS.LOCATION_ID, TBL_DATA_FILES.LNG_DATAFILE_NO,  TBL_DATA_FILES.DBL_SAMPLING_RATE_HZ ';
strFrm = 'FROM TBL_ELECTRODE_TYPES INNER JOIN ((((TBL_ANIMALS INNER JOIN ((TBL_AP_ML_ZERO INNER JOIN TBL_ELECTRODES ON TBL_AP_ML_ZERO.ELECTRODE_ID = TBL_ELECTRODES.ELECTRODE_ID) INNER JOIN TBL_PENETRATIONS ON TBL_AP_ML_ZERO.AP_ML_ZERO_ID = TBL_PENETRATIONS.AP_ML_ZERO_ID) ON TBL_ANIMALS.ANIMAL_ID = TBL_ELECTRODES.ANIMAL_ID) INNER JOIN TBL_LOCATIONS ON TBL_PENETRATIONS.PENETRATION_ID = TBL_LOCATIONS.PENETRATION_ID) INNER JOIN (TBL_LOC_DEPTHS INNER JOIN (TBL_DATA_FILES INNER JOIN TBL_LOCS_FILES_LINKS ON TBL_DATA_FILES.DATAFILE_ID = TBL_LOCS_FILES_LINKS.DATAFILE_ID) ON TBL_LOC_DEPTHS.LOC_DEPTH_ID = TBL_LOCS_FILES_LINKS.LOC_DEPTH_ID) ON TBL_LOCATIONS.LOCATION_ID = TBL_LOC_DEPTHS.LOCATION_ID) INNER JOIN TBL_GROUPS ON TBL_DATA_FILES.DATAFILE_ID = TBL_GROUPS.DATAFILE_ID) ON TBL_ELECTRODE_TYPES.ELECTRODE_TYPE_ID = TBL_ELECTRODES.ELECTRODE_TYPE_ID ';
strWhr = ['WHERE (((TBL_ELECTRODES.ELECTRODE_TYPE_ID)=2) AND ((TBL_GROUPS.GROUP_ID)=' num2str(Gid) '));'];

data = edbRunSelQry(hnd, [strSel strFrm strWhr]);
edbClose(hnd);

CatName = char(cat(1,data{:,1}));
PenetrID = double(cat(1,data{:,2}));
LocID = double(cat(1,data{:,3}));
FileNmb = double(cat(1,data{:,4}));
SampRateHz = double(cat(1,data{:,5}));

str_1 = [ CatName '   PenID: ' num2str(PenetrID) '   LocID: ' num2str(LocID) '   Gid: ' num2str(Gid)  ];
str_2 = ['File #: ' num2str(FileNmb) '  Cell# ' num2str(Cellnmb) '  N spikes: ' num2str(length(rel_t)) ];

suptitle_2([str_1 sprintf('\n') str_2]);

toc;            
      
