function [ret, tb, te] = dbMockMakeTbTe(hnd, Did, Glitch_check_flag)
% [ret, tb, te] = edbMakeTbTe(hnd, Did);
% WARNING:
% This function is designed to be called to write Tb Te into the database when syncs are read first time from a disk file
% or it may be used for TbTe editing (if some discrepancy in EDB is detected) 
% IT IS NOT INTENDED TO BE USED FOR OTHER PURPOSES
% 

ret = [];
tb = [];
te = [];

if nargin < 2
  return;
end

if nargin == 2
  Glitch_check_flag = 1;
end

if isempty(Did) == 1
   return;
end

% ================================ PREPARATIONS =============================================

% --- get synch group id -----------------------

syncs = dbGetSyncs('Did',Did, 'tick');

stimTableName = getDatabaseTableForDid(Did);

% --- get Sampling RATE and monitor RATE ----------------------

ticksPerSec = dbConvertTimeMeasures(Did, 1, 'sec', 'tick');

if checkDetails
    ticksPerSecCheck = unique( getFieldsFromDatabaseTable(hnd,'DBL_SAMPLING_RATE_HZ', 'TBL_DATA_FILES', {'DATAFILE_ID', Did}), [], 1);
    if length(ticksPerSecCheck) > 1
        error('Inconsistent sampling rate across presentations')
    end
    if (ticksPerSecCheck ~= ticksPerSec)
        error('Incorrect sampling rate');
    end
end

framesPerSec = 1/dbConvertTimeMeasures(Did, 1, 'sec', 'frame');
if checkDetails
    if length(framesPerSecCheck) > 1
        error('Inconsistent frame rate across presentations')
    end
    if (framesPerSecCheck ~= framesPerSec)
        error('Incorrect frame rate');
    end
end

% --- GET frame duration in ticks (approximate value)
ticksPerFrame = round(ticksPerSec/framesPerSec);

% --- GET list of presentations ---
Pnmb = getFieldsFromDatabaseTable(hnd, 'LNG_PRESENT_NO', stimTableName, {'DATAFILE_ID', Did});
if isempty(Pnmb) == 1
    return;
end
Npres_ord = length(Pnmb);

% LNG_N_SUSTAINED_FRM      = [];
% LNG_N_PRE_BLANK_FRM      = [];
% LNG_N_POST_BLANK_FRM     = [];
% LNG_N_FADEIN_FRM         = [];
% LNG_N_FADEOUT_FRM        = [];
% 
% LNG_N_DISPLAYED_FRM      = [];
% LNG_N_MISSED_FRM         = [];
% LNG_LAST_MISSED_FRM_NO   = [];
% 
% for k = 1:Npres_ord
%     tmp_info = edbGetInfo(hnd, tbl_name, 'LNG_PRESENT_NO', Pnmb(k));
%     if isempty(tmp_info) == 1
%         return;
%     end
%     LNG_N_SUSTAINED_FRM      = [LNG_N_SUSTAINED_FRM;edbParseInfo(tmp_info,'LNG_N_SUSTAINED_FRM')];
%     LNG_N_PRE_BLANK_FRM      = [LNG_N_PRE_BLANK_FRM;edbParseInfo(tmp_info,'LNG_N_PRE_BLANK_FRM')];
%     LNG_N_POST_BLANK_FRM     = [LNG_N_POST_BLANK_FRM;edbParseInfo(tmp_info,'LNG_N_POST_BLANK_FRM')];
%     LNG_N_FADEIN_FRM         = [LNG_N_FADEIN_FRM;edbParseInfo(tmp_info,'LNG_N_FADEIN_FRM')];
%     LNG_N_FADEOUT_FRM        = [LNG_N_FADEOUT_FRM;edbParseInfo(tmp_info,'LNG_N_FADEOUT_FRM')];
%     
%     LNG_N_DISPLAYED_FRM      = [LNG_N_DISPLAYED_FRM;edbParseInfo(pres_info,'LNG_N_DISPLAYED_FRM')];
%     LNG_N_MISSED_FRM         = [LNG_N_MISSED_FRM;edbParseInfo(pres_info,'LNG_N_MISSED_FRM')];
%     LNG_LAST_MISSED_FRM_NO   = [LNG_LAST_MISSED_FRM_NO;edbParseInfo(pres_info,'LNG_LAST_MISSED_FRM_NO')];
% end

% --- odered parameters ---
nSustained

LNG_N_SUSTAINED_FRM      = edbParseInfo(pres_info,'LNG_N_SUSTAINED_FRM');
LNG_N_PRE_BLANK_FRM      = edbParseInfo(pres_info,'LNG_N_PRE_BLANK_FRM');
LNG_N_POST_BLANK_FRM     = edbParseInfo(pres_info,'LNG_N_POST_BLANK_FRM');

% --- NOTE: LNG_FRAMES_PER_UPDATE exists in all tables except tbl_grating_pres only!
LNG_FRAMES_PER_UPDATE    = edbParseInfo(pres_info,'LNG_FRAMES_PER_UPDATE');
if isempty(LNG_FRAMES_PER_UPDATE) == 1
    LNG_FRAMES_PER_UPDATE = ones(length(LNG_N_SUSTAINED_FRM),1);
end

% --- NOTE: FADEIN/FADEOUT exist in tbl_grating_pres only!
LNG_N_FADEIN_FRM         = edbParseInfo(pres_info,'LNG_N_FADEIN_FRM');
if isempty(LNG_N_FADEIN_FRM) == 1
    LNG_N_FADEIN_FRM = zeros(length(LNG_N_SUSTAINED_FRM),1);
end
LNG_N_FADEOUT_FRM        = edbParseInfo(pres_info,'LNG_N_FADEOUT_FRM');
if isempty(LNG_N_FADEOUT_FRM) == 1
    LNG_N_FADEOUT_FRM = zeros(length(LNG_N_SUSTAINED_FRM),1);
end
% --- END OF EXCEPTION

% --- find the minimal pause and expected presentation time in frames
Minpause_frm = LNG_N_PRE_BLANK_FRM + LNG_N_POST_BLANK_FRM; 
Ptime_frm = LNG_N_FADEIN_FRM + LNG_N_SUSTAINED_FRM + LNG_N_FADEOUT_FRM;

if isempty(Minpause_frm) == 1 | isempty(Ptime_frm)
    return;
end

% --- returned parameters ---
LNG_N_DISPLAYED_FRM      = edbParseInfo(pres_info,'LNG_N_DISPLAYED_FRM');
% LNG_N_MISSED_FRM         = edbParseInfo(pres_info,'LNG_N_MISSED_FRM');
% LNG_LAST_MISSED_FRM_NO   = edbParseInfo(pres_info,'LNG_LAST_MISSED_FRM_NO');


% --- database information consistency check ---

Nfrm_ord = sum(Ptime_frm + Minpause_frm); % CHECK IF PRE- POST- BLANKS ARE REPORTED IN A FRAMES FILE
tmp = (Ptime_frm + Minpause_frm) - LNG_N_DISPLAYED_FRM; % + LNG_N_MISSED_FRM);

% --- IF IT IS NECESSARY??? ----
 Ndiscrep = length(find(tmp));
if Ndiscrep ~= 0
    disp('ERROR: inconsistent presentation parameters');
    return;
end
% ------------------------------ 
 
fpu = unique(LNG_FRAMES_PER_UPDATE);
if length(fpu) ~= 1
    disp('ERROR: Frames_per_update is not unique across presentations');
    return;
end

% --- Array of frame times ---
t = syncs(:,1);
% --- Array of between-frame intervals ---
dt = diff(t);

% --- GAP between successive presentations is 1.5 times greater than a regular interframe distance ---
gap_tick = 2.5 * T_tick * fpu;
%gap_tick = 2000;    % REMOVE IT!
% --- make indices for beg and end of presentations ---
inds = find(dt > gap_tick);
ind_b = [1;inds + 1];
ind_e = [inds;length(t)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tb = t(ind_b);
% te = t(ind_e);
% ind_skip = find((te-tb) == 275485)
% tb = [tb(1:ind_skip-1);tb(ind_skip+1:end)]    
% te = [te(1:ind_skip-1);te(ind_skip+1:end)]    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Number of presentations found ---
Npres_found = length(ind_b);

if Npres_found ~= Npres_ord
    % b = [t(ind_b) t(ind_e)]
%     figure;
%     plot([max(t(ind_e)); max(t(ind_e))],[0; 2],'-w');
%     hold on;
%     for ppp = 1:Npres_found
%         plot([t(ind_b(ppp));t(ind_b(ppp))],[0;0.6],'b-');
%         hold on;
%         plot([t(ind_e(ppp));t(ind_e(ppp))],[0.4;1],'r-');
%         hold on;
%     end
%     hold off;
% for ppp = 1:Npres_found-1
%     t(ind_b(ppp+1)) - t(ind_e(ppp))
% end    
    disp('ERROR: inconsistent data: number of presentation found is different from what is reported by VIDEO SERVER');
    return;
end


if Glitch_check_flag 
  inds = find(dt <= gap_tick);
  nfrms_expected = length(inds)/Npres_ord;        
  if fix(nfrms_expected) ~= fix(nfrms_expected)
%     length(inds)
%     Npres_ord
%     nfrms_expected = length(inds)/Npres_ord
    disp('ERROR: inconsistent data. Inseparable inter-presentation glitches or within presentation glitches found');
    return;
  end
end

% --- assign tb & te ---
tb = t(ind_b);
te = t(ind_e);

% Ngaps = length(tb);
% for pnmb = 1:Ngaps
%   iii = find(t >= tb(pnmb) & t < te(pnmb));
%   length(iii)
% end
% return;

% =================================== DB UPDATING =====================================
% --- Now writing them into the database ---

tblname = edbGetPresTbl(STIMULUS_TYPE_ID);
if isempty(tblname)
   return;
end

[tbInDB, teInDB] = getFieldsFromDatabaseTable(hnd, {'LNG_START_TICK', 'LNG_END_TICK'}, tblname, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');

output = [{'BeginTick: Now', 'BeginTick: DB', 'diff', 'EndTick: Now', 'EndTick: DB', 'diff'}; ...
           num2cell([tb, tbInDB, tb-tbInDB, te, teInDB, te-teInDB])];
disp(output);



pres_id_name = edbGetPresIdName(STIMULUS_TYPE_ID);
if isempty(pres_id_name) == 1
   return;
end

% setFieldsInDatabaseTable(hnd, [tb, te], {'LNG_START_TICK', 'LNG_END_TICK'}, tblname, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO')

% === END OF FUNCTION =================================================================================

