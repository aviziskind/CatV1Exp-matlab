% redoStats
global showWorking;
if isempty(showWorking)
    showWorking = 'oriSpf_rep';
end

Gid = redoStats_id(1);
cellId = redoStats_id(2);

varname = getName('celldata', Gid, cellId);
S = eval(varname);
if isfield(S, 'OSP') && isfield(S.OSP, 'stats')
    OSP = S.OSP;
    bckgRate = [];
    if isfield(S, 'PSTH') && isfield(S.PSTH, 'bckgRate') 
        bckgRate = S.PSTH.bckgRate;            
    end
    [Gid, R_full, R, oris] = deal(S.Gid, OSP.R_full, OSP.R, OSP.ori);
    f = 1;
    if isfield(OSP, 'R_full_factor')
        f = OSP.R_full_factor;
    end
    
    newstats = calcStatsFromOSPfull(R, R_full, f, Gid, bckgRate);
    OSP.stats = newstats;        
    eval([varname '.OSP = OSP;';]);       
end


% end