% cleanupAndFixVariables

% 1. clear variables with Gid and cellID names that don't match inner values
% s = who('-regexp', [varnamePrefix '*']);
% for si = 1:length(s)
%     nm = s{si};
%     [Gid_nm, cellId_nm] = getGidCellIdFromVarname(nm);
% 
%     v  = eval(nm);
%     if isstruct(v)
%         Gid_st = v.Gid;
%         cellId_st = v.cellId;
%         if (Gid_nm ~= Gid_st) || (cellId_nm ~= cellId_st)
%             disp(['clearing ' nm ' (incorrect)']);
%             clear(nm);
%         end
%     end
%     
% end
    

% 2. clear variables with Gid and cellID names that don't match stimulus type
if exist('stimType', 'var')
    disp(['Removing all non-' stimType ' cells from memory before saving...']);
    s = who([varnamePrefix '_*']);
    if exist('stimType', 'var')
        S = load(['cellsGroups_' stimType '_all_DB']);
        groupData = S.([stimType 'Groups_all']);
        GidsThisStimType = [groupData.Gid];

        for si = 1:length(s)
            nm = s{si};
            [Gid_nm, cellId_nm] = getGidCellIdFromVarname(nm);
            if ~any(Gid_nm == GidsThisStimType)
                disp(['clearing ' nm ' (not a ' stimType ' cell)']);
                clear(nm);
            end

        end

    end
end


%     
% for si = 1:length(s)
%     nm = s{si};
%     v = eval(nm);
%     if iscell(v)
%         v = v{1};
%         eval([nm ' = v;']);
%         disp(['shortened ' nm]);
%     end
% end
%             
%     


% % replace 'old' variable names with 'new' variable names
% s = who('-regexp', [varnamePrefix '*']);
% for si = 1:length(s)
%     oldname = s{si};
%     v = eval(oldname);
%     [G, C] = getGidCellIdFromVarname(oldname);
%     newname = getName('celldata', G, C);
%     eval([newname ' =  v;']);
%     clear(oldname);
% end

% ** update STA data: from {delays, windowSizes} --> {timeWindows}
% tic;
% disp('Starting');
% s = who('-regexp', [varnamePrefix '*']);
% if exist('stimType', 'var')
%     S = load(['cellsGroups_' stimType '_DB']);
%     groupData = S.([stimType 'Groups']);
%     GidsThisStimType = [groupData.Gid];
%     
%     for si = 1:length(s)
%         nm = s{si};
%         vr = eval(nm);
%         if isstruct(vr)
%             STAsC = vr.STAs;
%             if length(STAsC) == 3
%                 [STAs, delays_ms, windowSizes_ms] = elements(STAsC);
%                 timeWindows = [delays_ms(:), delays_ms(:)+windowSizes_ms(:)];
%                 STAsNew = {STAs, timeWindows};
%                 vr.STAs = STAsNew;
%                 eval([nm ' = vr']);
%             end
%         end
%     
%     end
%     
% end
% toc;




% 1. Fix OSP indices - if in radians - change to degrees
% if strcmp(stimType, 'movie');
%     phsSet_rad = [0; 1.570796; 3.141593; 4.712389];
%     s = who('-regexp', [varnamePrefix '*']);
% %     groupData = 
% %     S = load(['cellsGroups_' stimType '_DB']);
% %     groupData = S.([stimType 'Groups']);
% %     flashGratingGroupInds = findInStructArray(groupData, 'moviesType', 'Flashed_Gratings', @strcmp);
% %     fgGids = [groupData(flashGratingGroupInds).Gid];
% 
%     for si = 1:length(s)
%         nm = s{si};
%         [Gid, cellId] = getGidCellIdFromVarname(nm);
%         
%         v  = eval(nm);
%         if isstruct(v) && ~isempty(v.OSP)
%             [OSP, oris, sps, phs] = elements(v.OSP);
%             if equal(phs, phsSet_rad)
%                 phs = round(rad2deg(phsSet_rad));
%                 v.OSP{4} = phs;
%                 disp(['Fixed ' nm]);
%             end
%         end
%         eval([nm ' = v;']);
% 
%     end
% end

% change variable names from osp-->OSP, and make into struct.
% load indivCells_movie_fg

% s = who([varnamePrefix '_*']);
% progressBar('init-', length(s), 50);
% for si = 1:length(s)
%     progressBar(si);
%     nm = s{si};
%     celldata = eval(nm);
%     if isstruct(celldata)
%         if ~isstruct(celldata.PSTH)
%             if length(celldata.PSTH) == 5
%                 [PSTH_bins, PSTH_vals, frameLength_ms, bckgRate, timeWindow_ms] = elements(celldata.PSTH);
%                 celldata.PSTH = struct('bins', PSTH_bins(:), 'vals', PSTH_vals(:), 'frameLength_ms', frameLength_ms, 'bckgRate', bckgRate, 'timeWindow_ms', timeWindow_ms);
%             elseif length(celldata.PSTH) == 4
%                 [PSTH_bins, PSTH_vals, frameLength_ms, bckgRate] = elements(celldata.PSTH);
%                 celldata.PSTH = struct('bins', PSTH_bins(:), 'vals', PSTH_vals(:), 'frameLength_ms', frameLength_ms, 'bckgRate', bckgRate);
%             elseif length(celldata.PSTH) == 3
%                 [PSTH_bins, PSTH_vals, bckgRate] = elements(celldata.PSTH);
% %                 frameLength_ms = getFrameLength('Gid', celldata.Gid);
%                 celldata.PSTH = struct('bins', PSTH_bins(:), 'vals', PSTH_vals(:), 'bckgRate', bckgRate);
%             elseif (~isempty(celldata.PSTH))
%                 keyboard;
%             end
%         end
%         
%         if isfield(celldata, 'osp');
%             celldata = renameStructField(celldata, 'osp', 'OSP');
%             if ~isempty(celldata.OSP) && ~isstruct(celldata.OSP)
%                 [R, oris, sps, phs] = elements(celldata.OSP);
%                 celldata.OSP = struct('R', R, 'ori', oris, 'sp', sps, 'ph', phs);
%             end
%         end
%         
%         if ~isstruct(celldata.STAs)
%             STAC = celldata.STAs;
%             [STA, timeWindows] = elements(STAC);
%             celldata.STAs = struct('STA', STA, 'timeWindow_ms', timeWindows);
%         end        
%         eval([nm ' = celldata;']);
%     end
% 
% end
% save indivCells_movie_fg2 celldata_*
    
% 1. Fix OSP stats - add in spatialfrequencySelective,
% spatialfrequencySelectivePval fields
% if true
%     s = who('-regexp', [varnamePrefix '*']);
%     F = {'spatialfrequencySelective', 'spatialfrequencySelectivePval'};
%     
%     for si = 1:length(s)
%         nm = s{si};
%         [Gid, cellId] = getGidCellIdFromVarname(nm);        
%         v = eval(nm);
%         if isstruct(v) && ~isempty(v.OSP)
%             ospstats = v.OSP.stats;
%             if ~isfield(ospstats, F{1})
%                 ospstats.(F{1}) = NaN;
%             end
%             if ~isfield(ospstats, F{2})
%                 ospstats.(F{2}) = NaN;
%             end
%             ospstats = orderfields(ospstats, names);
%             v.OSP.stats = ospstats;            
%         end
%         eval([nm ' = v;']);
% 
%     end
% end


% 1. Fix OSP stats - make all true/false ones double class
% spatialfrequencySelectivePval fields
%{
if false
    s = who('-regexp', [varnamePrefix '*']);
    F = {'orientationSelective', 'spatialfrequencySelective', 'orientationReproducible', 'spatFreqReproducible', 'responseSize'};
    G = {'responseSizeFrac'};
    for si = 1:length(s)
        nm = s{si};
        [Gid, cellId] = getGidCellIdFromVarname(nm);        
        v = eval(nm);
        nl = false;
        if isstruct(v) && ~isempty(v.OSP)
            ospstats = v.OSP.stats;
            for fi = 1:length(F)         
                if islogical(ospstats.(F{1}))
                    ospstats.(F{1}) = double(ospstats.(F{1}));
                    fprintf('converted to double ... ');
                    nl = true;
                end
            end
            if ~isfield(ospstats, G{1})
                ospstats.(G{1}) = NaN;
                fprintf(' ... added field');
                nl = true;
            end
            if nl
                fprintf('\n');
            end                
            v.OSP.stats = ospstats;            
            eval([nm ' = v;']);
        end
        

    end
end

%}