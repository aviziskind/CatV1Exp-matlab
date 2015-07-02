
% function assignTunednessRanksManually

%     cellsToRedo = [2272, 5;  2841,1;  4486, 4;   4356, 3;    4730, 0;   2903, 0];
% cellsToRedo = [1875, 1; 1929, 0; 2066, 0; 1857, 1; 2214, 0; 2216, 0;  2288, 7; 2272,9; 2336,2;  2346,3; 2358,3; 2546, 1; 2544,3; 4878, 0; 5086, 2; 5114, 1; 5282,6; 5184,4];
% cellsToRedo = [...
%  Group 1494,  cell 2
% ];
    global showWorking

    gratingType = curGratingType('');  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
    ospDatafile = [CatV1Path gratingType 'GratingOSPs.mat'];    
    
    requireEnter = false;

    sortBy = 'rep_ori_sp_avPhase_pval';
    
    S = load(ospDatafile);
    allOSPs = S.allOSPs;    

    % restrict to nOri > 2 
    allOSPs = allOSPs(findInStructArray(allOSPs, 'ori', [], @(x) length(x) > 2));
    
    % sort by stats (response reproducible)
    allstats = [allOSPs.stats];
    rep_p = -log10([allstats.(sortBy)]);
    idx = ord(rep_p, 'descend');    
    allOSPs = allOSPs(idx);
    
    % gather list of group/cellids
    allGidCellIds = [[allOSPs(:).GroupId];  [allOSPs(:).cellId]]';
    nCells = length(allOSPs);

    % Retrieve existing rank info
    rankFilename = ['cellRanks_' gratingType '.mat']; 
    if exist(rankFilename, 'file')
        load(rankFilename);
    else
        ranks = repmat(struct('Gid', [], 'cellId', [], 'rank', []), nCells, 1);
    end
    
    % find out which ones have been done.
    GidsDone = [ranks(:).Gid]; if iscell(GidsDone), GidsDone = [GidsDone{:}]; end
    cellIdsDone = [ranks(:).cellId]; if iscell(cellIdsDone), cellIdsDone = [cellIdsDone{:}]; end
    GidCellIdsDone = unique([GidsDone(:), cellIdsDone(:)], 'rows');
    if exist('cellsToRedo', 'var') && ~isempty(cellsToRedo)
        [GidCellIdsDone, i] = setdiff(GidCellIdsDone, cellsToRedo, 'rows'); 
    end
            
    GidCellIdsToDo = allGidCellIds;
    assert(size(allGidCellIds,1) == size(unique(allGidCellIds, 'rows'), 1));
    if ~isempty(GidsDone)
        [tmp, idx] = setdiff(GidCellIdsToDo, GidCellIdsDone, 'rows');
        GidCellIdsToDo = GidCellIdsToDo(sort(idx),:);         
    end
    nDone = length(GidsDone);
    nToDo = size(GidCellIdsToDo, 1);

    
    
    
%     if nDone < 5000
%         progressBar('init=', nToDo);
        for i = 1:nToDo
%             progressBar(i);
            [Gid, cellId] = elements(GidCellIdsToDo(i,:));
            ind = findRows( [Gid cellId], allGidCellIds);        
            OSP = allOSPs(ind);
            [Gid2, cellId2] = deal(OSP.GroupId, OSP.cellId);
            assert(all([Gid, cellId] == [Gid2 cellId2]));
            [nOri, nSpf, nPh] = size(OSP.R);
            R = OSP.R;
            oris = OSP.ori - min(OSP.ori);            
            if strcmp(gratingType, 'flashed')
                figure(12); clf;
                [ori_pref] = getPreferredOriFromOSP(sum(squeeze(R),2), oris, 180, 1);                    
                ori_wgt_vector = cos( 2*( deg2rad( oris - ori_pref) ) );
                if nSpf == 1                    
                    subplot(1,2,1);
                    R_ori = sum(squeeze(R), 2);
                    R_js = R_ori .* ori_wgt_vector;
                    hist(R_js); drawVerticalLine(0, 'linestyle', '-', 'color', 'g', 'linewidth', 2);
                    m = mean(R_js); s = stderr(R_js);
                    drawVerticalLine(m, 'color', 'r');
                    drawVerticalLine([m-s, m+s], 'linestyle', ':', 'color', 'r');
%                     [ori_sel, ori_sel_pval] = ttest(R_js_m(:), 0, alpha, 'right'); % perform t-test to determine whether set of Rijs has mean signficantly > 0, with significance level alpha.
                    subplot(1,2,2);
                end                                
                imageOSP(OSP, 'mean:ph');
            elseif strcmp(gratingType, 'drifting')
%                 imageOSP(OSP, 'mean:ph');
                figure(12); clf;
                if any([nOri nSpf] == 1)
                    imageOSP(OSP); 
                else
                    subplot(1,2,1);
                    imageOSP(OSP, 'mean:ph', 'OSP', 'nolabels');  title('av over phases');
                    subplot(1,2,2);                    
                    imageOSP(OSP, 'pref:ori', 'SP', 'nolabels');  title('at pref. dir');
                    
                    if nOri <= 6
                        figure(13); clf;
                        imageOSP(OSP, 'subplots:horizontal', 'SPO');
                    end
%                     imageOSP(OSP, 'pref:sp', 'OP'); title('at pref. sp');                     
                end
                
                showWorking = 'rep_ori_sp_ph';
%                 nm = getName('celldata', Gid, cellId);
%                 S = getfield(eval(nm), 'OSP'); %#ok<GFLD>
%                 calcStatsFromOSPfull( S.R, S.R_full, S.R_full_factor, Gid, []);
                3;
            end
%             figure(12); title(sprintf('Group %d,  cell %d', Gid, cellId));
                                
            
            colorbar;
%             figure(13); clf;
%             images(

            fprintf('[%dx%dx%d][ %s ] Rank for Group %d,  cell %d  [0/1/2/3]: ', nOri, nSpf, nPh, outOf(nDone+i, nDone+nToDo), Gid, cellId);
            
            if requireEnter
                s = input('', 's');
                r = str2num(s);
            else
                k = -1;
                while k ~= 1
                    k = waitforbuttonpress;
                    s = get(gcf, 'CurrentCharacter');
                    r = str2num(s);     %#ok<ST2NM>
                    fprintf([s '\n']);
                end
            end
            
            if isempty(s) % skip
                continue;
            elseif ~isempty(r)  % user entered a number
                if any(r == [0 1 2 3])
                    idx = find(Gid == [ranks.Gid] & cellId == [ranks.cellId]);
                    if ~isempty(idx)     % check if cell is already ranked -- if so, replace rank
                        ranks(idx).rank = r;
                    else  % otherwise, add to end of array
                        idx2 = length(ranks)+1;
%                         assert(length(ranks) < idx2 || isempty(ranks(idx2).Gid));
                        ranks(idx2) = struct('Gid', Gid, 'cellId', cellId, 'rank', r);
                    end
                else
                    fprintf ('[input not accepted]\n');
                end
                
            elseif strcmpi(s, 'x');
                save(rankFilename, 'ranks');
                return;
            end        
            save(rankFilename, 'ranks');
        end
        
%         progressBar('done');
%     end    
    fprintf('Completed rank assignments for %s gratings...\n', gratingType);
    
    %sort the array, then resave
    [tmp, newOrder] = sortrows([[ranks.Gid]', [ranks.cellId]']);
    ranks = ranks(newOrder); %#ok<NASGU>
    save(rankFilename, 'ranks');
    
    % learn the parameter space:
    
%     function e = errorOnOSP(sOSP, params)
%         [sig_x, sig_y, ]
%         gaussian2D
%         
%         
%     end


% end

    
% 542 2 