% function testGetBestWindowFromPSTH

% before running this mfile:
%     load movieCells_PSTHs;
%     load cellsGroups_movie_fg;
    longFrameInds  = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x > 40);
    shortFrameInds  = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x < 25);
    allInds        = 1:length(movieGroups_fg);
    indsToDo = shortFrameInds; 
%   indsToDo = longFrameInds; 
    
    groupData = movieGroups_fg(indsToDo);    
    
	suppressionCells = [4466,1; 4476,2; 4482,3; 4482,6; 4494,1; 4506,1; 4520,3; 4522,1; 4712,2; 4712,3; 4796,1; 4798,1; 5084,2; 5112,2];
    noResponseCells =  [4476,3; 4482,4; 4502,1; 4520,1; 4706,2; 4718,4; 4732,1; 4744,3; 4890,1];
    borderlineCells =  [4462,1; 4462,3; 4466,3; 4470,5; 4488,1; 4506,0; 4506,3; 4508,1; 4508,2; 4522,2; 4712,4; 4718,1; 4726,1; 4732,2; 4732,3; 4882,2; 4890,0; 5112,3;];
    specialCases =     [4503,3; 4520,0; 4546,2; 4796,4; 4494,2; 4506,2];

    redoAll = true;
    pauseAfterEachPage = true;
    
    nCells = sum( cellfun(@length, {groupData.cellIds} ) );
%     categories = {'Classic', 'Suppressed', 'noResponse'};
    if redoAll || ~exist('PSTHs', 'var');
        clear('PSTHs')
        PSTHs = repmat(struct('Gid', [], 'cellId', [], 'bins', [], 'vals', [], 'frameLength_ms', [], 'bckgRate', [], 'timeWindow', [], 'msg', [], 'category', []), nCells, 1);
        ind = 1;

        % Stage 1: Gather data
        for Gid_i = 1:length(groupData)
            
            Gid = groupData(Gid_i).Gid;
%             groupDataNm = ['groupData_' num2str(Gid)];
%             if ~exist(groupDataNm, 'var')
%                 grpData = movieGroups_fg(findInStructArray(movieGroups_fg, 'Gid', Gid)); %#ok<NASGU>
%                 eval([grpDataNm ' = groupData;']);
%             else
%                 groupData = eval(groupDataNm);
%             end

            cellIds = groupData(Gid_i).cellIds;
            for cell_i = 1:length(cellIds)

                if (Gid == 5276) && (cellIds(cell_i) == 5)
                    3;
                end
                
                varname = [getName('celldata', Gid, cellIds(cell_i)) '_PSTH'];
                if ~exist(varname, 'var') || ~isstruct(eval(varname))                        
                    break;
                end
                psthData = eval(varname);
                [PSTH_bins, PSTH_vals, frameLength_ms, bckgRate] = elements(psthData);
    %             fprintf(['(' num2str(sub_fig) ') ']);
                [timeWindow, windowProfile, msg] = getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals );
                category = 'z';
                if any (  (Gid == suppressionCells(:,1)) & ( cellIds(cell_i) == suppressionCells(:,2)) )
                    category = 'Sup';
                end
                if any (  (Gid == noResponseCells(:,1)) & ( cellIds(cell_i) == noResponseCells(:,2)) )
                    category = 'NR';
                end
                if any (  (Gid == borderlineCells(:,1)) & ( cellIds(cell_i) == borderlineCells(:,2)) )
                    category = 'bor';
                end
                if any (  (Gid == specialCases(:,1)) & ( cellIds(cell_i) == specialCases(:,2)) )
                    category = 'spec';
                end

                PSTHs(ind) = struct('Gid', Gid, 'cellId', cellIds(cell_i), 'bins', PSTH_bins, 'vals', PSTH_vals, ...
                    'frameLength_ms', frameLength_ms, 'bckgRate', bckgRate, 'timeWindow', timeWindow, ...
                    'msg', msg, 'category', category);
                ind = ind + 1;            
            end


        end
        
        PSTHs(ind:end)=[];
        PSTHs = sortStructArray(PSTHs, 'msg', 'descend');
    %     PSTHs = sortStructArray(PSTHs, 'msg');

    end % of redoing all window calculations
    
    nCells = ind-1;
    disp(['Total of ' num2str(nCells) ' cells']);
    % Stage 2: Plot data
    
    fig_id = 1;
    mxPages = 10;
    gridSubPlot(3,5, [fig_id, mxPages])
    N = 105;
    curCat = '';
    numNs = ceil(length(PSTHs)/N);
    whichN = 5;
    
    for PSTH_i = (whichN-1)*N+1 : min(length(PSTHs), (whichN)*N);

    %     for PSTH_i = 1:length(PSTHs);
        [Gid, cellId, PSTH_bins, PSTH_vals, frameLength_ms, bckgRate, timeWindow, msg, category] = ...
            elements( PSTHs(PSTH_i) );

        lastOne = gridSubPlot;
        cla;

        plotThisPSTH(PSTH_bins, PSTH_vals, [], [], timeWindow);
        getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals, 'samePlot');
%         title([ '(' num2str(PSTH_i) ') ' num2str(Gid) ', ' num2str(cellId)  '  [' num2str(sub_fig) ']' iff(isempty(msg), '', ['*:' msg ])]);
        title([ '(' num2str(PSTH_i) ') ' num2str(Gid) ', ' num2str(cellId)  , ' [' num2str(msg) '] ' iff(isempty(category), '', [ ':' category ])]);

        if lastOne
            break;
        end
%         if any(PSTH_i == [27, 29, 51, 62, 67])
%           [timeWindow, windowProfile, msg] = getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals, bckgRate, PSTH_i);
%           delay_ms = timeWindow(1);
%           windowSize_ms = diff(timeWindow);
%         else        
%             getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals, bckgRate, 'flag');
%         end
        
%         if (PSTH_i == 7)
%             getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals, bckgRate, 'flag', 'flag2');
%         end
        
%         switch category
%             case {'', 'bor', 'spec'}
%                 drawVerticalLine([delay_ms], 'Color', 'c');
%                 drawVerticalLine([delay_ms + windowSize_ms], 'Color', 'r');                
%                 
%         end        
        
    end                  
                  
                  
                  
% end



%             if length(psthData.PSTH) == 3
%                 [PSTH_bins, PSTH_vals, bckgRate] = elements(psthData.PSTH);
%                 frameLength_ms = getFrameLength('Gid', Gid);
%                 psthData.PSTH = {PSTH_bins, PSTH_vals, frameLength_ms, bckgRate};
%                 eval([varname ' = psthData;']);
%             elseif length(psthData.PSTH) == 4
%                 [PSTH_bins, PSTH_vals, frameLength_ms, bckgRate] = elements(psthData.PSTH);
%             end
