% viewCalculationsForAllCells

longFrameGids = [ 4462, 4463, 4466, 4467, 4470, 4471, 4476, 4477, 4482, 4483, 4488, 4489, 4494, 4495, ...
                  4506, 4507, 4508, 4509, 4520, 4521, 4522, 4523, 4546, 4547, 4548, 4549, 4706, 4707, ...
                  4712, 4713, 4718, 4719, 4726, 4727, 4732, 4733, 4744, 4745, 4796, 4797, 4798, 4799, ...
                  4882, 4883, 4890, 4891, 5084, 5085, 5112, 5113 ];

    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    typesToView = [2];
    
    for ti = typesToView
%         if ~exist([stimTypes{ti} 'Groups'], 'var')
            S = load(['cellsGroups_' stimTypes{ti} '_DB']);
            groupData = S.([stimTypes{ti} 'Groups']);
            eval([stimTypes{ti} 'Groups = groupData;']);
%         end

        allGidsToDo = [groupData.Gid]';
%         allGidsToDo = [movieGroups.Gid]';
%         allGidsToDo = longFrameGids;
        ind =  find( allGidsToDo == 1481, 1, 'first');
        for Gid_i = ind:length(allGidsToDo);
            Gid = allGidsToDo(Gid_i);
            thisGroupData = groupData(findInStructArray(groupData, 'Gid', Gid) );
            cellIds = groupData.cellIds;
            if ( all(thisGroupData.presOK) )
            
                for cellId = cellIds;
                    varname = getName('celldata', Gid, cellId);
                    oktxt = iff(all(thisGroupData.presOK), '-', '!');
                    fprintf(['[' outOf(Gid_i, length(allGidsToDo)) '] [' oktxt ']']);

                    cellData = [];
                    if exist(varname, 'var')
                        cellData = eval(varname);
                    end

                    if isstruct(cellData)
                        plotCellData(cellData, 200);
    %                     disp(groupData.movieFiles);
                        input('');
                    end

                end  % of loop for cell      
            end
        end % of loop for each group
        
    end % of loop for each stimulus type
                       

% end

% possible noise ones with STA: 1481;0(?) 1481;1 (!!) 1881;0(?) 1881;1(?)  1883;0(?) 1883;1 (60-90)
% 2030;0/1(?)  2030;2(?)   2078;1 (!), 2111;0, 2112;0, 2124;0, 4189;1

% good psth: 916;2,  