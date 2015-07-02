function mergeNextGroup(i_in)

    persistent Gids SD dfDates uLocIds allLocIds locMaxNSpks_ordered nspks locIdx locId1 locIdList i locIdCount

    if isempty(Gids)

        Gids = sort(getAllGids );
        SD = siteDataFor('Gid', Gids, 1);
        allLocData = [SD.locationData];
        allLocIds = [allLocData.LocId];
        
        dfData = [SD.dataFileInfo];        
        dfDates = {dfData.dateCreated};
                
        nspks = arrayfun(@(s) sum(s.nSpikes(2:end)), SD);        

%         idx = ord(nspks, 'descend');
        [uLocIds, locIdList] = uniqueList(allLocIds);
        
%         GidsAtLocs  = cellfun(@(idx) Gids(idx), locIdList, 'un', 0);        
        dfDatesAtLocs = cellfun(@(idx) dfDates(idx), locIdList, 'un', 0);        
        newOrds = cellfun(@(G) ord(G), dfDatesAtLocs, 'un', 0);
%         dateStrFormat = 'yyyy/mm/dd HH:MM:SS';
%         datenums = cellfun(@(s) datenum(s, dateStrFormat), dfDates);
%         mds = cellfun(@(ids) max(pdist(datenums(ids)')), locIdList, 'un', 0);
        
        3;
        locIdList = cellfun(@(ids, idxs) ids(idxs), locIdList, newOrds, 'un', 0); 
               
        locIdCount = cellfun(@length, locIdList);
        locMaxNSpks = cellfun(@(idxs) max(nspks(idxs)), locIdList);

        [locMaxNSpks_ordered, loc_newOrder] = sort(locMaxNSpks, 'descend');        
        
        uLocIds = uLocIds(loc_newOrder);
        locIdCount = locIdCount(loc_newOrder);
        locIdList = locIdList(loc_newOrder);        
        
        locId1 = cumsum([1 locIdCount(1:end-1)]);
        
        
    end
    
    if isempty(i)
        i = 1;
    end    
    
    if (nargin == 1) && ~isempty(i_in)        
        i = i_in;             
    end
    
    firstI = true; fn = '';
    while firstI || exist(fn, 'file')            
        locIdx = find(i >= locId1, 1, 'last');
        locSubIdx = i - locId1(locIdx)+1;
        curGid = Gids(locIdList{locIdx}(locSubIdx));
        curDate = dfDates{locIdList{locIdx}(locSubIdx)};
        curLocMaxNSpk = locMaxNSpks_ordered(locIdx);
        fn = getFileName('clusterMerges', curGid);
        firstI = false;
        if exist(fn, 'file')
            i = i + 1;
        end
    end            
        
    curLocId = uLocIds(locIdx);    
    
    if locSubIdx > 1
        prevDate = dfDates{locIdList{locIdx}(locSubIdx-1)};
        dateStrFormat = 'yyyy/mm/dd HH:MM:SS';
        datenum_cur = datenum(curDate, dateStrFormat);
        datenum_prev = datenum(prevDate, dateStrFormat);
        diff_str = ['(+' sec2hms( (datenum_cur-datenum_prev)*(24*60*60) ) ') ']; %datestr(datenum_cur-datenum_prev, 'HH:MM:SS')
    else
        diff_str = '';
    end
    
    fprintf(' (%d/%d) LocId %d [max:%d] (%d/%d) Group %d (%d / %d at this location). [%s%s]\n', ...
        i, length(Gids), curLocId, curLocMaxNSpk,  locIdx, length(uLocIds), curGid, locSubIdx, locIdCount(locIdx), diff_str, curDate );

    findBestClusterMerges(curGid);


end