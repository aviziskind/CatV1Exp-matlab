function fsize = getRawDataFileSize(dfile)

    persistent dataFilesInfo
    redo = 0;
    
    dfsizesFileName = [CatV1Path 'MatLabDB_avi' filesep 'datafileSizes.mat'];
    
    if isempty(dataFilesInfo)        
        if exist(dfsizesFileName, 'file') && ~redo  % load if exists
            S = load(dfsizesFileName);
            dataFilesInfo = S.dataFilesInfo;
        
        else  % recalculate (make sure raw data hard drive is connected)        
            dataFilesInfo = calculateDataFilesInfo;            
            save(dfsizesFileName, 'dataFilesInfo');
        end
    end
            
    
    if isnumeric(dfile) 
        idx = find(dataFilesInfo.Dids == dfile, 1);
        if isempty(idx)
            error('No Did = %d found', dfile);
        end
        fsize = dataFilesInfo.dfsizes(idx);        
        
    elseif ischar(dfile)
        idx = find(strcmp(dfile, dataFilesInfo.dfnames),1);
        if isempty(idx)
            error('No datafile %s found', dfile);
        end
        fsize = dataFilesInfo.dfsizes(idx);                
        
    end    
    
end
                    

function dfInfo = calculateDataFilesInfo
    hnd = dbOpenExpDb;
    [tbl_Dids, tbl_dfnames] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'TXT_DATAFILE_NAME'}, 'TBL_DATA_FILES', [], 'DATAFILE_ID');
    nDataFiles = length(tbl_Dids);
    dfsizes = zeros(nDataFiles,1);
    progressBar('init-', nDataFiles);
    for di = 1:nDataFiles
        progressBar;
        fullpath = fullDatafilePath(tbl_dfnames{di});
        if exist(fullpath, 'file')
            dfsizes(di) = filesize(fullpath);
        else
            dfsizes(di) = nan;
        end
    end

    dfInfo = struct('Dids', tbl_Dids, 'dfnames', {tbl_dfnames}, 'dfsizes', dfsizes);
end

function s = fullDatafilePath(dfname)
    external_hd = 'G:\RawData\';
    catName = strtok(dfname, '_');
    s = [external_hd catName '\' dfname];
end