function s = experimentDBPath
    if ispc
        if isXPS15
            s = 'D:\ExperimentDB\';
        elseif isOldXPS
            s = 'C:\ExperimentDB\';
        end
    else
        if onNYUserver
            s = '~/';
        else
            s = '/media/avi/Storage/ExperimentDB/';
        end
    end

end