function s = experimentDBPath
    if ispc
        s = 'D:\ExperimentDB\';
%         s = 'C:\ExperimentDB\';
    else
        if onNYUserver
            s = '~/';
        else
            s = '/media/avi/Storage/ExperimentDB/';
        end
    end

end