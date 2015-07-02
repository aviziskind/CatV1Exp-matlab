function tableName = getDatabaseTableForDid(Did, stimType)

    if nargin < 2        
        stimType = getStimulusTypeForDid(Did);
    end
    tableName = ['TBL_' upper(stimType) '_PRES'];
    
end
