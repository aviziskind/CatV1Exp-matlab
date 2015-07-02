% function generateAllMseqStims
    
    [tapRegisterCol, nMseqBitsCol] = getFieldsFromDatabaseTable(dbOpenExpDb, {'LNG_TAP_REGISTER', 'LNG_MSEQ_BITS'}, 'TBL_MSEQ_PRES');
    
    [tapRegisterCol, nMseqBitsCol] = elements(unique([tapRegisterCol, nMseqBitsCol], 'rows'));

    for i = 1:length(tapRegisterCol)
        tapRegister = tapRegisterCol(i);
        nMseqBits = nMseqBitsCol(i);                
        generateMsequenceStim(tapRegister, nMseqBits);
    end
            
        
