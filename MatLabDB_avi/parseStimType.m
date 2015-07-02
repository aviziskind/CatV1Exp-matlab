function  [dims, repType, frmLength_ms] = parseStimType(stimType) 
    A = sscanf(stimType, '%dx%dx%d(%dx%d)')';
    dims = A(1:3);
    repType = A(4:5);
    [tmp, frm_str] = strtok(stimType, '_');
    frmLength_ms = str2double(frm_str(2:end));    
end