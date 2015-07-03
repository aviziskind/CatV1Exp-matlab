function hnd = dbOpenExpDb
%     return;
    if ~isOldXPS
        error('!');
    end
    dbAccount  = 'SYSDBA';
    dbPassword = 'masterkey';
    hnd = edbOpen(dbAccount, dbPassword);
    
end