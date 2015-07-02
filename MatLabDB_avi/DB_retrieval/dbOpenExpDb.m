function hnd = dbOpenExpDb
%     return;
    if isXPS15
        error('!');
    end
    dbAccount  = 'SYSDBA';
    dbPassword = 'masterkey';
    hnd = edbOpen(dbAccount, dbPassword);
    
end