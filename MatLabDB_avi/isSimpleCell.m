function [tf, F1oDC] = isSimpleCell(varargin)

    dbug = 0+false;
    
    [F1, DC] = getF1oDC(varargin{:});
        
    tf = F1 > DC;
    F1oDC = F1/DC;

    if dbug
        s = iff(tf, 'simple', 'complex');
        disp(['F1 = ' num2str(F1) ]);
        disp(['DC = ' num2str(DC) ]);
        disp(['F1/DC = ' num2str(F1/DC) ' (' s ' cell)'  ]);
    end;

end

