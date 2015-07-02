function exploreTcPowerSpectra


    [gratingType, gratingType_s] = curGratingType;      
    ospDatafile = [CatV1Path gratingType_s 'GratingCells_DB.mat'];    
    
    S = load(ospDatafile);
    allCells = S.allCells;
    
    nCells = length(allCells);

    nPh = arrayfun(@(s) length(s.ph), allCells);
    [uPh, phList] = uniqueList(nPh);
    pSpc = cell(1, length(uPh));
    
    skipDC = 0;
    figure(1); clf; hold on;
    curXmax = 0;
    maxHarmonic = 10;
    cols = ['br'];
    spc = .1;
    for ph_i = 1:length(uPh)
        pSpc{ph_i} = zeros(length(phList{ph_i}),  1);
        
        for ci = 1:length(phList{ph_i})
            C = allCells(phList{ph_i}(ci));
            ori_sp = C.oriSp_maxR_av;
            tc = C.R(ori_sp(1), ori_sp(2), :); tc = tc(:);        
        
            [f, p] = powerSpectrum(tc);
            pSpc{ph_i}(ci,1:length(p)) = p;
        end        
        
        
        PS = pSpc{ph_i};
%         n = size(PS,2);        
        
        if skipDC
            x_start = 1;
        else           
            x_start = 0;
        end
        if size(PS,2) > maxHarmonic
            PS = PS(:, 1:maxHarmonic+1);
        end        
        L = size(PS,2);  
        PS = PS(:, x_start+1:end);
        M = mean(PS,1); S = stderr(PS,1); 
        errorbar([x_start:L-1] - spc/2 + spc*(ph_i-1), M/M(1), S/M(1), [cols(ph_i) '.-'] );        
        if L-1 > curXmax;
            xlim([x_start-spc*2, L-1+spc*2]);
        end
        curXmax = L-1;
        ylims = ylim;
        ylim([0 ylims(2)]);        
        
    end
%     title(sprintf('# phases = %d', uPh(ph_i)));
    xlabel('harmonic'); ylabel('relative power');
    legend(legendarray('nPh = ', uPh));
    title([iff(skipDC, 'Without', 'With') ' DC']);


end