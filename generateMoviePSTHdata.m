% generateMoviePSTHdata

clear all
load indivCells_movie_fg
s = who;
for i = 1:length(s)
    celldata = eval(s{i});
    if isstruct(celldata)
        eval([s{i} '_PSTH = celldata.PSTH;']);
    end
    clear(s{i});
end
save movieCells_PSTHs celldata*


