function tmp_copyAllMatInSubdir
curdir = cd;
allDirs = dir('*');
allDirs = allDirs([allDirs.isdir]);
allDirs = allDirs( arrayfun(@(s) isempty( strfind(s.name, '.')), allDirs));

for i = 1:length(allDirs)
    cd([curdir '\' allDirs(i)]);
    copyfile('*.mat', '..');
end