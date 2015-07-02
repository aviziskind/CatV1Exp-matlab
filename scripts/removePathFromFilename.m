function fn = removePathFromFilename(filename)
    [~, fname, fext] = fileparts(filename);
    fn = [fname fext];

end