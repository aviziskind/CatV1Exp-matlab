function generateMsequenceStim(tapRegister, nMseqBits)
                    
    nBits = (2^nMseqBits - 1)*nMseqBits;

    msequenceGenerator('init', tapRegister, nMseqBits);    
    sequence = msequenceGenerator(nBits,1);

    filename = getName('mseqStimFile', tapRegister, nMseqBits);
    fileId = fopen(filename, 'w');
    fwrite(fileId, sequence, 'uint8');
    fclose(fileId);
    
end

% function generateMsequenceMovie(frameDims, nFrames, tapRegister, nMseqBits)
%                     
%     nCycleBits = (2^nMseqBits - 1)*nMseqBits;
%     nTotalPixels = prod(frameDims)*nFrames;
%     numSeqReps = ceil(nTotalPixels/nCycleBits);
% 
%     msequenceGenerator('init', tapRegister, nMseqBits);    
%     sequence = msequenceGenerator(nCycleBits,1);
% 
%     filename = getName('mseqStimFile', tapRegister, nMseqBits);
%     fileId = fopen(filename, 'w');
%     
%     progressBar('init-', nFrames);
%     for ri = 1:numSeqReps
%         progressBar(ri);        
%         fwrite(fileId, sequence, 'uint8');
%     end    
% 
%     fclose(fileId);
%     
% end