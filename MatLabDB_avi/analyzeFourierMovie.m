function analyzeFourierMovie(filename)
    sp_pix = [128   63.0000   31.5000   21.0000   15.7500   12.6000   10.5000    9.000    7.8750    7.0000    6.3000 ...
            5.7273    5.2500    4.8462    4.5000    4.2000    3.9375    3.7059    3.5000    3.3158    3.1500    3.0000 ...
            2.8636    2.7391    2.6250    2.5200    2.4231    2.3333    2.2500    2.1724    2.1000    2.0323    2.0000 ...
            2.0323    2.1000    2.1724    2.2500    2.3333    2.4231    2.5200    2.6250    2.7391    2.8636    3.0000 ...
            3.1500    3.3158    3.5000    3.7059    3.9375    4.2000    4.5000    4.8462    5.2500    5.7273    6.3000 ...
            7.0000    7.8750    9.0000   10.5000   12.6000   15.7500   21.0000   31.5000   63.0000];
    
    L = 64;
    Tr = dftmtx(L)./sqrt(L);  % Fourier basis

    function frm = generateFrameOrig(arg1, arg2, arg3)

        frm = zeros(L,L);
        frm(arg2,arg1) = arg3;
        frm = real(Tr'*frm*Tr)+imag(Tr'*frm*Tr);
%         frm = real(Tr'*frm*Tr) ;

    end
    
    function frm = generateFrame(x, y, ph)

        % get sp:
        sp_freq = norm([x-1, y-1]);
        sp_period = L/sp_freq;        
        
        % get ori:
        ori_deg = rad2deg (atan2(x-1,y-1));
        
        ph_deg = 45;%ph*90 + 47;  %[46 -> 51]; [43->48] --> [46-48]        
        
        frm = generateGratingFrame([L L], ori_deg, sp_period, ph_deg);       
        
%         frm = fliplr(flipud(frm));
    end
    
    if nargin < 1
        filename = upper('fgcomp_2x64x64x8192_1.raw');
    end

    hnd = dbOpenExpDb;
    [movieId, nCols, nRows, nFramesTotal] = getFieldsFromDatabaseTable(hnd, {'MOVIE_ID', 'LNG_N_COLUMNS', 'LNG_N_ROWS', 'LNG_N_FRAMES'}, 'TBL_MOVIES', ...
        {'TXT_MOVIE_FILE_NAME', upper(filename)});
    
    [spPhase, spFreq, ori] = getFlashGratingOSPs_fgmovie(hnd, movieId);
    Gid = 3091;
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid); 
    uPh = unique(spPhase);
    uSp = unique(spFreq);
    uOri = unique(ori);
    
    frms = find( (spPhase == uPh(1)) & (spFreq == uSp(64)) & (ori == uOri(64)) );
%     spFreqvalInds = ord ( spFreq(frms) );
%     frms = frms(spFreqvalInds);

    getSpatialFreqs = false;
    getOrientations = false;
    doInteractive = true;
% get list of spatial frequencies 
    if getSpatialFreqs
        sp_pix = zeros(1,length(uSp));
        for sp_i = 1:length(uSp)
            frmId = find( (spPhase == uPh(1)) & (spFreq == uSp(sp_i)) & (ori == uOri(1)) );
            figure(1); clf;
            thisFrame = getFrame(frmId);
            imagesc(thisFrame); colormap('gray');
            title(sprintf('(%d/%d)Frame #%d [%d, %d, %d]', sp_i, length(uSp), frmId, spPhase(frmId), spFreq(frmId), ori(frmId) ));
            set(gca, 'xtick', [], 'ytick', []);
            axis square;

            t = 1:64; f = thisFrame(:,1); f = f-mean(f);
            [strongestFrequencies] = findStrongestFrequencies(t,f, 1);
    %         strongestFrequencies = strongestFrequencies
            sp_pix(sp_i) = 1/strongestFrequencies;
        end
        getFrame('close'); 
        figure(2);
        plot(t, log(sp_pix), '.');
    end
    
    function frm = rescaledFrame(frm)
        mx = max(frm(:));
        mn = min(frm(:));
        d = mx-mn;
        if mx > mn
            frm = (frm - mn)/d;
        else
            frm = .5 * ones(size(frm));
        end        
    end
    
    if getOrientations
        sp_ind = 10;
        ph_ind = 1;
        ori_deg = zeros(1,length(uOri));
        for ori_i = 4:length(uOri)
            frmId = find( (spPhase == uPh(ph_ind)) & (spFreq == uSp(sp_ind)) & (ori == uOri(ori_i)) );
            figure(1); clf;
            thisFrame = getFrame(frmId);
            thisFrame = rescaledFrame(thisFrame);
            imagesc(thisFrame); colormap('gray'); 
            title(sprintf('(%d/%d)Frame #%d [%d, %d, %d]', ori_i, length(uOri), frmId, spPhase(frmId), spFreq(frmId), ori(frmId) ));
            set(gca, 'xtick', [], 'ytick', []);
            axis square xy;
%             c1 = thisFrame(1:32,1);
%             i1 = indmax(c1);
%             drawHorizontalLine(i1);
                        
            frm = generateFrameOrig(ori_i, sp_ind, uPh(ph_ind));
            frm = rescaledFrame(frm);
%             frm = frm / max(frm(:));
            
            figure(2);
            imagesc(frm); axis square xy; colormap('gray');
%             c2 = frm(1:32,1);
%             i2 = indmax(c2);
%             drawHorizontalLine(i2);
%             is = [i1, i2]
            
            figure(3); clf;
            imagesc(frm - thisFrame); axis square xy;  colormap('gray'); colorbar;
            
            return;
            3;            
        end
        getFrame('close'); 
        
    end
    
    
    if doInteractive
        [LEFT, RIGHT, UP, DOWN] = deal(28, 29, 30, 31);
        x = 1;
        y = 1;
        A = zeros(L, L);
        A(x,y) = 1;
        
        figure(1); clf;
        subplot(1,3,1);
        h_im1 = imagesc(A); axis ij square; 
        h_title = title(sprintf('(%d, %d) ', x, y));

        subplot(1,3,2);
        frmId = find( (spPhase == uPh(1)) & (spFreq == uSp(x)) & (ori == uOri(y)) );            
        thisFrame = getFrame(frmId); %#ok<FNDSB>
        % thisFrame = rescaledFrame(thisFrame);
        h_im2 = imagesc(thisFrame); colormap('gray'); axis xy square;

        subplot(1,3,3);        
        thisFrame = generateFrame(1,1,1);
        % thisFrame = rescaledFrame(thisFrame);
        h_im3 = imagesc(thisFrame); colormap('gray'); axis xy square;
        
        k = 0; 
        inputmode = 'keyboard';
        while ~isequal(k, 27)
            figure(1);
            if strcmp(inputmode, 'keyboard')
                figure(1);
                a = waitforbuttonpress;
            end
            if strcmp(inputmode, 'keyboard') && (a == 1)
                k = double(get(1, 'currentCharacter'));
                update = false;
                if isempty(k), continue; end;
                switch k
                    case UP,    if x > 1, x = x-1; update = true; end
                    case DOWN,  if x < L, x = x+1; update = true; end
                    case LEFT,  if y > 1, y = y-1; update = true; end 
                    case RIGHT, if y < L, y = y+1; update = true; end

                end
            elseif (a == 0)
                inputmode = 'mouse';
                [y1,x1] = ginput(1);
                x1 = round(x1); y1 = round(y1);
                if any([~ibetween(x1, [0, 64]), ~ibetween(y1, [0, 64])])
                    inputmode = 'keyboard';
                    continue;
                end
                if ~isequal([x1, y1], [x, y])
                    update = true;
                    x = x1;
                    y = y1;
                end
%                 inputmode = 'mouse';
                
            end
                
            
            if update
                A = zeros(L, L);
                A(x,y) = 1;
                set(h_im1, 'cdata', A);
                set(h_title, 'string', sprintf('(%d, %d) ', x, y));

                frmId = find( (spPhase == 1) & (spFreq == uSp(x)) & (ori == uOri(y)) );            
                movieFrame = getFrame(frmId); %#ok<FNDSB>
                % thisFrame = rescaledFrame(thisFrame);
                set(h_im2, 'cdata', movieFrame); 
                
                genFrame = fliplr(flipud( generateFrame(y,x,1) ));
                set(h_im3, 'cdata', genFrame); 
                
            end
            
        
        end        
        
    end
    
    
    
end



%     frms = find( (spPhase == uPh(1)) & (spFreq == uSp(64)) & (ori == uOri(64)) );
% %     spFreqvalInds = ord ( spFreq(frms) );
% %     frms = frms(spFreqvalInds);
% 
%     inds = ord ( spPhase(frms) );
%     frms = frms(inds);
%     
%     numFramesX = 4;
%     numFramesY = 2;
%     nFramesToShow = 100;
%     side_len = unique([nCols nRows]);
%     fid = fopen(filename);
%     figure(1); clf;
%     xy_i = 1;
%     gridSubPlot(numFramesY, numFramesX, [1 10]);
%     for i = 1:length(frms)
%         frmId = frms(i);
%         gridSubPlot;
% 
%         thisFrame = getFrame(frmId);% fread(fid, [side_len side_len]);        
%         imagesc(thisFrame); colormap('gray');
%         title(sprintf('(%d/%d)Frame #%d [%d, %d, %d]', i, length(frms), frmId, spPhase(frmId), spFreq(frmId), ori(frmId) ));
%         set(gca, 'xtick', [], 'ytick', []);
%         axis square
%         if ~mod(i, numFramesX *numFramesY)
%             input('');
%         end
%     end

%{
sp_pix = [128   63.0000   31.5000   21.0000   15.7500   12.6000   10.5000    9.0000    7.8750    7.0000    6.3000 ...
        5.7273    5.2500    4.8462    4.5000    4.2000    3.9375    3.7059    3.5000    3.3158    3.1500    3.0000 ...
        2.8636    2.7391    2.6250    2.5200    2.4231    2.3333    2.2500    2.1724    2.1000    2.0323    2.0000 ...
        2.0323    2.1000    2.1724    2.2500    2.3333    2.4231    2.5200    2.6250    2.7391    2.8636    3.0000 ...
        3.1500    3.3158    3.5000    3.7059    3.9375    4.2000    4.5000    4.8462    5.2500    5.7273    6.3000 ...
        7.0000    7.8750    9.0000   10.5000   12.6000   15.7500   21.0000   31.5000   63.0000];
%}