function nspk_mtx = edbGetMovieNspkPerFrame(hnd, Gid, Cellnmb)
  nspk_mtx = [];
	Did  = edbFindDID(hnd,Gid);
	Sid  = edbFindSid(hnd,Gid);
	Mids = edbFindMid(hnd,Gid)';
	[tb_tick, te_tick] = edbGetTbTe(hnd, 'mb', Did);
  if isnan(tb_tick(1)) | isempty(tb_tick)
    return;
  end
	spk = edbGetSpikes(hnd, Gid);
	if isempty(spk)
       return;
	end
	cspk = edbGetCellSpikes(spk,Cellnmb);
	if isempty(cspk)
       return;
	end
	sync = edbGetSyncs(hnd, Sid);
	if isempty(sync)
       return;
	end

	t_spk = cspk(:,1);
	t_syn = sync(:,1);

  hdrMids = [];
	uMids = unique(Mids);
	for mid = uMids
      inds = find(Mids == mid);
      Npres = length(inds);
      tb = tb_tick(inds);
      te = te_tick(inds);
      for pnmb = 1:Npres
        inds_spk = find(t_spk >= tb(pnmb) & t_spk < te(pnmb));
        inds_syn = find(t_syn >= tb(pnmb) & t_syn < te(pnmb));
        if isempty(inds_spk)
          nspk = zeros(size(inds_syn));
        else
%           length(t_syn(inds_syn))  
          nspk = histc(t_spk(inds_spk), t_syn(inds_syn));
        end
        nspk = nspk(:);
        if ~isempty(nspk_mtx) & (size(nspk_mtx,1) ~= size(nspk,1))
          nspk_mtx = [];
          return;
        end
        nspk_mtx = [nspk_mtx,nspk];
        hdrMids = [hdrMids,mid];
      end
	end

  nspk_mtx = [hdrMids;nspk_mtx];

  