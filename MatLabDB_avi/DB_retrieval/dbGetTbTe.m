function [tb_tick, te_tick] = dbGetTbTe(hnd, Did)
    
    stimTableName = getDatabaseTableForDid(Did);

    [tb_tick, te_tick] = getFieldsFromDatabaseTable(hnd, ...
        {'LNG_START_TICK', 'LNG_END_TICK'}, stimTableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');

%     if any( [isempty(tb_tick) isempty(te_tick)])
%         syncs = dbGetSyncs('Did', Did, 'tick');
%         [nPreBlankFrms, nPostBlankFrms, nSustainedFrms] = getFieldsFromDatabaseTable(hnd, ...
%             {'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM', 'LNG_N_SUSTAINED_FRM'}, stimTableName, {'DATAFILE_ID', Did});
%         nPres = length(nSustainedFrms);
%         tb_tick2 = zeros(1,nPres);
%         te_tick2 = zeros(1,nPres);
%         for pres_i = 1:nPres
%             presBaseTickId = sum(nPreBlankFrms(1:pres_i-1)) + sum(nSustainedFrms(1:pres_i-1)) + sum(nPostBlankFrms(1:pres_i-1)) ;
%             tb_tick2Id = presBaseTickId + nPreBlankFrms(pres_i);
%             te_tick2Id = presBaseTickId + nPreBlankFrms(pres_i) + nSustainedFrms(pres_i);
%             tb_tick2(pres_i) = syncs( tb_tick2Id );            
%             te_tick2(pres_i) = syncs( te_tick2Id );
%         end
% %         assert(length(syncs) == nPreBlankFrms+nSustainedFrms+nPostBlankFrms);        
% %     end
%     3;
    
end

%     if any( [isempty(tb_tick) isempty(te_tick)])
%         syncs = dbGetSyncs('Did', Did, 'tick');
%         [nPreBlankFrms, nPostBlankFrms, nSustainedFrms] = getFieldsFromDatabaseTable(hnd, ...
%             {'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM', 'LNG_N_SUSTAINED_FRM'}, stimTableName, {'DATAFILE_ID', Did});
%         tb_tick = syncs(nPreBlankFrms);
%         te_tick = syncs(nPreBlankFrms+nSustainedFrms);
%         assert(length(syncs) == nPreBlankFrms+nSustainedFrms+nPostBlankFrms);        
%     end
