function [cdata, curClims] = getSpikeDensityCData(spkDensity, pairId, useCells, densityType, densityScale, curClims, cmapSclFactor_cell, cmapSclFactor_grp, clustNSpikes, cellsGrouping, curCellZoom)
    nCells = length(useCells);
    densFieldName = [iff(strcmp(densityType, 'spike density'), 'abs', 'norm') '_' densityScale(1:3)];

    nBins = size(spkDensity.abs_lin{1}, 1);    
    mergeClustersOfSameCell = 0;

%     showAll_groupedByCluster = ~exist('cellsGrouping', 'var') || isempty(cellGrouping);

    zoomOnSingleCell = exist('curCellZoom', 'var') && ~isempty(curCellZoom);

    showAll_groupedByCell =   exist('cellsGrouping', 'var') && ~isempty(cellsGrouping) && ~ischar(cellsGrouping) && ~zoomOnSingleCell;        
    
    switch densityType
        case 'spike density',
            cdata = sum(spkDensity.abs_lin{pairId}(:,:,useCells), 3 );
            if strncmp(densityScale, 'log', 3)
                idx = cdata > 0;
                cdata(idx) = log10(cdata(idx));
                cdata(~idx) = nan;
            end
            curClims(1) = min([curClims(1); cdata(:)]);
            curClims(2) = max([curClims(2); cdata(:)]);

        case 'cluster density',            
                            
            mergeClustersOfSameCell_here = mergeClustersOfSameCell && showAll_groupedByCell;
            
            if mergeClustersOfSameCell_here                                  
                nCells = length(cellsGrouping);
                
                curSpkDs = zeros(nBins, nBins, nCells);
                for ci = 1:nCells                    
                    curClusts = cellsGrouping{ci};
                    curClusts = curClusts(useCells(curClusts));
                    n = length(curClusts);
                    
                    % add histograms of sub-clusters, weighted by # spikes.
                    props = clustNSpikes(curClusts); props = props/sum(props);
                    curSpkDs_cell = spkDensity.norm_lin{pairId}(:,:,curClusts);
                    curSpkDs_cell = sum(bsxfun(@times, curSpkDs_cell, reshape(props, [1,1,n])), 3);
                    
                    if strcmp(densityScale, 'log')
                        curSpkDs_cell = logNZ(1+curSpkDs_cell);
                        curSpkDs_cell = curSpkDs_cell / max(curSpkDs_cell(:));                        
                    end                    
                    curSpkDs(:,:,ci) = curSpkDs_cell /(cmapSclFactor_cell);
                    
                end                    
                                
            else  % just retrieve 
                curSpkDs = spkDensity.(densFieldName){pairId}(:,:,useCells)/(cmapSclFactor_cell*1.2);
                                
            end
            
            if isempty(curSpkDs)
                cdata = zeros(nBins, nBins, 1);
                return;
            end

%             curSpkDs = spkDensity.(densFieldName){pairId}(:,:,useCells)/(cmapSclFactor_cell);
            
            % find which cluster has the highest density.
            [maxDens, cell_used_idx] = max(curSpkDs, [], 3);
            if ~mergeClustersOfSameCell_here
                
                if zoomOnSingleCell
                    useCellsIdx = find(useCells(cellsGrouping{curCellZoom}));
                else
                    useCellsIdx = find(useCells);
                end                
                cellsIdx = useCellsIdx(cell_used_idx);                
                
            end
            
            if showAll_groupedByCell && ~mergeClustersOfSameCell_here
                    
                    clustIdx = cellsIdx;
                    for cell_i = 1:length(cellsGrouping)
                        clustIds_thisCell = cellsGrouping{cell_i}(:)';
                        for clust_id = clustIds_thisCell
                            cellsIdx(clustIdx == clust_id) = cell_i;
                        end
                    end
                    
            end
                
            
            cdata = maxDens + cellsIdx-1; 
            
            cdata(maxDens==0) = nan; % give color of background (black) and not just of the first cell.
            cdata = cdata * (1 - cmapSclFactor_grp) + cmapSclFactor_grp*nCells;
            assert(all(ibetween(cdata(~isnan(cdata)), 0, nCells)));
    end
    
end

function X = logNZ(X)
    idx_nz = X > 0;
    X(idx_nz) = log10(X(idx_nz));    
end

%         case 'cluster density',
%             nm = ['abs_' densityScale(1:3)];
%             % [log], then [normalize], then rescale, then take max, then final rescale
%             curSpkDs = spkDensity.(densFieldName){pairId}(:,:,useCells)/(cmapSclFactor_cell);
%             %                                 if strncmp(densityScale, 'log', 3)
%             %                                     idx = curSpkDs > 0;
%             %                                     curSpkDs = log10(curSpkDs);
%             %                                     curSpkDs(~idx) = 0;
%             %                                 end
%             %                                 for i = 1:nnz(showCells)
%             %                                     curSpkDs(:,:,i) = curSpkDs(:,:,i) / (maxElement(curSpkDs(:,:,i))*(cmapSclFactor_cell)) ;
%             %                                 end
%             %                                 curSpkDs = curSpkDs / cmapSclFactor_cell;
%             [maxDens, cell_used_idx] = max(curSpkDs, [], 3);
%             cdata = maxDens + (useCells(cell_used_idx)-1);
%             cdata(maxDens==0) = nan; % and just not the first cell
%             cdata = cdata * (1 - cmapSclFactor_grp) + cmapSclFactor_grp*nCells;
%             assert(all(ibetween(cdata(~isnan(cdata)), 0, nCells)));
