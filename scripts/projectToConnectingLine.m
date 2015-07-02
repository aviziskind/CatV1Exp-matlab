function [x1_s_nrm, x2_s_nrm, x1_proj, x2_proj, pca_comps] = projectToConnectingLine(x1, x2, usePCAforProj, pca_comps_prev)
    %%
    x1_orig = x1;
    x2_orig = x2;
    %%
    
    if usePCAforProj
%         x1 = x1_orig;
%         x2 = x2_orig;
%         mX1 = mean(x1, 1); mX2 = mean(x2,1);   
%         x1_ms = bsxfun(@minus, x1, mX1);
%         x2_ms = bsxfun(@minus, x2, mX1);
%         x1_ms = bsxfun(@minus, x1, mX1);
%         x2_ms = bsxfun(@minus, x2, mX1);
        
        [coeff, pca_comps, ~, ~, eig_vals] = doPCA([x1;x2], 2);
        if ~isempty(pca_comps_prev)
            cc1 = pearsonR(pca_comps(:,1), pca_comps_prev(:,1));
            if cc1 < 0
                coeff(:,1) = -coeff(:,1);
            end
            cc2 = pearsonR(pca_comps(:,2), pca_comps_prev(:,2));
            if cc2 < 0
                coeff(:,2) = -coeff(:,2);
            end            
            
        end
        
        n1 = size(x1,1); n2 = size(x2,1);
%         x1_s = coeff(1:n1);
%         x2_s = coeff(n1+1:n1+n2);
        f = 1/sqrt(eig_vals(1));
        coeff = coeff * f;

        x1_proj = coeff(:,1:n1)';
        x2_proj = coeff(:,n1+1:n1+n2)';
        
        x1_s_nrm = x1_proj(:,1);
        x2_s_nrm = x2_proj(:,1);
        
        if mean(x1_s_nrm) > 0
            x1_proj = -coeff(:,1:n1)';
            x2_proj = -coeff(:,n1+1:n1+n2)';

            x1_s_nrm = x1_proj(:,1);
            x2_s_nrm = x2_proj(:,1);            
        end
        3;
    else
        %%
        x1 = x1_orig;
        x2 = x2_orig;
        if size(x1, 2) > 16
            x1 = x1';
        end
        if size(x2, 2) > 16
            x2 = x2';
        end
        
        mX1 = mean(x1, 1); mX2 = mean(x2,1);   
        v = mX2 - mX1;
        w = v(:)/(norm(v));
        D = norm(v);
        
        
        x1 = bsxfun(@minus, x1, mX1)/D;
        x2 = bsxfun(@minus, x2, mX1)/D;

%         M = w*w';    
%         
%         x1_proj1 = (M * x1'); 
%         x2_proj1 = (M * x2');
%         
%         x1_s = w'*x1_proj1;
%         x2_s = w'*x2_proj1;
        
        W = GramSchmidt(w);

%         x1_proj1 = (w' * x1')'; 
%         x2_proj1 = (w' * x2')';        
        
        x1_proj = (W' * x1')'; 
        x2_proj = (W' * x2')';
        
%         x1_s = w'*x1_proj1;
%         x2_s = w'*x2_proj1;
        
        x1_s_nrm = x1_proj(:,1); 
        x2_s_nrm = x2_proj(:,1); 
        
%         rescale = 1;
%         if rescale
%             
%             
%         end
%         i
        
    end
    
    show = 0;
    if show
        figure(57); clf;
        hist2({x1_s_nrm, x2_s_nrm}, 50, 'line');    

        %% 
        figure(55);  clf; hold on;
        plot(x1(:,1), x1(:,2), 'bo');
        plot(x2(:,1), x2(:,2), 'ro');
        plot(mean(x1(:,1),1), mean(x1(:,2),1), 'ko');
        plot(mean(x2(:,1),1), mean(x2(:,2),1), 'ks');
    
        %%
        figure(56); clf; hold on;
        plot(x1_proj(:,1), x1_proj(:,2), 'b.');
        plot(x2_proj(:,1), x2_proj(:,2), 'r.');

        

%         figure(58); clf;
%         hist2({x1_s2, x2_s2}, 50, 'stacked');
    end
end
