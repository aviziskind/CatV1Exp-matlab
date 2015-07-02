% function fig_equations
% correlation coefficient:
cc_s = '$$r = \frac{(x-\bar{x})^T(y-\bar{y})}{||x-\bar{x}||\ ||y-\bar{y}||} = \frac{\sum_i(x_i-\bar{x})\sum_i(y_i-\bar{y}) }{ \sqrt{\sum_{i}(x_i-\bar{x})^2} \sqrt{\sum_{i}(y_i-\bar{y})^2} }$$' ;

% rho (spearman's rho)
rho_s = '$$\rho = \frac{(R-\bar{R})^T(S-\bar{S})}{||R-\bar{R}||\ ||S-\bar{S}||} = \frac{\sum_i(R_i-\bar{R})\sum_i(S_i-\bar{S}) }{ \sqrt{\sum_{i}(R_i-\bar{R})^2} \sqrt{\sum_{i}(S_i-\bar{S})^2} }$$' ;
rho_exp_s = '$$x_i \rightarrow R_i, \quad y_i \rightarrow S_i$$';

% n choose 2
nchoose2 = '$$\left(\!\! \begin{array}{c}n\\2\end{array} \!\!\right)$$';



str = cc_s;
figure(10); clf;
text('interpreter', 'latex', 'position', [.05,.5], 'string', str, 'fontsize', 20 );
axis(axis);
set(gca, 'position', get(gca, 'outerposition'));
set(gca, 'visible', 'off')