function [gene_order] = zavit(M,shift_order)

M = zscore(M')'; 
[~, score] = pca(M); 
X = score(:,1); Y = score(:,2); 
zavit = mod(atan2d(X,Y)+shift_order, 360) ;
[~, gene_order] = sort(zavit); %gene order is the zavit ordering
