%function [pVals, H,ints] = MF_Enrichment_LA(gene_list, p_thresh,GO_matrix,names)
function [pVals, H,ints] = MF_Enrichment_LA(gene_list, p_thresh,GO_matrix,names)
import bioma.*;    import bioma.data.*;    import statistics.*
A = GO_matrix;
G = zeros( length(A) , 1);
G(gene_list,1) = 1;
all_genes = length(A);
min_group_size = 5;
% all genes, const
pVals = zeros(size(G,2),size(A,2));
ints = zeros(size(A,2),1);
for j = 1:size(G,2) % do for all custom groups
    for k = 1:size(A,2)
        a = length(intersect(find(G(:,j)),find(A(:,k))));
        b = all_genes;
        c = sum(G(:,j));
        d = sum(A(:,k));
        ints(k) = a;
        if (d < min_group_size)
            pVals(j,k) = 1;
        else
            pVals(j,k) = hygecdf(c-a-1,b,c,b-d);        
        end
    end
end
i=find(pVals < p_thresh);
H = DataMatrix(pVals(:,i),{''},names(i));
end