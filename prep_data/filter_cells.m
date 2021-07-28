function [good_cells,ribo_mito] = filter_cells(raw_exp_matrix, gene_names, min_umi,ribo_fraction,mito_fraction)
%filter cells with low number of umi, and high ribo and mito fraction  
%the expression matrix will be the raw matrix befor normalization
%gene_names is a list of gene names in order as in the matrix

    cond_1     = sum(raw_exp_matrix)>min_umi; 
    ribo_genes = startsWith(gene_names,'RPL')|startsWith(gene_names,'RPS');
    mito_genes = startsWith(gene_names,'MT-')|startsWith(gene_names,'MTRNR');
    cond_2     = sum(raw_exp_matrix(ribo_genes,:))./sum(raw_exp_matrix) < ribo_fraction;
    cond_3     = sum(raw_exp_matrix(mito_genes,:))./sum(raw_exp_matrix) < mito_fraction;
    good_cells = (cond_1 & cond_2 & cond_3); 
    ribo_mito  = ribo_genes|mito_genes;
    
end

