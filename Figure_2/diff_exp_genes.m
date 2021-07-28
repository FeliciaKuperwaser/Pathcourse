function [p_values, eff_vec, sign_genes, eff_genes, eff_sign_genes, ind_p_values] = diff_exp_genes(gene_names,cluster_vec, tpm_mat)
% This function calculate differentially expressed genes between each cluster to the other clusters
%   gene_names is a vectors of the gene names in the size of your matrix
%   cluster_vec is a vector in the size of the number of samples with cluster number for each cell
%   tpm_mat is the tpm expression matrix
    
    p_values   = NaN(length(gene_names),max(cluster_vec));
    eff_vec    = NaN(length(gene_names),max(cluster_vec));
    sign_genes = NaN(length(gene_names),max(cluster_vec));
    
    warning('off');
    for g = 1:length(gene_names)
        for c = 1:max(cluster_vec)
            p_values(g,c)   = ranksum(tpm_mat(g,cluster_vec==c),tpm_mat(g,~(cluster_vec==c)));
            eff_temp        = mes(tpm_mat(g,cluster_vec==c),tpm_mat(g,~(cluster_vec==c)),'hedgesg');
            eff_vec(g,c)    = eff_temp.hedgesg;
            sign_genes(g,c) = (mean(tpm_mat(g,cluster_vec==c))-mean(tpm_mat(g,~(cluster_vec==c))))>0;
        end
    end
    warning('on');
    
    eff_genes = abs(eff_vec)>0.2;
    eff_sign_genes = eff_genes&sign_genes;
    p_values(isnan(p_values))=10; % to rm the NaN to have a bad p value
    ind_p_values = NaN(length(gene_names),max(cluster_vec));
    
    for c = 1:max(cluster_vec)
        [~, ind_p_values(:,c)] = sort(p_values(:,c) + ~eff_sign_genes(:,c));
    end
    
    
end

