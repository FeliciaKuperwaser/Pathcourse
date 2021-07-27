function [info_genes] = informative_genes(tpm_mat,mean_thresh,fano_factor_thresh, ribo_mito_vec)
% This function generates a vector of 1/0 of the informative genes
% the tpm_mat is your tpm matrix
% the mean and fano thresh are the treshold for the mean and fano factor
% ribo_mito_vec is a vector of 1/0 in the size of all the genes in your matrix
max_read_per_cell = max(tpm_mat,[],2);
m =max_read_per_cell>=mean_thresh*mean(max_read_per_cell);
fano_factor = var(tpm_mat,[],2)./mean(tpm_mat,2);
f = fano_factor>fano_factor_thresh*nanmean(fano_factor);
info_genes = f&m&(~ribo_mito_vec);
end

