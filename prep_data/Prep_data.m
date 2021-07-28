load_ABCD          = 0;
normalize_ABCD     = 0; %general info is saved here
macs_ABCD          = 0;
organize_GBS_unexp = 0;
GBS_macs           = 0;
load_mut           = 0;
normalize_mut      = 0;
macs_mut           = 0;

if load_ABCD
    A   = csvread('filtered_feature_bc_matrix_inf_A.csv', 1, 1);
    B1  = csvread('filtered_feature_bc_matrix_inf_B1.csv', 1, 1);
    B2  = csvread('filtered_feature_bc_matrix_inf_B2.csv', 1, 1);
    C1  = csvread('filtered_feature_bc_matrix_inf_C1.csv', 1, 1);
    C2  = csvread('filtered_feature_bc_matrix_inf_C2.csv', 1, 1);
    D1  = csvread('filtered_feature_bc_matrix_inf_D1.csv', 1, 1);
    D2  = csvread('filtered_feature_bc_matrix_inf_D2.csv', 1, 1);
    
    all_inf = [A,B1,B2,C1,C2,D1,D2];
    inf_vec = repelem(1:7,[size(A,2),size(B1,2),size(B2,2),size(C1,2),size(C2,2),size(D1,2),size(D2,2)]);
    
    features = readtable('features'); % all features are the same between samples
    % there are double names in features(:,2) in most of them only one is expressed:
    proteins = readtable('protein_coding_genes_GRCh38.p12.txt');
    proteins_vec = ismember(features{:,2}, proteins{:,1});
    n=(features{:,2});
    [ii, jj, kk] = unique(n);
    double_names = ii(accumarray(kk, 1) > 1);
    lia= ismember(n,double_names);
    figure; imagesc(log10(1+A(lia,:))) % I do it just on A, but it's the same for all them
    set(gca,'ytick',1:sum(lia)); set(gca,'yticklabel',n(lia));
    
    % I'll remove the doublets, the one that has lower sum exp
    sc_no_dup = NaN(size(n));
    for g = 1:length(n)
        if lia(g) < 1
            sc_no_dup(g,1) = 1;
        else
            ind_tmp = strmatch(n(g),n,'exact');
            ind_sum = sum(A(ind_tmp,:),2);
            cond = ind_sum(1)>ind_sum(2);
            if cond >0 %means that the first has higher exp
                sc_no_dup(ind_tmp(1),1) = 1;
                sc_no_dup(ind_tmp(2),1) = 0;
            else
                sc_no_dup(ind_tmp(1),1) = 0;
                sc_no_dup(ind_tmp(2),1) = 1;
            end
        end
    end
    % sc_no_dup is a vec without the duplicates (zeros)
    protein_no_dup = proteins_vec&logical(sc_no_dup); % i will use the same vec in the other datasets as well
    gene_names = n(protein_no_dup);
    
    cell_barcodes_A  = readtable('barcodes_inf_A');
    cell_barcodes_B1 = readtable('barcodes_inf_B1');
    cell_barcodes_B2 = readtable('barcodes_inf_B2');
    cell_barcodes_C1 = readtable('barcodes_inf_C1');
    cell_barcodes_C2 = readtable('barcodes_inf_C2');
    cell_barcodes_D1 = readtable('barcodes_inf_D1');
    cell_barcodes_D2 = readtable('barcodes_inf_D2');
    
    antibody_barcode_A  = readtable('infA_hashing_1_4_2');
    antibody_barcode_B1 = readtable('infB1_hashing_1_4_2');
    antibody_barcode_B2 = readtable('infB2_hashing_1_4_2');
    antibody_barcode_C1 = readtable('infC1_hashing_1_4_2');
    antibody_barcode_C2 = readtable('infC2_hashing_1_4_2');
    antibody_barcode_D1 = readtable('infD1_hashing_1_4_2');
    antibody_barcode_D2 = readtable('infD2_hashing_1_4_2');
    
    barcodes_vec_A  = ismember(cell_barcodes_A{:,1},antibody_barcode_A{:,1});
    barcodes_vec_B1 = ismember(cell_barcodes_B1{:,1},antibody_barcode_B1{:,1});
    barcodes_vec_B2 = ismember(cell_barcodes_B2{:,1},antibody_barcode_B2{:,1});
    barcodes_vec_C1 = ismember(cell_barcodes_C1{:,1},antibody_barcode_C1{:,1});
    barcodes_vec_C2 = ismember(cell_barcodes_C2{:,1},antibody_barcode_C2{:,1});
    barcodes_vec_D1 = ismember(cell_barcodes_D1{:,1},antibody_barcode_D1{:,1});
    barcodes_vec_D2 = ismember(cell_barcodes_D2{:,1},antibody_barcode_D2{:,1});
    
    barcodes_vec = [barcodes_vec_A ; barcodes_vec_B1 ; barcodes_vec_B2 ; ...
        barcodes_vec_C1 ; barcodes_vec_C2 ; barcodes_vec_D1 ; barcodes_vec_D2];
    %some antibodies barcodes are not present in cell barcodes
    antibody_vec_A  = ismember(antibody_barcode_A{:,1},cell_barcodes_A{:,1});
    antibody_vec_B1 = ismember(antibody_barcode_B1{:,1},cell_barcodes_B1{:,1});
    antibody_vec_B2 = ismember(antibody_barcode_B2{:,1},cell_barcodes_B2{:,1});
    antibody_vec_C1 = ismember(antibody_barcode_C1{:,1},cell_barcodes_C1{:,1});
    antibody_vec_C2 = ismember(antibody_barcode_C2{:,1},cell_barcodes_C2{:,1});
    antibody_vec_D1 = ismember(antibody_barcode_D1{:,1},cell_barcodes_D1{:,1});
    antibody_vec_D2 = ismember(antibody_barcode_D2{:,1},cell_barcodes_D2{:,1});
    
    ab_A  = antibody_barcode_A{:,2};
    ab_B1 = antibody_barcode_B1{:,2};
    ab_B2 = antibody_barcode_B2{:,2};
    ab_C1 = antibody_barcode_C1{:,2};
    ab_C2 = antibody_barcode_C2{:,2};
    ab_D1 = antibody_barcode_D1{:,2};
    ab_D2 = antibody_barcode_D2{:,2};
    
    antibody_ABCD = [ab_A(antibody_vec_A);ab_B1(antibody_vec_B1);ab_B2(antibody_vec_B2);...
        ab_C1(antibody_vec_C1);ab_C2(antibody_vec_C2);ab_D1(antibody_vec_D1);ab_D2(antibody_vec_D2);];
    
    singlets_ABCD = all_inf(protein_no_dup,barcodes_vec);
    inf_vec_ABCD = inf_vec(barcodes_vec);
    
    save('inf_ABCD_data.mat', 'singlets_ABCD', 'gene_names', 'antibody_ABCD', 'inf_vec_ABCD', 'protein_no_dup', '-v7.3')
end

if normalize_ABCD
    load inf_ABCD_data.mat
    
    [good_cells, ribo_mito] = filter_cells(singlets_ABCD,gene_names,1000,0.3,0.15);
    
    good_singlets_ABCD = singlets_ABCD(:,good_cells);
    good_inf_vec_ABCD  = inf_vec_ABCD(good_cells);
    good_antibody_ABCD = antibody_ABCD(good_cells);
    
    ABCD_tpm = median(sum(good_singlets_ABCD))*bsxfun(@rdivide, good_singlets_ABCD, sum(good_singlets_ABCD,1));
    ABCD_ft  = sqrt(ABCD_tpm)+sqrt(ABCD_tpm+1);
    
    bac_colors    = [0.1300 0.6700 0.8000;0.2000 0.6275 0.1725;0.8902 0.1020 0.1098;1.0000 0.4980   0;0.5333 0.2549 0.6157;1.0000 0.8510 0.1843;0.9059 0.5412 0.7647;0.5 0.5 0.5];
    rep_colors    = [0.2549 0.6706 0.3647;1.0000 0.8510 0.1843;  0 0.4275 0.1725;0.9373 0.2314 0.1725;0.6471 0.0588 0.0824;0.1137 0.5686 0.7529;0.1451 0.2039 0.5804];
    state_colors  = [0.6510 0.8078 0.8902 ; 0.1216 0.4706 0.7059 ; 0.6980 0.8745 0.5412 ; 0.8902 0.1020 0.1098 ; 0.9922 0.7490 0.4353];
    module_colors = [0.6196 0.0039 0.2588 ; 0.9569 0.4275 0.2627 ; 0.4000 0.7608 0.6471 ; 0.3686 0.3098 0.6353];
    bac_names     = {'Salmonella','Yersinia','Shigella','Enterococcus','Staph','Listeria','GBS','Unexposed'};
    
    save('filter_inf_ABCD.mat', 'good_singlets_ABCD','good_inf_vec_ABCD', 'good_antibody_ABCD', '-v7.3')
    save('norm_ABCD.mat', 'ABCD_tpm','ABCD_ft','-v7.3')
    save general_info.mat ribo_mito gene_names bac_colors rep_colors bac_names protein_no_dup state_colors module_colors

end

if macs_ABCD 
    load norm_ABCD.mat
    load general_info.mat ribo_mito
    load filter_inf_ABCD.mat good_singlets_ABCD good_antibody_ABCD good_inf_vec_ABCD
     
    info_ABCD = informative_genes(ABCD_tpm,2,2,ribo_mito); sum(info_ABCD)
    cluster_id_ABCD = cluster(linkage(pdist(ABCD_ft(info_ABCD,:)'),'ward'),'maxclust',2)';
    
    macs_vec_ABCD = false(size(cluster_id_ABCD));
    sum_cluster_1 = mean(sum(good_singlets_ABCD(:,cluster_id_ABCD == 1)));
    sum_cluster_2 = mean(sum(good_singlets_ABCD(:,cluster_id_ABCD == 2)));
    if sum_cluster_1 > sum_cluster_2 %if cluster 1 is macs
        macs_vec_ABCD(cluster_id_ABCD == 1) = true; 
    else 
        macs_vec_ABCD(cluster_id_ABCD ~= 1) = true;
    end
    
    % macs_vec equals 1 are macs
    macs_tpm_ABCD = ABCD_tpm(:,macs_vec_ABCD);
    macs_ft_ABCD  = ABCD_ft(:,macs_vec_ABCD);
    macs_ab_ABCD  = good_antibody_ABCD(macs_vec_ABCD);
    macs_donor_ABCD = good_inf_vec_ABCD(macs_vec_ABCD);
    
    save macs_ABCD.mat cluster_id_ABCD macs_tpm_ABCD macs_ft_ABCD macs_ab_ABCD macs_donor_ABCD macs_vec_ABCD

end

if organize_GBS_unexp
    load general_info.mat
    load filter_inf_ABCD.mat
    
    GBS_mat     = csvread('filtered_GBS.csv', 1, 1);
    GBS_protein = GBS_mat(protein_no_dup,:);
    [good_cells, ribo_mito] = filter_cells(GBS_protein,gene_names,1000,0.3,0.15);
    good_GBS = GBS_protein(:,good_cells);

    GBS = 7;
    unexp = 8;
    GBS_unexp_cells = good_antibody_ABCD == GBS | good_antibody_ABCD == unexp; 
    %join mats
    GBS_unexp_mat = [good_GBS,good_singlets_ABCD(:,GBS_unexp_cells)];
    GBS_unexp_tpm = median(sum(GBS_unexp_mat))*bsxfun(@rdivide, GBS_unexp_mat, sum(GBS_unexp_mat,1));
    GBS_unexp_ft  = sqrt(GBS_unexp_tpm)+sqrt(GBS_unexp_tpm+1);
    
    GBS_unexp_vec = [ones(1,sum(good_cells)),good_antibody_ABCD(GBS_unexp_cells)'];
    GBS_unexp_vec(GBS_unexp_vec == 7) = 1; % GBS is one 
    GBS_unexp_vec(GBS_unexp_vec == 8) = 2; % unexp is 2

    % reaction vec
    reaction_vec = [zeros(1,sum(good_cells)),good_inf_vec_ABCD(GBS_unexp_cells)];
    reaction_vec = reaction_vec +1; % GBS only is 1 and the rest is now 2-8
    
    save filtered_GBS_unexp.mat GBS_unexp_mat GBS_unexp_tpm GBS_unexp_ft GBS_unexp_vec reaction_vec
    
end

if GBS_macs
    load filtered_GBS_unexp.mat
    load general_info.mat

    info_GBS_unexp = informative_genes(GBS_unexp_tpm,2,2,ribo_mito); sum(info_GBS_unexp)
    cluster_id_GBS_unexp = cluster(linkage(pdist(GBS_unexp_ft(info_GBS_unexp,:)'),'ward'),'maxclust',2)';
    
    macs_vec = false(size(cluster_id_GBS_unexp));
    sum_cluster_1 = mean(sum(GBS_unexp_mat(:,cluster_id_GBS_unexp == 1)));
    sum_cluster_2 = mean(sum(GBS_unexp_mat(:,cluster_id_GBS_unexp == 2)));
    if sum_cluster_1 > sum_cluster_2 %if cluster 1 is macs
        macs_vec(cluster_id_GBS_unexp == 1) = true; 
    else 
        macs_vec(cluster_id_GBS_unexp ~= 1) = true;
    end
    
    % macs_vec equals 1 are macs
    macs_tpm_GBS_unexp   = GBS_unexp_tpm(:,macs_vec);
    macs_ft_GBS_unexp    = GBS_unexp_ft(:,macs_vec);
    macs_GBS_unexp_vec   = GBS_unexp_vec(macs_vec); % the names are confusing!
    macs_reaction_vec    = reaction_vec(macs_vec);
    
    save GBS_macs.mat macs_tpm_GBS_unexp macs_ft_GBS_unexp macs_GBS_unexp_vec macs_reaction_vec cluster_id_GBS_unexp macs_vec
end

if load_mut
    load general_info.mat

    mut_A = csvread('filtered_feature_bc_matrix_mut_A.csv', 1, 1);
    mut_B = csvread('filtered_feature_bc_matrix_mut_B.csv', 1, 1);
    
    all_mut = [mut_A,mut_B];
    mut_vec = repelem(1:2,[size(mut_A,2),size(mut_B,2)]);
    
    cell_barcodes_mut_A = readtable('barcodes_mut_A');
    cell_barcodes_mut_B = readtable('barcodes_mut_B');
    
    antibody_barcode_mut_A = readtable('mut_A_hashing');
    antibody_barcode_mut_B = readtable('mut_B_hashing');

    barcodes_vec_mut_A = ismember(cell_barcodes_mut_A{:,1},antibody_barcode_mut_A{:,1});
    barcodes_vec_mut_B = ismember(cell_barcodes_mut_B{:,1},antibody_barcode_mut_B{:,1});

    barcodes_vec = [barcodes_vec_mut_A ; barcodes_vec_mut_B];
    %some antibodies barcodes are not present in cell barcodes
    antibody_vec_mut_A = ismember(antibody_barcode_mut_A{:,1},cell_barcodes_mut_A{:,1});
    antibody_vec_mut_B = ismember(antibody_barcode_mut_B{:,1},cell_barcodes_mut_B{:,1});
    
    ab_mut_A  = antibody_barcode_mut_A{:,2};
    ab_mut_B = antibody_barcode_mut_B{:,2};
    
    antibody_mut = [ab_mut_A(antibody_vec_mut_A);ab_mut_B(antibody_vec_mut_B)];
    
    singlets_mat = all_mut(protein_no_dup,barcodes_vec);
    mut_vec_singlets = mut_vec(barcodes_vec);
    
    save('mut_AB_data.mat', 'singlets_mat', 'antibody_mut', 'mut_vec_singlets', '-v7.3')

end

if normalize_mut
    load mut_AB_data.mat
    load general_info.mat
    [good_cells, ribo_mito] = filter_cells(singlets_mat,gene_names,1000,0.3,0.15);
    
    good_singlets_mut = singlets_mat(:,good_cells);
    good_inf_vec_mut  = mut_vec_singlets(good_cells);
    good_antibody_mut = antibody_mut(good_cells);
    % order the mut to be WT cpsE cylE unexp
    re_antibody_mut = ones(size(good_antibody_mut));
    re_antibody_mut(good_antibody_mut == 1) = 1;
    re_antibody_mut(good_antibody_mut == 2) = 4;
    re_antibody_mut(good_antibody_mut == 3) = 2;
    re_antibody_mut(good_antibody_mut == 4) = 3;

    mut_tpm = median(sum(good_singlets_mut))*bsxfun(@rdivide, good_singlets_mut, sum(good_singlets_mut,1));
    mut_ft  = sqrt(mut_tpm)+sqrt(mut_tpm+1);
    
    mut_colors = [0.1300 0.6700 0.8000;0.0314 0.1882 0.4196;0.1373 0.5451 0.2706;0.5000 0.5000 0.5000];
    mut_names = {'WT','cpsE','cylE','Unexposed'};
    
    save('filter_mut.mat', 'good_singlets_mut','good_inf_vec_mut', 're_antibody_mut', '-v7.3')
    save('norm_mut.mat', 'mut_tpm','mut_ft','mut_colors','mut_names','-v7.3')

end

if macs_mut
    load norm_mut.mat
    load general_info.mat ribo_mito
    load filter_mut.mat good_singlets_mut good_inf_vec_mut re_antibody_mut
    
    info_mut = informative_genes(mut_tpm,2,2,ribo_mito); sum(info_mut)
    cluster_id_mut = cluster(linkage(pdist(mut_ft(info_mut,:)'),'ward'),'maxclust',3)';
    
    % in this data set, there are three clusters, cluster 1 is the mono
    macs_vec_mut = false(size(cluster_id_mut));
    sum_cluster_1 = mean(sum(good_singlets_mut(:,cluster_id_mut == 1)));
    sum_cluster_2 = mean(sum(good_singlets_mut(:,cluster_id_mut == 2)));
    sum_cluster_3 = mean(sum(good_singlets_mut(:,cluster_id_mut == 3)));
    macs_vec_mut(cluster_id_mut ~= 1) = true;

    % macs_vec equals 1 are macs
    macs_tpm_mut = mut_tpm(:,macs_vec_mut);
    macs_ft_mut  = mut_ft(:,macs_vec_mut);
    macs_ab_mut  = re_antibody_mut(macs_vec_mut);
    macs_donor_mut = good_inf_vec_mut(macs_vec_mut);
    
    save macs_mut.mat cluster_id_mut macs_vec_mut macs_tpm_mut macs_ft_mut macs_ab_mut macs_donor_mut
    
end

if load HIJKLM
    load general_info.mat

    H  = csvread('filtered_feature_bc_matrix_inf_H.csv', 1, 1);
    I  = csvread('filtered_feature_bc_matrix_inf_I.csv', 1, 1);
    J  = csvread('filtered_feature_bc_matrix_inf_J.csv', 1, 1);
    K  = csvread('filtered_feature_bc_matrix_inf_K.csv', 1, 1);
    L  = csvread('filtered_feature_bc_matrix_inf_L.csv', 1, 1);
    M  = csvread('filtered_feature_bc_matrix_inf_M.csv', 1, 1);
    
    inf_HIJKLM = [H,I,J,K,L,M];
    inf_vec = repelem(1:6,[size(H,2),size(I,2),size(J,2),size(K,2),size(L,2),size(M,2)]);
  
    cell_barcodes_H  = readtable('barcodes_inf_H');
    cell_barcodes_I = readtable('barcodes_inf_I');
    cell_barcodes_J = readtable('barcodes_J');
    cell_barcodes_K = readtable('barcodes_K');
    cell_barcodes_L = readtable('barcodes_L');
    cell_barcodes_M = readtable('barcodes_M');
    
    antibody_barcode_H = readtable('inf_H_hashing');
    antibody_barcode_I = readtable('inf_I_hashing');
    antibody_barcode_J = readtable('inf_J_hashing');
    antibody_barcode_K = readtable('inf_K_hashing');
    antibody_barcode_L = readtable('inf_L_hashing');
    antibody_barcode_M = readtable('inf_M_hashing');
    
    barcodes_vec_H = ismember(cell_barcodes_H{:,1},antibody_barcode_H{:,1});
    barcodes_vec_I = ismember(cell_barcodes_I{:,1},antibody_barcode_I{:,1});
    barcodes_vec_J = ismember(cell_barcodes_J{:,1},antibody_barcode_J{:,1});
    barcodes_vec_K = ismember(cell_barcodes_K{:,1},antibody_barcode_K{:,1});
    barcodes_vec_L = ismember(cell_barcodes_L{:,1},antibody_barcode_L{:,1});
    barcodes_vec_M = ismember(cell_barcodes_M{:,1},antibody_barcode_M{:,1});
    
    barcodes_vec = [barcodes_vec_H ; barcodes_vec_I ; barcodes_vec_J ; ...
                        barcodes_vec_K ; barcodes_vec_L ; barcodes_vec_M];
    %some antibodies barcodes are not present in cell barcodes
    antibody_vec_H = ismember(antibody_barcode_H{:,1},cell_barcodes_H{:,1});
    antibody_vec_I = ismember(antibody_barcode_I{:,1},cell_barcodes_I{:,1});
    antibody_vec_J = ismember(antibody_barcode_J{:,1},cell_barcodes_J{:,1});
    antibody_vec_K = ismember(antibody_barcode_K{:,1},cell_barcodes_K{:,1});
    antibody_vec_L = ismember(antibody_barcode_L{:,1},cell_barcodes_L{:,1});
    antibody_vec_M = ismember(antibody_barcode_M{:,1},cell_barcodes_M{:,1});
    
    ab_H = antibody_barcode_H{:,2};
    ab_I = antibody_barcode_I{:,2};
    ab_J = antibody_barcode_J{:,2};
    ab_K = antibody_barcode_K{:,2};
    ab_L = antibody_barcode_L{:,2};
    ab_M = antibody_barcode_M{:,2};
    
    antibody_HIJKLM = [ab_H(antibody_vec_H);ab_I(antibody_vec_I);ab_J(antibody_vec_J);...
                          ab_K(antibody_vec_K);ab_L(antibody_vec_L);ab_M(antibody_vec_M)];
    
    singlets_HIJKLM = inf_HIJKLM(protein_no_dup,barcodes_vec);
    inf_vec_HIJKLM = inf_vec(barcodes_vec);
    
    save('inf_HIJKLM_data.mat', 'singlets_HIJKLM', 'antibody_HIJKLM', 'inf_vec_HIJKLM', '-v7.3')

end

if normalize_HIJKLM
    load inf_HIJKLM_data.mat
    load general_info.mat
    
    [good_cells, ribo_mito] = filter_cells(singlets_HIJKLM,gene_names,1000,0.3,0.15);
    
    good_singlets_HIJKLM = singlets_HIJKLM(:,good_cells);
    good_inf_vec_HIJKLM  = inf_vec_HIJKLM(good_cells);
    good_antibody_HIJKLM = antibody_HIJKLM(good_cells);
    
    HIJKLM_tpm = median(sum(good_singlets_HIJKLM))*bsxfun(@rdivide, good_singlets_HIJKLM, sum(good_singlets_HIJKLM,1));
    HIJKLM_ft  = sqrt(HIJKLM_tpm)+sqrt(HIJKLM_tpm+1);
    
    two_four_h_vec = ones(size(good_inf_vec_HIJKLM)); 
    two_four_h_vec(good_inf_vec_HIJKLM == 2 | good_inf_vec_HIJKLM == 5 |good_inf_vec_HIJKLM == 6) = 2; % 1 is 2h, 2 is 4h  
    two_four_h_colors = [0.9882 0.5529 0.3843; 0.5529 0.6275 0.7961];
    
    save('filter_inf_HIJKLM.mat', 'good_singlets_HIJKLM','good_inf_vec_HIJKLM', 'good_antibody_HIJKLM', '-v7.3')
    save('norm_HIJKLM.mat', 'HIJKLM_tpm','HIJKLM_ft', 'two_four_h_vec','two_four_h_colors','-v7.3')

end

if macs_HIJKLM
    load norm_HIJKLM.mat
    load general_info.mat ribo_mito
    load filter_inf_HIJKLM.mat good_singlets_HIJKLM good_antibody_HIJKLM good_inf_vec_HIJKLM
     
    info_HIJKLM = informative_genes(HIJKLM_tpm,2,2,ribo_mito); sum(info_HIJKLM)
    cluster_id_HIJKLM = cluster(linkage(pdist(HIJKLM_ft(info_HIJKLM,:)'),'ward'),'maxclust',2)';
    
    macs_vec_HIJKLM = false(size(cluster_id_HIJKLM));
    sum_cluster_1 = mean(sum(good_singlets_HIJKLM(:,cluster_id_HIJKLM == 1)));
    sum_cluster_2 = mean(sum(good_singlets_HIJKLM(:,cluster_id_HIJKLM == 2)));
    if sum_cluster_1 > sum_cluster_2 %if cluster 1 is macs
        macs_vec_HIJKLM(cluster_id_HIJKLM == 1) = true; 
    else 
        macs_vec_HIJKLM(cluster_id_HIJKLM ~= 1) = true;
    end
    
    % macs_vec equals 1 are macs
    macs_tpm_HIJKLM = HIJKLM_tpm(:,macs_vec_HIJKLM);
    macs_ft_HIJKLM  = HIJKLM_ft(:,macs_vec_HIJKLM);
    macs_ab_HIJKLM  = good_antibody_HIJKLM(macs_vec_HIJKLM);
    macs_donor_HIJKLM = good_inf_vec_HIJKLM(macs_vec_HIJKLM);
    
    save macs_HIJKLM.mat cluster_id_HIJKLM macs_vec_HIJKLM macs_tpm_HIJKLM macs_ft_HIJKLM macs_ab_HIJKLM macs_donor_HIJKLM

end

