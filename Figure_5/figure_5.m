% fig 5a will be a schematic
fig_5b  = 0; % PCA macs colored by 2/4h and unexposed with a bar stack of 2/4h of exposed
fig_5c  = 0; % stackbar 2/4 h with error bar
fig_5d  = 0; % stakebar modules over 2/4h
%fig 5e will be the model

if fig_5b
    load ../prepare_the_data/norm_HIJKLM.mat two_four_h_colors two_four_h_vec HIJKLM_ft HIJKLM_tpm
    load ../prepare_the_data/macs_HIJKLM.mat macs_vec_HIJKLM macs_ab_HIJKLM
    load ../prepare_the_data/general_info.mat gene_names ribo_mito
    
    [info_2_4]     = informative_genes(HIJKLM_tpm,2,2, ribo_mito);
    [~,score_2_4]  = pca(HIJKLM_ft(info_2_4,:)');
    PC2            = score_2_4(:,2);
    PC3            = score_2_4(:,3);
    PC2_macs       = PC2(macs_vec_HIJKLM);
    PC3_macs       = PC3(macs_vec_HIJKLM);
    macs_2_4_h_vec = two_four_h_vec(macs_vec_HIJKLM);
    %generate vec for two/four/unexp --> 1/2/3
    unexp_two_four_vec = macs_2_4_h_vec;
    unexp_two_four_vec(macs_ab_HIJKLM == 8) = 3;
    rand_vec = randperm(RandStream('mt19937ar','Seed',62),length(unexp_two_four_vec));
    unexp_two_four_colors = [two_four_h_colors; [0.5 0.5 0.5]];
    % bin to ten bins by PC1, without unexp and calc proportions of 2 and 4h in each bin
    unexp_vec = macs_ab_HIJKLM == 8;
    X_m_exp = PC2_macs(~unexp_vec);
    Y_m_exp = PC3_macs(~unexp_vec);
    number_of_intervals = 10;
    min_cord = min(X_m_exp);
    max_cord = max(X_m_exp)+eps; % I added eps for the mat_interval
    intervals = linspace(min_cord, max_cord, number_of_intervals + 1); % 10 bins have 11 edges
    categories = discretize(X_m_exp, intervals);
    bin_2_4 = NaN(max(macs_2_4_h_vec),number_of_intervals);
    for i = 1:number_of_intervals
        bin_2_4(1,i) = sum(categories' == i & macs_2_4_h_vec(~unexp_vec) == 1);
        bin_2_4(2,i) = sum(categories' == i & macs_2_4_h_vec(~unexp_vec) == 2);
    end
    bin_prop = bin_2_4 ./ sum(bin_2_4);
    time_names = {'two h','four h'};
    
    figure;
    subplot(10,1,1:8)
    scatter(PC2_macs(rand_vec),PC3_macs(rand_vec),20,unexp_two_four_vec(rand_vec),'filled');
    colormap(unexp_two_four_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(10,1,9:10);
    h= bar(bin_prop','stacked');
    h(1).FaceColor = two_four_h_colors(1,:);
    h(2).FaceColor = two_four_h_colors(2,:);
    legend(time_names);
    
    
end

if fig_5c
    load ../prepare_the_data/macs_ABCD.mat macs_tpm_ABCD macs_ft_ABCD
    load ../prepare_the_data/general_info.mat ribo_mito state_colors bac_names
    load ../Figure_4/ABCD_states.mat state_vec_new
    load ../prepare_the_data/macs_HIJKLM.mat macs_ft_HIJKLM macs_vec_HIJKLM macs_ab_HIJKLM macs_donor_HIJKLM
    load ../prepare_the_data/norm_HIJKLM.mat two_four_h_vec
    
    %infer state for 2 and 4 hours samples based on ABCD
    [info_macs] = informative_genes(macs_tpm_ABCD,2,2, ribo_mito);
    m = fitcknn(macs_ft_ABCD(info_macs,:)',state_vec_new);
    c = predict(m,macs_ft_HIJKLM(info_macs,:)');
    two_four_vec_macs = two_four_h_vec(macs_vec_HIJKLM);
    %reorder donor vec
    donor_macs = zeros(size(macs_donor_HIJKLM));
    donor_macs(macs_donor_HIJKLM == 1 | macs_donor_HIJKLM == 2) = 1;
    donor_macs(macs_donor_HIJKLM == 3 | macs_donor_HIJKLM == 5) = 2;
    donor_macs(macs_donor_HIJKLM == 4 | macs_donor_HIJKLM == 6) = 3;
    %separete the 2 and 4 hours to to mats
    two_h = 1; four_h = 2;
    ab_2h = macs_ab_HIJKLM(two_four_vec_macs == two_h);
    d_2h  = donor_macs(two_four_vec_macs == two_h);
    s_2h  = c(two_four_vec_macs == two_h);
    ab_4h = macs_ab_HIJKLM(two_four_vec_macs == four_h);
    d_4h  = donor_macs(two_four_vec_macs == four_h);
    s_4h  = c(two_four_vec_macs == four_h);
    %count number of cells per bac/state/donor 2h
    bac_state_donor_2h = NaN(max(ab_2h),max(s_2h),max(d_2h));
    for ab = 1:max(ab_2h)
        for s = 1: max(s_2h)
            for d = 1:max(d_2h)
                bac_state_donor_2h(ab,s,d) = sum((ab_2h == ab)&(s_2h == s)&(d_2h' == d));
            end
        end
    end
    perc_state_ab_d_2h = bac_state_donor_2h ./ sum(bac_state_donor_2h,2);
    std_d_2h  = std(perc_state_ab_d_2h,[],3);
    ste_d_2h  = std(perc_state_ab_d_2h,[],3)./max(d_2h);
    bac_state_2h = sum(bac_state_donor_2h,3); % all the donors together
    perc_bac_state_2h = bac_state_2h./sum(bac_state_2h,2);
    %count number of cells per bac/state/donor 4h
    bac_state_donor_4h = NaN(max(ab_4h),max(s_4h),max(d_4h));
    for ab = 1:max(ab_4h)
        for s = 1: max(s_4h)
            for d = 1:max(d_4h)
                bac_state_donor_4h(ab,s,d) = sum((ab_4h == ab)&(s_4h == s)&(d_4h' == d));
            end
        end
    end
    perc_state_ab_d_4h = bac_state_donor_4h ./ sum(bac_state_donor_4h,2);
    std_d_4h  = std(perc_state_ab_d_4h,[],3);
    ste_d_4h  = std(perc_state_ab_d_4h,[],3)./max(d_4h);
    bac_state_4h = sum(bac_state_donor_4h,3); % all the donors together
    perc_bac_state_4h = bac_state_4h./sum(bac_state_4h,2);
    %merge the mats
    perc_bac_state_2_4 = [perc_bac_state_2h ; perc_bac_state_4h];
    ste_d_2_4          = [ste_d_2h ; ste_d_4h];
    bac_vec            = repmat(1:length(bac_names),1,2)';
    
    figure;
    for ab = 1:max(macs_ab_HIJKLM)
        subplot(2,4,ab)
        prop_ab_temp = perc_bac_state_2_4(bac_vec == ab,:); hold on;
        h= bar(prop_ab_temp,'stacked');
        h(1).FaceColor = state_colors(1,:);
        h(2).FaceColor = state_colors(2,:);
        h(3).FaceColor = state_colors(3,:);
        h(4).FaceColor = state_colors(4,:);
        h(5).FaceColor = state_colors(5,:);
        title(bac_names(ab)); ylim([0 1]);
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        errorbar(repmat(1:sum(bac_vec == ab),max(state_vec_new),1),cumsum(perc_bac_state_2_4(bac_vec == ab,:),2)',ste_d_2_4(bac_vec == ab,:)','.','color','k');
    end
    
    %precent state 5 in 2/4h
    perc_5_2h = perc_state_ab_d_2h(1:7,5,:); 
    perc_5_4h = perc_state_ab_d_4h(1:7,5,:); 
    [~,pval] = ttest2(perc_5_2h(:),perc_5_4h(:));
    save two_four_states.mat c m two_four_vec_macs donor_macs
end

if fig_5d
    load ../prepare_the_data/macs_HIJKLM.mat macs_tpm_HIJKLM macs_ft_HIJKLM macs_vec_HIJKLM
    load ../prepare_the_data/general_info.mat ribo_mito gene_names
    load two_four_states.mat c two_four_vec_macs donor_macs
    load ../prepare_the_data/norm_HIJKLM.mat two_four_h_colors
    
    %PCA on state 5 cells
    state_5_two_four = two_four_vec_macs(c == 5);
    state_5_donor = donor_macs(c == 5);
    state_5_tpm = macs_tpm_HIJKLM(:,c == 5);
    state_5_ft  = macs_ft_HIJKLM(:,c == 5);
    info_5 = informative_genes(state_5_tpm,2,2,ribo_mito); sum(info_5)
    [coeff,score_5] = pca(state_5_ft(info_5,:)');
    X_5 = score_5(:,1);
    Y_5 = score_5(:,2);
    %binning
    number_of_intervals = 10;
    binned_mat = NaN(length(gene_names),number_of_intervals);
    min_cord = min(X_5);
    max_cord = max(X_5)+eps; % I added eps to the last interval
    intervals = linspace(min_cord, max_cord, number_of_intervals + 1);
    categories = discretize(X_5, intervals);
    for i = 1:length(intervals)-1
        binned_mat(:,i) = mean(state_5_tpm(:,categories==i),2);
    end
    z_binned = zscore(binned_mat,0,2);
    dyn_mat = binned_mat(info_5,:);
    %kmeans
    bin_num = 4;
    rng (17);
    km_HIJKLM = kmeans(zscore(dyn_mat,0,2),bin_num);
    %         figure;
    %         for k = 1:max(km_HIJKLM)
    %             subplot(1,4,k)
    %             k_mat = zscore(dyn_mat(km_HIJKLM == k,:),0,2);
    %             mean_k = mean(k_mat);
    %             plot(1:number_of_intervals,mean_k)
    %         end
    km_2_4 = zeros(size(km_HIJKLM));
    km_2_4(km_HIJKLM == 1) = 4;
    km_2_4(km_HIJKLM == 2) = 1;
    km_2_4(km_HIJKLM == 3) = 2;
    km_2_4(km_HIJKLM == 4) = 3;
    %genes vec per cluster
    km_genes_2_4 = zeros(length(info_5),max(km_2_4));
    km_info  = zeros(size(info_5));
    km_info(info_5) = km_2_4;
    for k = 1:max(km_2_4)
        km_genes_2_4(:,k) = km_info == k;
    end
    %mean expression per module
    z_state_5 = zscore(state_5_tpm,0,2);
    z_mean_exp_k = NaN(size(km_genes_2_4,2),size(z_state_5,2));
    for k = 1:max(km_2_4)
        z_mean_exp_k(k,:) = mean(z_state_5(logical(km_genes_2_4(:,k)),:));
    end
    %normalize modules
    min_mean = z_mean_exp_k-min(z_mean_exp_k);
    norm_mean = min_mean./sum(min_mean);
    norm_mean(isnan(norm_mean)) = 0;
    [~,max_norm] = max(norm_mean);
    %count 2_4/module/donor
    two_four_km_donor = NaN(max(state_5_two_four),max(max_norm),max(state_5_donor)); %do not include unexposed
    for t = 1:max(state_5_two_four)
        for s = 1: max(max_norm)
            for d = 1:max(state_5_donor)
                two_four_km_donor(t,s,d) = sum((state_5_two_four == t)&(max_norm == s)&(state_5_donor == d));
            end
        end
    end
    perc_state_two_four_d = two_four_km_donor ./ sum(two_four_km_donor,2);
    std_d  = std(perc_state_two_four_d,[],3);
    ste_d  = std(perc_state_two_four_d,[],3)./max(state_5_donor);
    bac_state = sum(two_four_km_donor,3); % all the donors together
    perc_two_four_state = bac_state./sum(bac_state);
    
    time_names = {'two h','four h'};
    
    figure;
    h= bar(perc_two_four_state','stacked'); hold on;
    h(1).FaceColor = two_four_h_colors(1,:);
    h(2).FaceColor = two_four_h_colors(2,:);
    ylim([0 1]);
    errorbar(repmat(1:max(max_norm),length(time_names),1),cumsum(perc_two_four_state),ste_d,'.','color','k');
    
end

%sup figures:
fig_s5a  = 0; %two culsters colored by 2/4h
fig_s5bc  = 0; %poodle colored by state and 2/4

if fig_s5a
    load ../prepare_the_data/norm_HIJKLM.mat HIJKLM_ft HIJKLM_tpm two_four_h_vec two_four_h_colors
    load ../prepare_the_data/general_info.mat ribo_mito gene_names rep_colors
    load ../prepare_the_data/filter_inf_HIJKLM.mat good_antibody_HIJKLM good_inf_vec_HIJKLM
    
    info_genes = informative_genes(HIJKLM_tpm,2,2,ribo_mito); sum(info_genes)
    [~,score_A] = pca(HIJKLM_ft(info_genes,:)');
    X_A    = score_A(:,1);
    Y_A    = score_A(:,2);
    [~,Xi] = sort(X_A);
      
    donor_vec = zeros(size(good_inf_vec_HIJKLM));
    donor_vec(good_inf_vec_HIJKLM == 1 | good_inf_vec_HIJKLM == 2) = 1;
    donor_vec(good_inf_vec_HIJKLM == 3 | good_inf_vec_HIJKLM == 5) = 2;
    donor_vec(good_inf_vec_HIJKLM == 4 | good_inf_vec_HIJKLM == 6) = 3;
  
    %mono/macs ratio
    mono_gene_list = {'B2M','HLA-A','HLA-B','HLA-C'};
    macs_gene_list = {'HLA-DRA','HLA-DRB1','CD4','APOE','CD68'};
    mono_vec = ismember(gene_names, mono_gene_list);
    macs_vec = ismember(gene_names, macs_gene_list);
    mono_sum = sum(HIJKLM_tpm(mono_vec,:));
    macs_sum = sum(HIJKLM_tpm(macs_vec,:));
    mono_macs_log = log2(mono_sum./macs_sum);
    
    rand_vec = randperm(RandStream('mt19937ar','Seed',7),length(two_four_h_vec));  
    figure;
    subplot(10,1,1:8)
    scatter(X_A(rand_vec),Y_A(rand_vec),20,two_four_h_vec(rand_vec),'filled')
    colormap(gca,two_four_h_colors)
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(10,1,9);
    imagesc(mono_macs_log(Xi));
    colormap(gca,'Parula');
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    subplot(10,1,10);
    imagesc(donor_vec(Xi));
    %colormap(gca,[rep_colors(1,:);rep_colors(3,:);rep_colors(5,:)]);
    colormap(gca,rep_colors);
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    
end

if fig_s4bc
    load two_four_states.mat c two_four_vec_macs
    load ../prepare_the_data/general_info.mat state_colors ribo_mito gene_names
    load ../prepare_the_data/norm_HIJKLM.mat two_four_h_colors 
    load ../prepare_the_data/macs_ABCD.mat macs_tpm_ABCD macs_ft_ABCD
    load ../prepare_the_data/macs_HIJKLM.mat macs_ft_HIJKLM
    
    coord_import_x = readtable('../Figure_4/coord_x_ABCD');
    coord_x        = coord_import_x{:,1};
    coord_import_y = readtable('../Figure_4/coord_y_ABCD');
    coord_y        = coord_import_y{:,1};
    
    [info_macs] = informative_genes(macs_tpm_ABCD,2,2, ribo_mito);
    m_x = fitcknn(macs_ft_ABCD(info_macs,:)',coord_x);
    c_x = predict(m_x,macs_ft_HIJKLM(info_macs,:)');
    m_y = fitcknn(macs_ft_ABCD(info_macs,:)',coord_y);
    c_y = predict(m_y,macs_ft_HIJKLM(info_macs,:)');
    rand_vec = randperm(RandStream('mt19937ar','Seed',62),length(two_four_vec_macs));

    figure;
    subplot(1,2,1);
    plot(coord_x,coord_y,'.','MarkerSize',20,'color',[0.75 0.75 0.75]); hold on;
    scatter(c_x,c_y,20,c,'filled');
    colormap(gca,state_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(1,2,2);
    plot(coord_x,coord_y,'.','MarkerSize',20,'color',[0.75 0.75 0.75]); hold on;
    scatter(c_x(rand_vec),c_y(rand_vec),20,two_four_vec_macs(rand_vec),'filled');
    colormap(gca,two_four_h_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end

