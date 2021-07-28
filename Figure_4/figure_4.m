%fig_4A will be the tree
fig_4b            = 0; % gumballs two clusters with bars- donor, mono-macs
fig_4c            = 0; % macs TNF
figu_4d           = 0; % macs histogram
data_for_monocole = 0;
fig_4e            = 0; % monocole states
fig_4f            = 0; % stack bar of states per bacteria
fig_4g            = 0; % modules profile
fig_4h            = 0; % stack bar modules

if fig_4b
    load ../prepare_the_data/general_info.mat
    load ../prepare_the_data/norm_ABCD.mat
    load ../prepare_the_data/filter_inf_ABCD.mat good_inf_vec_ABCD good_antibody_ABCD
    
    %PCA on all cells
    info_genes = informative_genes(ABCD_tpm,2,2, ribo_mito);
    [~,score_A] = pca(ABCD_ft(info_genes,:)');
    X_A = score_A(:,1);
    Y_A = score_A(:,2);
    [~,Xi]     = sort(X_A);     
    %diffrentiation score
    mono_gene_list = {'B2M','HLA-A','HLA-B','HLA-C'};
    macs_gene_list = {'HLA-DRA','HLA-DRB1','CD4','APOE','CD68'};
    mono_vec = ismember(gene_names, mono_gene_list);
    macs_vec = ismember(gene_names, macs_gene_list);
    mono_sum = sum(ABCD_tpm(mono_vec,:));
    macs_sum = sum(ABCD_tpm(macs_vec,:));
    mono_macs_log = log2(mono_sum./macs_sum);
    %reorder donor vec
    donor_vec = zeros(size(good_inf_vec_ABCD));
    donor_vec(good_inf_vec_ABCD == 1) = 1;
    donor_vec(good_inf_vec_ABCD == 2 | good_inf_vec_ABCD == 3) = 2;
    donor_vec(good_inf_vec_ABCD == 4 | good_inf_vec_ABCD == 5) = 3;
    donor_vec(good_inf_vec_ABCD == 6 | good_inf_vec_ABCD == 7) = 4;
   
    rand_vec = randperm(RandStream('mt19937ar','Seed',7),length(good_antibody_ABCD));  
    figure;
    subplot(10,1,1:8)
    scatter(X_A(rand_vec),Y_A(rand_vec),20,good_antibody_ABCD,'filled')
    colormap(bac_colors)
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(10,1,9);
    imagesc(mono_macs_log(Xi));
    colormap(gca,'Parula');
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    subplot(10,1,10);
    imagesc(donor_vec(Xi));
    colormap(gca,[rep_colors(1:2,:);rep_colors(5:6,:)]);
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    
    save PCA_data_ABCD.mat X_A Y_A score_A info_genes
end

if fig_4c
    load PCA_data_ABCD.mat score_A
    load ../prepare_the_data/general_info.mat gene_names
    load ../prepare_the_data/norm_ABCD.mat ABCD_tpm 
    load ../prepare_the_data/macs_ABCD.mat macs_vec_ABCD
    load cmap_color_blind.mat
    
    PC2 = score_A(:,2);
    PC3 = score_A(:,3);
    figure;
    TNF = ismember(gene_names,'TNF');
    scatter(PC2(macs_vec_ABCD),PC3(macs_vec_ABCD),20,zscore(ABCD_tpm(TNF,macs_vec_ABCD),0,2),'filled');
    colormap(cmap_color_blind); caxis([-1 1]);
    title('TNF');
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end

if fig_4d 
    load PCA_data_ABCD.mat Y_A
    load ../prepare_the_data/general_info.mat bac_colors bac_names
    load ../prepare_the_data/macs_ABCD.mat macs_vec_ABCD
    load ../prepare_the_data/filter_inf_ABCD.mat good_antibody_ABCD
    
    macs_ab  = good_antibody_ABCD(macs_vec_ABCD);
    Y_A_mac  = Y_A(macs_vec_ABCD);
    figure;
    for ab = 1:max(good_antibody_ABCD)
        subplot(8,1,ab);
        histogram(Y_A_mac(macs_ab == ab),100,'EdgeColor', bac_colors(ab,:),'FaceColor',bac_colors(ab,:)); hold on;
        xlim([-20 100]); % same axis as TNF above
        set(gca,'xtick',[]);    set(gca,'ytick',[]); box off;
       % title(bac_names(ab));
    end
    
end

if data_for_monocole
    load ../prepare_the_data/filter_inf_ABCD.mat good_singlets_ABCD
    load ../prepare_the_data/macs_ABCD.mat macs_vec_ABCD
    load ../prepare_the_data/general_info.mat ribo_mito gene_names
    load ../prepare_the_data/macs_ABCD.mat macs_tpm_ABCD
    
    macs_raw = good_singlets_ABCD(:,macs_vec_ABCD);
    dlmwrite('macs_mat.txt',macs_raw,'delimiter','\t')
    
    info_macs = informative_genes(macs_tpm_ABCD,2,2, ribo_mito);
    %manually save gene_names(info_macs)

    
end

if fig_4e
    load ../prepare_the_data/general_info.mat state_colors
    state_import      = readtable('state_vec_R_ABCD');
    state_vec         = state_import{:,1};
    coord_import_x    = readtable('coord_x_ABCD');
    coord_x           = coord_import_x{:,1};
    coord_import_y    = readtable('coord_y_ABCD');
    coord_y           = coord_import_y{:,1};
     
    % merge states 3+4+7 
    % to be consistente with fig 1, we will change the states to be:
    % 1 - non-activated macs (1)
    % 2 - mito (2)
    % 3 - antigene presentation (3+4+7)
    % 4 - apoptosis (5)
    % 5 - inflammatory (6)
    
    state_vec_new = zeros(size(state_vec));
    state_vec_new(state_vec == 1) = 1;
    state_vec_new(state_vec == 2) = 2;
    state_vec_new(state_vec == 3) = 3;
    state_vec_new(state_vec == 4) = 3;
    state_vec_new(state_vec == 5) = 4;
    state_vec_new(state_vec == 6) = 5;
    state_vec_new(state_vec == 7) = 3;
 
    figure;
    gscatter(coord_x,coord_y,state_vec_new,state_colors,'.',15);    
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    save ABCD_states.mat state_vec_new coord_x coord_y

end

if fig_4f
    load ABCD_states.mat state_vec_new
    load ../prepare_the_data/macs_ABCD.mat macs_ab_ABCD macs_donor_ABCD
    load ../prepare_the_data/general_info.mat bac_names state_colors 
  
    %reorder donor vec
    donor_macs = zeros(size(macs_donor_ABCD));
    donor_macs(macs_donor_ABCD == 1) = 1;
    donor_macs(macs_donor_ABCD == 2 | macs_donor_ABCD == 3) = 2;
    donor_macs(macs_donor_ABCD == 4 | macs_donor_ABCD == 5) = 3;
    donor_macs(macs_donor_ABCD == 6 | macs_donor_ABCD == 7) = 4;
    %count number of cells per bac/state/donor
    bac_state_donor = NaN(max(macs_ab_ABCD),max(state_vec_new),max(donor_macs));
    for ab = 1:max(macs_ab_ABCD)
        for s = 1: max(state_vec_new)
            for d = 1:max(donor_macs)
                bac_state_donor(ab,s,d) = sum((macs_ab_ABCD == ab)&(state_vec_new == s)&(donor_macs' == d));
            end
        end
    end
    perc_state_ab_d = bac_state_donor ./ sum(bac_state_donor,2);
    std_d  = std(perc_state_ab_d,[],3);
    ste_d  = std(perc_state_ab_d,[],3)./max(donor_macs);
    bac_state = sum(bac_state_donor,3); % all the donors together
    perc_bac_state = bac_state./sum(bac_state,2);
    
    figure;
    h= bar(perc_bac_state,'stacked'); hold on;
    h(1).FaceColor = state_colors(1,:);
    h(2).FaceColor = state_colors(2,:);
    h(3).FaceColor = state_colors(3,:);
    h(4).FaceColor = state_colors(4,:);
    h(5).FaceColor = state_colors(5,:);
    set(gca,'xtick',1:length(bac_names)); set(gca,'xticklabel',bac_names);
    ylim([0 1]);  xtickangle(45)
    errorbar(repmat(1:length(bac_names),max(state_vec_new),1),cumsum(perc_bac_state,2)',ste_d','.','color','k');
     
end

if fig_4g
   load ../prepare_the_data/macs_ABCD.mat macs_ft_ABCD macs_tpm_ABCD
   load  ABCD_states.mat state_vec_new
   load ../prepare_the_data/general_info.mat gene_names ribo_mito module_colors
   
   %PCA on state 5 cells
   state_5_tpm = macs_tpm_ABCD(:,state_vec_new == 5);
   state_5_ft  = macs_ft_ABCD(:,state_vec_new == 5);
   info_5 = informative_genes(state_5_tpm,3,2,ribo_mito); sum(info_5)
   [coeff,score_5] = pca(state_5_ft(info_5,:)');
   X_5 = score_5(:,1);
   Y_5 = score_5(:,2);
   %binning
   number_of_intervals = 10;
   binned_mat = NaN(length(gene_names),number_of_intervals);
   min_cord = min(X_5);
   max_cord = max(X_5)+eps; % I added eps to the last interval
   intervals_12 = linspace(min_cord, max_cord, number_of_intervals + 3); % i'm adding 3 to get 12 intervals
   %manually merge the last 3 intervals together to have more than 10 cells per bin
   intervals = [intervals_12(1:10),intervals_12(end)];
   categories = discretize(X_5, intervals);
   for i = 1:length(intervals)-1
       binned_mat(:,i) = mean(state_5_tpm(:,categories==i),2);
   end
   z_binned = zscore(binned_mat,0,2);
   dyn_mat = binned_mat(info_5,:);
   %kmeans
   bin_num = 4;
   rng (17);
   km_ABCD = kmeans(zscore(dyn_mat,0,2),bin_num);
   km_new_ABCD = NaN(size(km_ABCD));
   km_new_ABCD(km_ABCD == 1) = 1;
   km_new_ABCD(km_ABCD == 2) = 3;
   km_new_ABCD(km_ABCD == 3) = 2;
   km_new_ABCD(km_ABCD == 4) = 4;
   %genes vec per cluster
   km_genes_ABCD = zeros(length(info_5),max(km_new_ABCD));
   km_info  = zeros(size(info_5));
   km_info(info_5) = km_new_ABCD;
   for k = 1:max(km_new_ABCD)
       km_genes_ABCD(:,k) = km_info == k;
   end  
   %mean exp profile per cluster
   mean_exp_k = NaN(max(km_new_ABCD),size(binned_mat,2));
   for k = 1:max(km_new_ABCD)
       k_mat = zscore(binned_mat(logical(km_genes_ABCD(:,k)),:),0,2);
       mean_exp_k(k,:) = mean(k_mat);
   end
   
   figure;
   for k = 1:max(km_new_ABCD)
       plot(mean_exp_k(k,:)','color',module_colors(k,:),'LineWidth',4); hold on;
       set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
       xlim([1 10])
   end
   
   save ABCD_state_5.mat state_5_tpm state_5_ft X_5 Y_5 info_5 number_of_intervals intervals categories binned_mat z_binned
   save ABCD_modules.mat dyn_mat km_new_ABCD km_genes_ABCD 
       
end

if fig_4h
    load ABCD_state_5.mat state_5_tpm
    load ABCD_modules.mat km_new_ABCD km_genes_ABCD
    load ../prepare_the_data/macs_ABCD.mat macs_ab_ABCD macs_donor_ABCD
    load ABCD_states.mat state_vec_new
    load ../prepare_the_data/general_info.mat module_colors bac_names
    
    %reorder donor vec
    donor_macs = zeros(size(macs_donor_ABCD));
    donor_macs(macs_donor_ABCD == 1) = 1;
    donor_macs(macs_donor_ABCD == 2 | macs_donor_ABCD == 3) = 2;
    donor_macs(macs_donor_ABCD == 4 | macs_donor_ABCD == 5) = 3;
    donor_macs(macs_donor_ABCD == 6 | macs_donor_ABCD == 7) = 4;
    state_5_d  = donor_macs(state_vec_new == 5);
    state_5_ab = macs_ab_ABCD(state_vec_new == 5);
    %mean expression per module
    z_state_5 = zscore(state_5_tpm,0,2);
    z_mean_exp_k = NaN(size(km_genes_ABCD,2),size(z_state_5,2));
    for k = 1:max(km_new_ABCD)
        z_mean_exp_k(k,:) = mean(z_state_5(logical(km_genes_ABCD(:,k)),:));
    end
    %normalize modules
    min_mean = z_mean_exp_k-min(z_mean_exp_k);
    norm_mean = min_mean./sum(min_mean);
    norm_mean(isnan(norm_mean)) = 0;
    [~,max_norm] = max(norm_mean);
    %count bac/module/donor
    bac_km_donor = NaN(max(state_5_ab)-1,max(max_norm),max(state_5_d)); %do not include unexposed
    for ab = 1:max(state_5_ab)-1
        for s = 1: max(max_norm)
            for d = 1:max(state_5_d)
                bac_km_donor(ab,s,d) = sum((state_5_ab' == ab)&(max_norm == s)&(state_5_d == d));
            end
        end
    end
    perc_state_ab_d = bac_km_donor ./ sum(bac_km_donor,2);
    std_d  = std(perc_state_ab_d,[],3);
    ste_d  = std(perc_state_ab_d,[],3)./max(state_5_d);
    bac_state = sum(bac_km_donor,3); % all the donors together
    perc_bac_state = bac_state./sum(bac_state,2);
    
    figure;
    h= bar(perc_bac_state,'stacked'); hold on;
    h(1).FaceColor = module_colors(1,:);
    h(2).FaceColor = module_colors(2,:);
    h(3).FaceColor = module_colors(3,:);
    h(4).FaceColor = module_colors(4,:);
    set(gca,'xtick',1:length(bac_names)-1); set(gca,'xticklabel',bac_names);
    ylim([0 1]);  xtickangle(45)
    errorbar(repmat(1:length(bac_names)-1,max(max_norm),1),cumsum(perc_bac_state,2)',ste_d','.','color','k');
    
end

%sup figures:
fig_s4a  = 0; %FACS 50%
fig_s4b  = 0; %pie chart bac/donor
fig_s4c  = 0; %two clusters
fig_s4d  = 0; %GBS/all bac states
fig_s4e  = 0; %heatmap GO terms
fig_s4f  = 0; %state_5 gradient
fig_s4g  = 0; %GBS/all modules

if fig_s4a
   load ../prepare_the_data/general_info.mat bac_colors bac_names rep_colors
   precent_inf = xlsread('precent_inf.xlsx');

   d_names = {'donor A','donor B','donor C','donor D'};
   d_colors = [rep_colors(1:2,:);rep_colors(5:6,:)];

   figure;
   for d = 1:length(d_names)
       plot(1:length(bac_names),precent_inf(:,d),'.','color',d_colors(d,:),'MarkerSize',20); hold on;
   end
   xlim([0.5 8.5]); ylim([0 100]);
   ylabel('precent infection')
   set(gca,'xtick',1:length(bac_names)); set(gca,'xticklabel',bac_names); xtickangle(45)
   legend(d_names)
    
end

if fig_s4b
   load ../prepare_the_data/filter_inf_ABCD.mat good_inf_vec_ABCD good_antibody_ABCD
   load ../prepare_the_data/general_info.mat bac_colors bac_names
   
   donor_vec = zeros(size(good_inf_vec_ABCD));
   donor_vec(good_inf_vec_ABCD == 1) = 1;
   donor_vec(good_inf_vec_ABCD == 2 | good_inf_vec_ABCD == 3) = 2;
   donor_vec(good_inf_vec_ABCD == 4 | good_inf_vec_ABCD == 5) = 3;
   donor_vec(good_inf_vec_ABCD == 6 | good_inf_vec_ABCD == 7) = 4;
   
   sum_donor_bac_mat = NaN(max(good_antibody_ABCD),max(donor_vec));
   for ab = 1:max(good_antibody_ABCD)
       for d = 1:max(donor_vec)
           sum_donor_bac_mat(ab,d) = sum(good_antibody_ABCD' == ab & donor_vec == d);
       end
   end
   prop_bac_donor = sum_donor_bac_mat./sum(sum_donor_bac_mat);
   d_names = {'donor A','donor B','donor C','donor D'};

   figure;
   for d = 1:max(donor_vec)
      subplot(1,4,d);
      pie(prop_bac_donor(:,d));
      colormap(bac_colors);
      title(d_names(d));
   end
   
end

if fig_s4c
   load ../prepare_the_data/norm_ABCD.mat ABCD_ft ABCD_tpm
   load ../prepare_the_data/macs_ABCD.mat cluster_id_ABCD
   load ../prepare_the_data/general_info.mat ribo_mito
   
   info_ABCD = informative_genes(ABCD_tpm,2,2,ribo_mito); sum(info_ABCD)
   [~,score_A] = pca(ABCD_ft(info_ABCD,:)');
   X_A = score_A(:,1);
   Y_A = score_A(:,2);

   two_colors =[0.9769 0.9839 0.0805; 0.2422 0.1504 0.6603]; 

   figure;
   subplot(8,1,1:3)
   [~,~,s_order] = dendrogram(linkage(pdist(ABCD_ft(info_ABCD,:)'),'ward'),0);
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
   subplot(8,1,4);
   imagesc(cluster_id_ABCD(s_order));
   colormap(two_colors)
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
   subplot(8,1,5:8);
   scatter(X_A,Y_A,20,cluster_id_ABCD,'filled');
   colormap(two_colors);
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
 
    
end

if fig_s4d
   load ABCD_states.mat state_vec_new
   load ../prepare_the_data/macs_ABCD.mat macs_tpm_ABCD
   load ../prepare_the_data/general_info.mat ribo_mito gene_names
   load ../Figure_2/states_diff_genes.mat genes_up_P
   names_no_ribo_mito = gene_names(~ribo_mito);
   no_ribo_mito_tpm = macs_tpm_ABCD(~ribo_mito,:);
   
   %[p_ABCD, ~, ~, ~, ~, ind_p_ABCD] = diff_exp_genes(names_no_ribo_mito,state_vec_new, no_ribo_mito_tpm);
   %load diff_exp_genes_ABCD.mat p_ABCD ind_p_ABCD
   Pval_thersh = 0.005;
   P_genes_ABCD = p_ABCD < Pval_thersh; sum(P_genes_ABCD)
   % take top X genes
   gene_count = 100;
   genes_up_ABCD   = ind_p_ABCD(1:gene_count,:);
   genes_up_P_ABCD = genes_up_ABCD.*(P_genes_ABCD(genes_up_ABCD));
   
   %compare GBS and all bac state genes
   state_names = {'state 1','state 2','state 3','state 4','state 5'};
   GBS_state_list = NaN(length(names_no_ribo_mito),length(state_names));
   for s = 1:length(state_names)
       list_temp = ismember(1:length(names_no_ribo_mito),genes_up_P(:,s));
       GBS_state_list(:,s) = list_temp;
   end
   bac_state_list = NaN(length(names_no_ribo_mito),length(state_names));
   for s = 1:length(state_names)
       list_temp = ismember(1:length(names_no_ribo_mito),genes_up_P_ABCD(:,s));
       bac_state_list(:,s) = list_temp;
   end
   pvals = NaN(length(state_names),length(state_names));
   for n = 1:length(state_names)
       pvals(:,n) = Enrichment(logical(bac_state_list(:,n)),0.05,GBS_state_list,state_names);
   end
   
   figure;
   imagesc(-log10(pvals))
   
    save diff_exp_genes_ABCD.mat p_ABCD ind_p_ABCD
end

if fig_s4e
    load diff_exp_genes_ABCD.mat
    load ABCD_states.mat state_vec_new
    load GO_human_mat_10X.mat GO_mat_human_50_10X GO_names_50_10X
    load ../prepare_the_data/general_info.mat ribo_mito
    
    Pval_thersh = 0.005;
    P_genes_ABCD = p_ABCD < Pval_thersh; sum(P_genes_ABCD)
    gene_count = 100;
    genes_up_ABCD   = ind_p_ABCD(1:gene_count,:);
    genes_up_P_ABCD = genes_up_ABCD.*(P_genes_ABCD(genes_up_ABCD));
    GO_no_ribo_mito = GO_mat_human_50_10X(~ribo_mito,:);

    p_thresh_go = 0.000001; % this p_thresh is not relevant at all because it's in the function
    Pvals = NaN(length(GO_names_50_10X),max(state_vec_new));
    for c = 1:max(state_vec_new)
        genes_to_take = genes_up_P_ABCD(:,c);
        genes_to_take(genes_to_take == 0) = [];
        Pvals(:,c) = Enrichment(genes_to_take, p_thresh_go, GO_no_ribo_mito, GO_names_50_10X);
    end
    
    number_of_GOs = 5;
    p_ind    = NaN(number_of_GOs,max(state_vec_new));
    for c = 1:max(state_vec_new)
        [~,xi] = sort(Pvals(:,c));
        p_ind(:,c) = xi(1:number_of_GOs);
    end
    p_ind_join = p_ind(:);
    p_val_ind  = Pvals(p_ind_join,:);
    state_names = {'state 1','state 2','state 3','state 4','state 5'};
    
    figure;
    imagesc(-log10(p_val_ind)); 
    set(gca,'ytick',1:size(p_ind_join,1)); set(gca,'yticklabel',GO_names_50_10X(p_ind_join));
    set(gca,'xtick',1:length(state_names)); set(gca,'xticklabel',state_names);
    colorbar; caxis([0 10]);
    
end

if fig_s4f
   load ABCD_state_5.mat state_5_tpm X_5 Y_5 intervals
   load ../prepare_the_data/general_info.mat gene_names
   load cmap_color_blind.mat
    
   figure;
   gene = ismember(gene_names,'IL1B');
   scatter(X_5,Y_5,20,zscore(state_5_tpm(gene,:),0,2),'filled'); hold on;
   colormap(cmap_color_blind); caxis([-1 1]);
   for i = 1:length(intervals)
      line([intervals(i) intervals(i)],[-40 30],'Color','k','LineStyle','--') 
   end
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;

    
end

if fig_s4g
   load ABCD_modules.mat km_genes_ABCD
   load ../Figure_3/GBS_modules.mat km_genes
   
   clust_names = {'cluster 1','cluster 2','cluster 3','cluster 4'};
   pvals = NaN(length(clust_names),length(clust_names));
   for n = 1:length(clust_names)
       pvals(:,n) = Enrichment(logical(km_genes_ABCD(:,n)),0.00005,km_genes,clust_names);
   end
   figure;
   imagesc(-log10(pvals))
   caxis([0 80])
   
end
