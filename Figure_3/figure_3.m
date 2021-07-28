fig_3a     = 1; %kmeans gene modules
fig_3b     = 0; %GO terms km 
fig_3c     = 0; %PCA colored by modules
fig_3d     = 0; %stack bar modules in sc
fig_3e     = 0; %corr sc on modules
fig_3f     = 0; %boxplot pdist

if fig_3a
   load ../prepare_the_data/GBS_macs.mat macs_ft_GBS_unexp macs_tpm_GBS_unexp
   load ../figure_2/GBS_unexp_states.mat state_vec_new
   load ../prepare_the_data/general_info.mat gene_names ribo_mito module_colors
   
   %PCA on state 5 cells
   state_5_tpm = macs_tpm_GBS_unexp(:,state_vec_new == 5);
   state_5_ft  = macs_ft_GBS_unexp(:,state_vec_new == 5);
   info_5 = informative_genes(state_5_tpm,3,2,ribo_mito); sum(info_5)
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
   rng (27);
   km_GBS = kmeans(zscore(dyn_mat,0,2),bin_num);
   km_ordered = zeros(size(km_GBS));
   km_ordered(km_GBS == 1) = 4;
   km_ordered(km_GBS == 2) = 1;
   km_ordered(km_GBS == 3) = 2;
   km_ordered(km_GBS == 4) = 3;
   %genes vec per cluster
   km_genes = zeros(length(info_5),max(km_ordered));
   km_info  = zeros(size(info_5));
   km_info(info_5) = km_ordered;
   for k = 1:max(km_ordered)
       km_genes(:,k) = km_info == k;
   end
   
   figure;
   for k = 1:max(km_ordered)
       subplot(1,4,k)
       k_mat = zscore(binned_mat(logical(km_genes(:,k)),:),0,2);
       mean_exp_k = mean(k_mat);
       plot(k_mat','color',[0.9 0.9 0.9]); hold on;
       plot(mean_exp_k','color',module_colors(k,:),'LineWidth',4);
       set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
       xlim([0.5 8.5])
   end
   
   save GBS_state_5.mat state_5_tpm state_5_ft X_5 Y_5 info_5 number_of_intervals intervals categories binned_mat z_binned
   save GBS_modules.mat dyn_mat km_ordered km_genes
end

if fig_3b %empty
    
end

if fig_3c
    load GBS_state_5.mat X_5 Y_5 state_5_tpm info_5 
    load GBS_modules.mat km_ordered km_genes
    load cmap_color_blind.mat
    
    %mean exp genes per cluster
    z_mean_exp_k = NaN(max(km_ordered),size(state_5_tpm,2));
    for k = 1:max(km_ordered)
        z_mean_exp_k(k,:) = mean(zscore(state_5_tpm(logical(km_genes(:,k)),:),0,2));
    end
     
    figure;
    for k = 1:max(km_ordered)
        subplot(1,4,k)
        scatter(X_5,Y_5,20,z_mean_exp_k(k,:),'filled');
        colormap(cmap_color_blind); 
        caxis([-0.5 0.5]);
        set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    end
        
end

if fig_3d
    load GBS_state_5.mat state_5_tpm info_5 X_5 Y_5
    load GBS_modules.mat km_ordered km_genes
    load ../prepare_the_data/general_info.mat module_colors
    
    %mean expression per module
    z_state_5 = zscore(state_5_tpm,0,2);
    z_mean_exp_k = NaN(size(km_genes,2),size(z_state_5,2));
    for k = 1:max(km_ordered)
        z_mean_exp_k(k,:) = mean(z_state_5(logical(km_genes(:,k)),:));
    end
    %normalize modules
    min_mean = z_mean_exp_k-min(z_mean_exp_k);
    norm_mean = min_mean./sum(min_mean);
    norm_mean(isnan(norm_mean)) = 0;
    [~, ind_coord] = sort(X_5+Y_5);
    norm_mean_sorted = norm_mean(:,ind_coord);
    
    figure; 
    h= bar(norm_mean_sorted',1,'stacked','edgecolor','none');
    h(1).FaceColor = module_colors(1,:);
    h(2).FaceColor = module_colors(2,:);
    h(3).FaceColor = module_colors(3,:);
    h(4).FaceColor = module_colors(4,:);
    ylim([0 1]);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    save km_single_cells.mat norm_mean_sorted
end

if fig_3e
    load km_single_cells.mat

    corr_color = flip(brewermap(64,'RdBu'));
    figure;
    imagesc(corrcoef(norm_mean_sorted(:,:)))
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    caxis([-1 1]);
    colormap(corr_color);
    
end

if fig_3f
    load km_single_cells.mat
    
    %running window on pdist betwwen individual cells
    win_length = 10;
    step = 5;
    clear P; count = 0;
    for k = 1:step:size(norm_mean_sorted,2)
        left_window_start = k - (win_length+1);
        left_window_end = k -1;
        right_window_start = k;
        right_window_end = k + win_length;
        
        if(min(left_window_start)<1 || right_window_end > size(norm_mean_sorted,2))
            continue
        end
        count = count+1;
        P(count) = pdist([mean(norm_mean_sorted(:,(left_window_start:left_window_end)),2)';mean(norm_mean_sorted(:,(right_window_start:right_window_end)),2)']);
    end
     
    num_box = 5;
    P_box = NaN(num_box,length(P));
    P_4 = repelem(1:num_box,round(length(P)/num_box)+1); P_4 = P_4(1:length(P));
    for n = 1:num_box
        P_box(n,1:sum(P_4 == n)) = P(P_4 == n);
    end
    
    figure;
    boxplot(P_box','color','k','Symbol','o','OutlierSize',2)
    
    ranksum(P_box(1,:),P_box(3,:))
    ranksum(P_box(5,:),P_box(3,:))
    
end

%sup figures:
fig_s3a    = 0; %PCA IL1B with bins
fig_s3b    = 0; %zavit with kmeans bar
fig_s3c    = 0; %silhouette
fig_s3d    = 0; %kmeans NMF comp
fig_s3e    = 0; %gene gene corr

if fig_s3a
    load ../prepare_the_data/general_info.mat gene_names
    load GBS_state_5.mat state_5_tpm X_5 Y_5 intervals
    load cmap_color_blind.mat
    
    figure;
    gene = ismember(gene_names,'IL1B');
    scatter(X_5,Y_5,20,zscore(state_5_tpm(gene,:),0,2),'filled'); hold on;
    colormap(cmap_color_blind); caxis([-1 1]);
    for i = 1:length(intervals)
        line([intervals(i) intervals(i)],[-60 60],'Color','k','LineStyle','--')
    end
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end

if fig_s3b
    load GBS_state_5.mat binned_mat info_5
    load GBS_modules.mat km_ordered
    load cmap_color_blind.mat
    load ../prepare_the_data/general_info.mat module_colors
    %zavit on dyn genes
    dyn_mat     = binned_mat(info_5,:);
    z_order     = zavit(dyn_mat,0);
    
    figure;
    subplot(1,10,1:9);
    imagesc(zscore(dyn_mat(z_order,:),0,2));
    colormap(cmap_color_blind); caxis([-2 2]);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(1,10,10);
    imagesc(km_ordered(z_order));
    colormap(gca,module_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    
end

if fig_s3c
    load GBS_modules.mat km_ordered dyn_mat
    load ../prepare_the_data/general_info.mat module_colors
    
    figure;
    [s,h]= silhouette(zscore(dyn_mat,0,2),km_ordered);
    b = h.CurrentAxes.Children(1);
    invalues = false; % Are we currently in "values" or in nans.
    group = 0;
    c = ones(length(b.XData), 1);
    for x = 1:length(b.XData)
        if invalues
            if isnan(b.YData(x))
                invalues = false;
            end
        elseif ~isnan(b.YData(x))
            % We are out of the nan area, back to the values area.
            group = group + 1;
            invalues = true;
        end
        c(x) = group;
    end
    b.FaceColor = 'flat';
    c(c==0) = 1;
    b.CData = module_colors(c,:);
    
    
end

if fig_s3d
    load GBS_state_5.mat z_binned info_5
    load GBS_modules.mat dyn_mat km_genes
    
    num_clust = 4;
    rng(7)
    [nmf_bin_genes,~] = nnmf(zscore(dyn_mat,0,2),num_clust);
    [~,ind_max] = max(nmf_bin_genes,[],2);
    %reorder the clusters
    NMF_ordered = zeros(size(ind_max));
    NMF_ordered(ind_max == 1) = 1;
    NMF_ordered(ind_max == 2) = 2;
    NMF_ordered(ind_max == 3) = 4;
    NMF_ordered(ind_max == 4) = 3;
    %genes vec per cluster
    NMF_genes = zeros(length(info_5),max(NMF_ordered));
    NMF_info  = zeros(size(info_5));
    NMF_info(info_5) = NMF_ordered;
    for k = 1:max(NMF_ordered)
        NMF_genes(:,k) = NMF_info == k;
    end
    %compare nmf and km clusters
    clust_names = {'cluster 1','cluster 2','cluster 3','cluster 4'};
    pvals = NaN(length(clust_names),length(clust_names));
    for n = 1:length(clust_names)
        pvals(:,n) = Enrichment(logical(NMF_genes(:,n)),0.00005,km_genes,clust_names);
    end
    figure;
    imagesc(-log10(pvals))
    
end

if fig_s3e
    load GBS_state_5.mat state_5_tpm info_5 z_binned
    load GBS_modules.mat km_ordered
    load ../prepare_the_data/general_info.mat module_colors
    
    z_info     = z_binned(info_5,:);
    c_info     = corrcoef(z_info');
    corr_color = flip(brewermap(64,'RdBu'));
    [~,sort_km] = sort(km_ordered);
    
    figure;
    subplot(6,1,1:5)
    imagesc(c_info(sort_km,sort_km));
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    caxis([-1 1]);
    colormap(corr_color);
    subplot(6,1,6)
    imagesc(km_ordered(sort_km)');
    colormap(gca,module_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end
