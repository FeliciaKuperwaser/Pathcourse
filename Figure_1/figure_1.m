%fig_1a is schematic 
%fig_1b is FACS histograms 
fig_1c         = 0; % PCA two clusters
fig_1d         = 0; % genes PC1
fig_1ef        = 0; % macs TNF and histogram GBS/unexp on PC2
fig_1g         = 0; % FACS TNF
fig_1h         = 0; % heatmap binned expression

if fig_1c
    load ../prepare_the_data/general_info.mat 
    load ../prepare_the_data/filtered_GBS_unexp.mat
    
    % PCA on informative genes
    [info_genes] = informative_genes(GBS_unexp_tpm,3,2, ribo_mito);
    sum(info_genes)
    [coeff,score_G,~,~,explained,~] = pca(GBS_unexp_ft(info_genes,:)');
    X_G = score_G(:,1);
    Y_G = score_G(:,2);
     
    % color by GBS or unexp
    GBS_unexp_color = bac_colors(7:8,:);
    rand_vec = randperm(RandStream('mt19937ar','Seed',7),length(GBS_unexp_vec));
    figure;
    scatter(X_G(rand_vec),Y_G(rand_vec),20,GBS_unexp_vec(rand_vec),'filled');
    colormap(GBS_unexp_color);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    rectangle('position',[-10 -38 85 155],'LineWidth',1,'LineStyle','--');
    
    save fig_1C.mat coeff score_G info_genes X_G Y_G GBS_unexp_color
   
end

if fig_1d
    load ../prepare_the_data/general_info.mat 
    load ../prepare_the_data/filtered_GBS_unexp.mat GBS_unexp_tpm reaction_vec GBS_unexp_vec
    load fig_1C.mat
    load cmap_color_blind.mat
   
    % genes that contribute the most to PC1
    info_genes_names = gene_names(info_genes);
    thresh_coeff = sqrt(1/(sum(info_genes)));
    PC = 1;
    genes_coeff = abs(coeff(:,PC)) > thresh_coeff; 
    genes_coeff_names = info_genes_names(genes_coeff);
    genes_coeff_vec = zeros(size(info_genes));
    genes_coeff_vec(info_genes) = genes_coeff;
    mat_inform = GBS_unexp_tpm(info_genes,:);
    z_coeff    = zscore(mat_inform(genes_coeff,:),0,2);
    [z_order]  = zavit(mat_inform(genes_coeff,:),-60);
    [~,Xi]     = sort(X_G);     
    
    %donor vec
    donor_vec = zeros(size(reaction_vec));
    donor_vec(reaction_vec == 1 | reaction_vec == 2) = 1;
    donor_vec(reaction_vec == 3 | reaction_vec == 4) = 2;
    donor_vec(reaction_vec == 5 | reaction_vec == 6) = 3;
    donor_vec(reaction_vec == 7 | reaction_vec == 8) = 4;
    
    %mono/macs ratio
    mono_gene_list = {'B2M','HLA-A','HLA-B','HLA-C'};
    macs_gene_list = {'HLA-DRA','HLA-DRB1','CD4','APOE','CD68'};
    mono_vec = ismember(gene_names, mono_gene_list);
    macs_vec = ismember(gene_names, macs_gene_list);
    mono_sum = sum(GBS_unexp_tpm(mono_vec,:));
    macs_sum = sum(GBS_unexp_tpm(macs_vec,:));
    mono_macs_log = log2(mono_sum./macs_sum);
    
    figure;
    subplot(14,1,1:10);
    imagesc(z_coeff(z_order,Xi)); %order cells by PCA coordinates
    set(gca,'ytick',1:length(genes_coeff_names)); set(gca,'yticklabel',genes_coeff_names(z_order),'FontSize',8);
    set(gca,'xtick',[]);
    colormap(cmap_color_blind);
    caxis([-2 2]);
    subplot(14,1,11);
    imagesc(X_G(Xi)');
    colormap(gca,'Jet'); 
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    subplot(14,1,12);
    imagesc(GBS_unexp_vec(Xi));
    colormap(gca,GBS_unexp_color);
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    subplot(14,1,13);
    imagesc(donor_vec(Xi));
    colormap(gca,[rep_colors(1:2,:);rep_colors(5:6,:)]);
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    subplot(14,1,14);
    imagesc(mono_macs_log(Xi));
    colormap(gca,'Parula');
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    
end

if fig_1ef
    load ../prepare_the_data/general_info.mat gene_names ribo_mito 
    load ../prepare_the_data/filtered_GBS_unexp.mat GBS_unexp_vec GBS_unexp_tpm
    load ../prepare_the_data/GBS_macs.mat macs_vec macs_GBS_unexp_vec
    load fig_1C.mat GBS_unexp_color score_G
    load cmap_color_blind.mat
    
    PC2 = score_G(:,2);
    PC3 = score_G(:,3);
 
    figure;
    subplot(5,1,1:4)
    TNF = ismember(gene_names,'TNF');
    scatter(PC2(macs_vec),PC3(macs_vec),20,zscore(GBS_unexp_tpm(TNF,macs_vec),0,2),'filled');
    colormap(gca,cmap_color_blind);
    caxis([-1 1]);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(5,1,5);
    PC2_macs = PC2(macs_vec);
    histogram(PC2_macs(macs_GBS_unexp_vec == 1),100,'EdgeColor', GBS_unexp_color(1,:),'FaceColor',GBS_unexp_color(1,:)); hold on; % GBS first
    histogram(PC2_macs(macs_GBS_unexp_vec == 2),100,'EdgeColor', GBS_unexp_color(2,:),'FaceColor',GBS_unexp_color(2,:)); % unexp second
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end

if fig_1g
    GBS_4A  = csvread('TNF_FACS/export_GBS 4A  _Data Source - 1_Single Cells.csv', 1);
  % set negative values to zero
    log_GBS_A   = log10(1+ (GBS_4A > 0) .* GBS_4A);

    PE_Cy7  = 10; PE = 7;
    inf_thresh = 3; % the threshold based on unexp (sup)
   
    figure;
    plot(log_GBS_A(:,PE),log_GBS_A(:,PE_Cy7),'.','color','k');
    xlim([0 6]);ylim([0 6]);
    line([3 3],[0 6],'color','k');
    line([0 6],[3 3],'color','k');
    xlabel('GBS infection'); ylabel('TNF exp'); 
    title('GBS donor A');
   
end

if fig_1h
    load ../prepare_the_data/GBS_macs.mat
    load ../prepare_the_data/general_info.mat gene_names
    load cmap_color_blind.mat
    load GO_human_mat_10X.mat GO_mat_human_50_10X GO_names_50_10X
    %load fig_1D.mat
    load fig_1C.mat GBS_unexp_color score_G
    
    PC2 = score_G(:,2);
    PC2_macs = PC2(macs_vec);
    %bin genes in GBS only
    number_of_intervals = 10;
    binned_mat = NaN(length(gene_names),number_of_intervals);
    GBS_mat    = macs_tpm_GBS_unexp(:,macs_GBS_unexp_vec == 1); % take only GBS cells
    GBS_X_Gm   = PC2_macs(macs_GBS_unexp_vec == 1);
    min_cord = min(GBS_X_Gm);
    max_cord = max(GBS_X_Gm)+eps; % I added eps for the mat_interval
    intervals = linspace(min_cord, max_cord, number_of_intervals + 1); % 10 bins have 11 edges
    categories = discretize(GBS_X_Gm, intervals);
    for i = 1:number_of_intervals
        binned_mat(:,i) = mean(GBS_mat(:,categories==i),2);
    end
    z_binned = zscore(binned_mat,0,2);
    
    %dyn genes zavit
    fano_factor = var(binned_mat,[],2)./mean(binned_mat,2);
    dyn_genes   = fano_factor>4*nanmean(fano_factor); sum(dyn_genes)
    dyn_names   = gene_names(dyn_genes);
    dyn_mat     = binned_mat(dyn_genes,:);
    z_order     = zavit(dyn_mat,0);
    figure;
    subplot(1,5,1:4);
    imagesc(zscore(dyn_mat(z_order,:),0,2));
    colormap(cmap_color_blind); caxis([-2 2]);
    set(gca,'ytick',[]);
    EML_vec = ones(size(dyn_names));
    EML_vec(140:165) = 2;
    EML_vec(166:end) = 3;
    subplot(1,5,5);
    imagesc(EML_vec); 
    colormap(gca,'Jet');
    set(gca,'ytick',1:size(dyn_names,1)); set(gca,'yticklabel',dyn_names(z_order));
    set(gca,'ytick',[]); set(gca,'xtick',[]);

    dyn_names_orderd = dyn_names(z_order);
 
    save EML.mat dyn_names_orderd EML_vec

    % GO terms per state
    % change EML and p thresh
    EML = 3;
    p_thresh_go = 0.000005; 
    names_to_take = dyn_names_orderd(EML_vec == EML);
    genes_to_take = ismember(gene_names,names_to_take);
    [~,GO_up] = Enrichment(genes_to_take, p_thresh_go, GO_mat_human_50_10X, GO_names_50_10X);
   
end

% sup figures;
fig_s1a        = 0; % FACS 50% inf 
fig_s1b        = 0; % dendrogram and PCA two clusters
fig_s1cd       = 0; % FACS TNF unexp and boxplot
fig_s1e        = 0; % PC2 binned exp
fig_s1f        = 0; % GO terms heatmap E/M/L

if fig_s1a
   load ../prepare_the_data/general_info.mat rep_colors
   precent_inf = xlsread('precent_inf.xlsx');

   d_names = {'donor A','donor B','donor C','donor D'};
   d_colors = [rep_colors(1:2,:);rep_colors(5:6,:)];

   cond_names = {'GBS','Unexposed'};
   figure;
   for d = 1:length(d_names)
       plot(1:length(cond_names),precent_inf(7:8,d),'.','color',d_colors(d,:),'MarkerSize',20); hold on; %take only GBS and unexposd
   end
   xlim([0.5 2.5]); ylim([0 100]);
   ylabel('precent infection')
   set(gca,'xtick',1:length(cond_names)); set(gca,'xticklabel',cond_names); xtickangle(45)
   legend(d_names); box off;
end

if fig_s1b
    load ../prepare_the_data/filtered_GBS_unexp.mat GBS_unexp_ft GBS_unexp_tpm
    load ../prepare_the_data/GBS_macs.mat cluster_id_GBS_unexp
    load ../prepare_the_data/general_info.mat ribo_mito gene_names
    
    two_colors =[0.9769 0.9839 0.0805; 0.2422 0.1504 0.6603];
    [info_genes] = informative_genes(GBS_unexp_tpm,3,2, ribo_mito);
    sum(info_genes)
    [coeff,score_G,~,~,explained,~] = pca(GBS_unexp_ft(info_genes,:)');
    X_G = score_G(:,1);
    Y_G = score_G(:,2);
    
    figure;
    subplot(10,1,1:4);
    [~,~,cells_order] = dendrogram(linkage(pdist(GBS_unexp_ft(info_genes,:)'),'ward'),0);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    subplot(10,1,5);
    imagesc(cluster_id_GBS_unexp(cells_order));
    colormap(gca,two_colors)
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    subplot(10,1,6:10);
    scatter(X_G,Y_G,20,cluster_id_GBS_unexp,'filled')
    colormap(two_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    
end

if fig_s1cd
    unexp_A = csvread('TNF_FACS/export_unexposed A_Data Source - 1_Single Cells.csv', 1);
    GBS_4A   = csvread('TNF_FACS/export_GBS 4A  _Data Source - 1_Single Cells.csv', 1);

    % set negative values to zero
    log_unexp_A = log10(1+ (unexp_A > 0) .* unexp_A);
    log_GBS_A   = log10(1+ (GBS_4A > 0) .* GBS_4A);
    PE_Cy7  = 10; PE = 7;
    inf_thresh = 3; % the threshold based on unexp
    
    figure;
    subplot(1,5,1:3);
    plot(log_unexp_A(:,PE),log_unexp_A(:,PE_Cy7),'.','color','k');
    xlim([0 6]);ylim([0 6]);
    line([3 3],[0 6],'color','k');
    line([0 6],[3 3],'color','k');
    xlabel('GBS infection'); ylabel('TNF exp');
    title('unexp donor A');
    subplot(1,5,4:5);
    TNF_box = NaN(size(log_GBS_A,1),3);
    TNF_box(1:sum(log_unexp_A(:,PE)<= inf_thresh),1) = log_unexp_A(log_unexp_A(:,PE)<= inf_thresh,PE_Cy7);
    TNF_box(1:sum(log_GBS_A(:,PE)<= inf_thresh),2) = log_GBS_A(log_GBS_A(:,PE)<= inf_thresh,PE_Cy7);
    TNF_box(1:sum(log_GBS_A(:,PE)> inf_thresh),3) = log_GBS_A(log_GBS_A(:,PE)> inf_thresh,PE_Cy7);
    boxplot(TNF_box,'color','k','outliersize',2,'Symbol','ko')
    ylabel('TNF exp')
    names = {'unexposed','non-infected','infected'};
    set(gca,'xtick',1:length(names)); set(gca,'xticklabel',names); xtickangle(45)
    [~,p_12] = ttest2(TNF_box(:,1),TNF_box(:,2));
    [~,p_13] = ttest2(TNF_box(:,1),TNF_box(:,3));
    [~,p_23] = ttest2(TNF_box(:,2),TNF_box(:,3));
    
end

if fig_s1e
   load fig_1C.mat GBS_unexp_color score_G
   load ../prepare_the_data/general_info.mat gene_names
   load ../prepare_the_data/GBS_macs.mat macs_vec macs_GBS_unexp_vec macs_tpm_GBS_unexp
   load cmap_color_blind.mat
   
   PC2      = score_G(:,2);
   PC3      = score_G(:,3);
   PC2_macs = PC2(macs_vec);
   PC3_macs = PC3(macs_vec);
   
   number_of_intervals = 10;
   GBS_X_Gm   = PC2_macs(macs_GBS_unexp_vec == 1);
   min_cord = min(GBS_X_Gm);
   max_cord = max(GBS_X_Gm)+eps; % I added eps for the mat_interval
   intervals = linspace(min_cord, max_cord, number_of_intervals + 1); % 10 bins have 11 edges
   
   figure;
   TNF = ismember(gene_names,'TNF');
   scatter(PC2_macs,PC3_macs,20,zscore(macs_tpm_GBS_unexp(TNF,:),0,2),'filled'); hold on;
   colormap(cmap_color_blind); caxis([-1 1]);
   for i = 1:length(intervals)
      line([intervals(i) intervals(i)],[-30 50],'Color','k','LineStyle','--') 
   end
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;

end

if fig_s1f
    load EML.mat
    load GO_human_mat_10X.mat GO_mat_human_50_10X GO_names_50_10X
    load ../prepare_the_data/general_info.mat gene_names
    load cmap_color_blind.mat
    
    Pvals = NaN(length(GO_names_50_10X),max(EML_vec));
    for e = 1:max(EML_vec) 
        names_EML = ismember(gene_names,dyn_names_orderd(EML_vec == e));
        Pvals(:,e) = Enrichment(names_EML, 0.0000000001,GO_mat_human_50_10X,GO_names_50_10X);
    end
    
    number_of_GOs = 5;
    p_ind    = NaN(number_of_GOs,max(EML_vec));
    for c = 1:max(EML_vec)
        [~,xi] = sort(Pvals(:,c));
        p_ind(:,c) = xi(1:number_of_GOs);
    end
    p_ind_join = p_ind(:);
    p_val_ind  = Pvals(p_ind_join,:);
    figure;
    imagesc(-log10(p_val_ind)); 
    colorbar; caxis([0 10]);
    set(gca,'ytick',1:size(p_ind_join,1)); set(gca,'yticklabel',GO_names_50_10X(p_ind_join));
       
end

