%In vivo analysis

identify_cell_types        = 0; %Step 1 %also generated part A for Sup Figure
control_cell_analysis      = 1; %for SUP  /// FIRST RUN: human_mouse_comparison
classify_the_infected_macs = 0; %Part B of figure
gene_heatmap_for_figure    = 0; %Part C
Profile_Human_modules      = 0; %Part D
mid_infection_analysis     = 0; %Part E


load('/Users/yanaii01/Dropbox/Itai/Projects/Pep38_GBS_IY_access/IY_analysis/cmap_color_blind.mat')

try
    in_vivo_tpm;
catch
    load in_vivo_data.mat
    load Lincs.mat
    load gene_names_human.mat gene_names_human
    load GO_mouse_mat_10X.mat
    A=GO_mat_mouse_10X;
    load GO_human_mouse_compatible GO_human_R GO_mouse_R GO_names mouse_gene_names human_gene_names;
end

if identify_cell_types
    %make more groups: cell_annot
    %define the top 200 genes for each cluster: gene2cluster
    %look at the GO categories: cluster2GO
    %compare with human GO categories: 1. get modules, 2. make GO, 3. compare

    %identify_red_blood_cells
    [info_genes] = informative_genes(in_vivo_tpm,2,2, ribo_mito+Linc_annot');
    sum(info_genes)
    [coeff,score_G,~,~,explained,~] = pca(in_vivo_ft(info_genes,:)');
    X_G = score_G(:,1); Y_G = score_G(:,2);
    
    %color by gene
    gene = ismember(gene_names,'Hbb-bt');
    figure;     subplot(1,3,1);
    
    scatter(X_G,Y_G,30,zscore(in_vivo_tpm(gene,:),0,2),'filled');
    caxis([0 2]); %colorbar;
    title(gene_names(gene)); xlabel('PC1'); ylabel('PC2');
    colormap(cmap_color_blind);
  
    info_genes_names = gene_names(info_genes);
    thresh_coeff = sqrt(1/(sum(info_genes)));
    PC = 1;
    genes_coeff = abs(coeff(:,PC)) > thresh_coeff;
    genes_coeff_names = info_genes_names(genes_coeff);
    genes_coeff_vec = zeros(size(info_genes));
    genes_coeff_vec(info_genes) = genes_coeff;
    mat_inform = in_vivo_tpm(info_genes,:);
    z_coeff    = zscore(mat_inform(genes_coeff,:),0,2);
    [z_order]  = zavit(mat_inform(genes_coeff,:),-60);
    [~,Xi]     = sort(X_G);
    
    %     subplot(1,2,2); imagesc(z_coeff(z_order,Xi));
    %     set(gca,'ytick',1:length(genes_coeff_names)); set(gca,'yticklabel',genes_coeff_names(z_order));
    %     xlabel('PC1-ordered');
    
    %filter out high PC1
    cell_annot(X_G>50) = 1;
    filtered_cells = X_G<50;
    filtered_tpm = in_vivo_tpm(:,filtered_cells);
    
    save cell_annot.mat cell_annot filtered_tpm;

    % identify_lymphocytes
    cell_ind = setdiff(1:length(cell_annot),find(cell_annot==1));
    filtered_ft = in_vivo_ft(:,cell_ind);
    [info_genes_f] = informative_genes(filtered_tpm,2,2, ribo_mito+Linc_annot');
    sum(info_genes_f)
    [coeff_f,score_f,~,~,explained,~] = pca(filtered_ft(info_genes_f,:)');
    X_f = score_f(:,1);
    Y_f = score_f(:,2);
    
    %     figure;    plot(X_f,Y_f,'.','MarkerSize',20)
    %     set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    group_A = find(X_f<-43);
    group_B = find(X_f>=-43);
    
    cell_annot(cell_ind(find(X_f<-43))) = 2;
    save cell_annot.mat cell_annot;
    
    %     gene = ismember(gene_names,'Nfkb1');
    %     figure; scatter(X_f,Y_f,30,zscore(in_vivo_tpm(gene,cell_ind),0,2),'filled');
    %     caxis([-1 4]); colorbar; title(gene_names(gene)); xlabel('PC1'); ylabel('PC2');
    
    gene = ismember(gene_names,'Cd2');
    subplot(1,3,2);
    scatter(X_f,Y_f,30,zscore(in_vivo_tpm(gene,cell_ind),0,2),'filled');
    caxis([-1 4]); %colorbar;
    title(gene_names(gene)); xlabel('PC1'); ylabel('PC2');
    colormap(cmap_color_blind);
  
    info_genes = find(informative_genes(in_vivo_tpm,1,1, ribo_mito+Linc_annot'));
    P = ones(length(gene_names),1);
    
    c = cell_ind(group_A);
    non_c = cell_ind(group_B);
    for j = 1:length(info_genes)
        [P(info_genes(j),1),~] = ranksum(in_vivo_ft(info_genes(j),c),in_vivo_ft(info_genes(j),non_c), 'tail', 'right');
    end
    [x,y] = sort(P);
    
    [i,j]=Enrichment(y(1:50),0.00005,GO_mouse_R,GO_names)
    [x,y] = sort(j);
    GO_names(y(1:10))
    
    %identify_apoptotic_cells
    
    bins = 2;
    
    load cell_annot.mat cell_annot;
    cell_ind     = setdiff(1:length(cell_annot),union(find(cell_annot==1),find(cell_annot==2)));
    filtered_ft  = in_vivo_ft(:,cell_ind);
    filtered_tpm = in_vivo_tpm(:,cell_ind);
    
    [info_genes_f] = informative_genes(filtered_ft,4,4, ribo_mito+Linc_annot');
    sum(info_genes_f)
    [coeff_f,score_f,~,~,explained,~] = pca(filtered_ft(info_genes_f,:)');
    X_f = score_f(:,1); Y_f = score_f(:,2);
    
    cell_annot(cell_ind(find(X_f<0))) = 3; %these are the presumed appoptotic
    %     cell_annot(cell_ind(find(Y_f<-80))) = 3;
    %     cell_annot(cell_ind(intersect(find(X_f>=0),find(Y_f>=-80))))  = 4; %these are the presumed infected macs
    
    gene = ismember(gene_names,'Nfkb1');
    gene = ismember(gene_names,'Lars2');
    %     figure; scatter(X_f,Y_f,30,zscore(in_vivo_tpm(gene,cell_ind),0,2),'filled');
    %     caxis([-1 4]); colorbar; title(gene_names(gene)); xlabel('PC1'); ylabel('PC2');
    
    %     figure; scatter(X_f,Y_f,30,cell_annot(cell_ind),'filled');
    
    mito = strmatch('mt-',gene_names);
    mito_exp = mean(in_vivo_tpm(mito,:));
    
    subplot(1,3,3); scatter(X_f,Y_f,30,mito_exp(cell_ind),'filled');
    title('Mitochondrial expression'); xlabel('PC1'); ylabel('PC2');
    
    print -painters -depsc 'Fig_S6_A.pdf'
end

if (control_cell_analysis)
    top_genes_in_cluster = 200;
    Y_or_X_for_split     = 0; %0 is X a-xis
    
    try
        control_ft;
    catch
        load control_data.mat control_ft control_tpm
    end
    
    %first filter out the lympocytes
    [info_genes_control] = informative_genes(control_tpm,2,2, ribo_mito+Linc_annot');
    sum(info_genes_control)
    [coeff,score_C,~,~,explained,~] = pca(control_ft(info_genes,:)');
    X_C = score_C(:,1); Y_C = score_C(:,2); 
    
    %color by gene
    gene = ismember(gene_names,'Cd74');
    figure; subplot(1,3,1); scatter(X_C,Y_C,30,zscore(control_tpm(gene,:),0,2),'filled');
    caxis([-1 3]); title(gene_names(gene)); xlabel('PC1'); ylabel('PC2');    
    % c = colorbar; set(c,'ytick',[]);
    
    cell_annot_control = zeros(1,size(control_tpm,2));
    %filter out high PC1
    cell_annot_control(X_C<-10) = 1;
    
    cells_select = find(cell_annot_control==0);
    [info_genes_control] = informative_genes(control_tpm(:,cells_select),2,2, ribo_mito+Linc_annot');
    sum(info_genes_control)
    [coeff,score_C,~,~,explained,~] = pca(control_ft(info_genes,cells_select)');
    X_C = score_C(:,1); Y_C = score_C(:,2); Z_C = score_C(:,3);
    
    %color by gene
    gene = find(ismember(gene_names,'Il6'));    
    subplot(1,3,2); 
    scatter(X_C+Y_C,Z_C,30,zscore(control_tpm(gene,cells_select),0,2),'filled');
    caxis([-1 1]); title(gene_names(gene)); xlabel('PC1+PC2'); ylabel('PC2');
        
    [i, xi] = sort(X_C+Y_C);
    [i, Z_rank_C] = sort(xi);
    [i, Z_rank_C_i] = sort(Z_rank_C);
     Z_rank_C_i = Z_rank_C_i(end:-1:1);
    
    %heatmap
     clear Gene_exp;
     genes_to_show = {'Il1b','Ccl4','Gadd45a','Nfkbiz','Ier3','Zfp36','Ppp1r15a','Ctsc','Ier2','Tnfaip3','Nfkbia','Fn1','Tmsb4x','Cybb','C1qa','C1qb','Cyba','Aplp2','Ftl1','Fabp4','Psap','Fth1','Ctss',...
        'Il6','Cd14','Il1a', 'Fcer1g'};
    genes_to_show_module_number = [4 4 4 4 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 1 ...
            4 1 4 2];

    for s = 1:length(genes_to_show)
        gene_i = strmatch(genes_to_show(s),gene_names,'exact');        
        Gene_score = control_ft(gene_i,:);
        Gene_exp(s,:) = movmean(Gene_score(cells_select(Z_rank_C_i)),10);
    end
     A = zscore(Gene_exp,0,2);
    [i,z] = sort(genes_to_show_module_number); %z = z(end:-1:1);

    subplot(1,3,3); %subplot('position',[0.1 0.4 0.8 0.5])
    imagesc(A(z,:),[-3,3]);
    set(gca,'ytick',1:length(genes_to_show));
    set(gca,'yticklabel',genes_to_show(z));
    colormap(cmap_color_blind);
    print -painters -depsc 'Fig_S6_B.pdf'

end

if (classify_the_infected_macs)
    bins = 2;
    
    %recompute the PCA
    cell_ind = setdiff(1:length(cell_annot),union(find(cell_annot==3),union(find(cell_annot==1),find(cell_annot==2))));
    filtered_ft = in_vivo_ft(:,cell_ind);
    filtered_tpm = in_vivo_tpm(:,cell_ind);
    
    [info_genes_f] = informative_genes(filtered_tpm,2,2,ribo_mito+Linc_annot');
    sum(info_genes_f)
    [coeff_f,score_f,~,~,explained,~] = pca(filtered_tpm(info_genes_f,:)');
    X_f = score_f(:,1); Y_f = score_f(:,2); Z_f = score_f(:,3);
    
    [i, xi] = sort(X_f+Y_f);
    [i, Z_rank] = sort(xi);
    [i, Z_rank_i] = sort(Z_rank);

    figure; scatter(X_f+Y_f,Z_f,30,Z_rank,'filled'); 
    colormap(gca,cmocean('haline'))
    set(gca,'xtick',[]);    set(gca,'ytick',[]);   
    cb = colorbar; set(cb,'ytick',[]);
    xlabel('PC1+PC2','fontsize',20);ylabel('PC3','fontsize',20);
    print -painters -depsc 'Fig6_B1.pdf'

    [i, xi] = sort(Z_rank);
    sign_posts = i(1:floor(length(Z_rank)/bins):length(Z_rank));
    sign_posts(length(sign_posts)+1) = 100000; %limit to last sign post
    
    for i = 1:bins
        c = intersect(find(Z_rank>=sign_posts(i)),find(Z_rank<=sign_posts(i+1)));
        cell_annot(cell_ind(c)) = i+3;
    end
    
    gene_i = 'Il6';
    gene = strmatch(gene_i,gene_names,'exact');
    exp  = zscore(in_vivo_tpm(gene,cell_ind),0,2);
    figure; %subplot(1,3,1);
    scatter(X_f+Y_f,Z_f,30,exp,'filled');caxis([-1,.5]);
    xlabel('PC1+PC2','fontsize',20);ylabel('PC3','fontsize',20);
    colormap(cmap_color_blind);
    set(gca,'xtick',[]);    set(gca,'ytick',[]);    colorbar;
    cb = colorbar; set(cb,'ytick',[]);
    title(gene_i,'fontsize',20);
    print -painters -depsc 'Fig6_B2.pdf'

end

if (Profile_Human_modules)
    load module_mouse_genes.mat %these are the orthologs
    clear Terms_exp;
    for sets = 1:4
        genes_s = find(modules_mouse(:,sets));
        GS_score = sum(zscore(in_vivo_ft(genes_s,:),0,2));
        GS_score = median((in_vivo_ft(genes_s,:)));
        Terms_exp(sets,:) = movmean(GS_score(cell_ind(Z_rank_i)),100);
    end
    A = zscore(Terms_exp,0,2);
    figure; imagesc(A,[-2,2]);
    set(gca,'ytick',1:4);    set(gca,'yticklabel',1:4);
    ylabel('Modules','fontsize',16);
    set(gca,'xtick',[]); c=colorbar; set(c,'ytick',[]);
    xlabel('Cell trajectory','fontsize',16);
    colormap(cmap_color_blind);     
    print -painters -depsc 'Fig6_D.pdf'
    
end

if (gene_heatmap_for_figure)
    clear Gene_exp;
    module_colors = [0.6196 0.0039 0.2588 ; 0.9569 0.4275 0.2627 ; 0.4000 0.7608 0.6471 ; 0.3686 0.3098 0.6353];
    
    genes_to_show = {'Il1b','Ccl4','Gadd45a','Nfkbiz','Ier3','Zfp36','Ppp1r15a','Ctsc','Ier2','Tnfaip3','Nfkbia','Fn1','Tmsb4x','Cybb','C1qa','C1qb','Cyba','Aplp2','Ftl1','Fabp4','Psap','Fth1','Ctss',...
        'Il6','Cd14','Il1a', 'Fcer1g'};
    
    genes_to_show_module_number = [4 4 4 4 3 3 3 3 3 3 3 2 2 2 1 1 1 1 1 1 1 1 1 ...
        4 1 4 2];
    
    Gene_exp = [];
    for s = 1:length(genes_to_show)
        gene_i = strmatch(genes_to_show(s),gene_names,'exact');        
        Gene_score = in_vivo_ft(gene_i,:);
        Gene_exp(s,:) = movmean(Gene_score(cell_ind(Z_rank_i)),20);
    end
    A = zscore(Gene_exp,0,2);
    z = zavit(A,0); %z = z(end:-1:1);
    
    [i,z] = sort(genes_to_show_module_number); %z = z(end:-1:1);

    height = 0.4; space = 0.025; left = 0.4; width = 0.5;
    top = 0.975 - height;
    %%subplot('position',[left bottom width height])
    figure;  subplot('position',[left top width height])
    imagesc(A(z,:),[-3,3]);
    set(gca,'xtick',[]);
    set(gca,'ytick',1:length(genes_to_show));
    set(gca,'yticklabel',genes_to_show(z),'fontsize',8);
    colormap(gca,cmap_color_blind); h = colorbar; set(h,'ytick',[]);
    m = min(find(cell_annot(cell_ind(Z_rank_i))==5));
    
    subplot('position',[0.9 top 0.05 height]) 
    imagesc(genes_to_show_module_number(z)'); 
    set(gca,'xtick',[]);     set(gca,'ytick',[]); 
    colormap(gca,module_colors); 

    %What are the GO terms for state 5 modulews
    %what are the genes
    terms = {'lysosome','extracellular exosome','endosome',...
        'cellular response to tumor necrosis factor','positive regulation of ERK1 and ERK2 cascade',...
        'inflammatory response' ,...
        'cellular response to interleukin-1','cytokine activity'};
    for s = 1:length(terms)
        GO_s = strmatch(terms(s),GO_names,'exact');
        genes_s = find(GO_mouse_R(:,GO_s));
        GO_score = sum(zscore(in_vivo_ft(genes_s,:),0,2));
        GO_score = sum((in_vivo_ft(genes_s,:)));
        Terms_exp(s,:) = movmean(GO_score(cell_ind(Z_rank_i)),50);
    end
    
    top = top - height - space;
    subplot('position',      [left top width height]);
    A = zscore(Terms_exp,0,2);
    z = zavit(A,0); z = z(end:-1:1);
    imagesc(A(z,:),[-2,2]);
    set(gca,'ytick',1:length(terms));    set(gca,'xtick',[]);    set(gca,'yticklabel',terms(z))       
    colormap(gca,cmap_color_blind); h = colorbar; set(h,'ytick',[]);

    height2 = 0.075;
    top = top - height2 - space;        
    subplot('position',        [left top width height2])
    imagesc(sort(Z_rank_i)'); colorbar;
    set(gca,'xtick',[]);    set(gca,'ytick',[]);
    xlabel('Cell trajectory'); set(gca,'xtick',[]); 
    colormap(gca,cmocean('haline')); h = colorbar; set(h,'ytick',[]);
    
    print -painters -depsc 'Fig6_C.pdf'
end

if (mid_infection_analysis)
    [info_genes] = informative_genes(in_vivo_tpm(:,cell_ind),3,3,ribo_mito+Linc_annot');
    sum(info_genes)

    B = nan(2000,5);
    limits = round([1:size(cell_ind,2)/5:size(cell_ind,2), size(cell_ind,2)]);
    for i = 1:length(limits)-1
        cells = limits(i):(limits(i+1));
        C = corrcoef(in_vivo_tpm(info_genes,cell_ind(Z_rank_i(cells))));
        C_u = triu(C,1);
        b = C_u(find(C_u));
        B(1:length(b),i) = 1-b;
    end
    figure;         boxplot(B,'color','k','Symbol','o','OutlierSize',2)
    set(gca,'xtick',[])
    xlabel('binned single cells');
    ylabel('pairwise distance');
        
    ranksum(B(1:(min(find(isnan(B(:,1))))-1),1),B(1:min(find(isnan(B(:,3))))-1,3))
    ranksum(B(1:(min(find(isnan(B(:,5))))-1),5),B(1:min(find(isnan(B(:,3))))-1,3))

    print -painters -depsc 'Fig6_E.pdf'
end