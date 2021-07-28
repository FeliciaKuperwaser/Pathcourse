data_for_monocole = 0;
fig_2a            = 0; % states in GBS
fig_2b            = 0; % GBS unexp stack bar
fig_2c            = 0; % heatmap genes states + GO terms
fig_2def          = 0; % FACS data
fig_2g            = 0; % GBS_mut_states bar


if data_for_monocole
    load ../prepare_the_data/filtered_GBS_unexp.mat GBS_unexp_mat
    load ../prepare_the_data/GBS_macs.mat macs_vec macs_tpm_GBS_unexp
    load ../prepare_the_data/general_info.mat ribo_mito gene_names
    
    GBS_unexp_macs_raw = GBS_unexp_mat(:,macs_vec);
    dlmwrite('GBS_unexp_macs_mat.txt',GBS_unexp_macs_raw,'delimiter','\t')
    
    [info_genes] = informative_genes(macs_tpm_GBS_unexp,2,2, ribo_mito);  
    %manually save gene_names(info_genes)
    
end

if fig_2a
    load ../prepare_the_data/general_info.mat state_colors
    
    state_import      = readtable('state_vec_R_GBS_unexp');
    state_vec         = state_import{:,1};
    coord_import_x    = readtable('coord_x_GBS_unexp');
    coord_x           = coord_import_x{:,1};
    coord_import_y    = readtable('coord_y_GBS_unexp');
    coord_y           = coord_import_y{:,1};
    
    % state 2 has few cells, so i'll manually marge them with state 7 
    % to be consistente with fig 3, we will change the states to be:
    % 1 - non-activated macs (1)
    % 2 - mito (7)
    % 3 - antigene presentation (4)
    % 4 - apoptosis (5)
    % 5 - inflammatory (6)
    % 6 - stress (3)
    state_vec_new = zeros(size(state_vec));
    state_vec_new(state_vec == 1) = 1;
    state_vec_new(state_vec == 2) = 1;
    state_vec_new(state_vec == 3) = 1;
    state_vec_new(state_vec == 4) = 3;
    state_vec_new(state_vec == 5) = 4;
    state_vec_new(state_vec == 6) = 5;
    state_vec_new(state_vec == 7) = 2;
    
    figure;
    gscatter(coord_x,coord_y,state_vec_new,state_colors,'.',15);    
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    save GBS_unexp_states.mat state_vec_new coord_x coord_y

    
end

if fig_2b
    load ../prepare_the_data/GBS_macs.mat  macs_GBS_unexp_vec
    load GBS_unexp_states.mat state_vec_new
    load ../prepare_the_data/general_info.mat state_colors
    
    GBS_unexp_state = NaN(max(macs_GBS_unexp_vec),max(state_vec_new));
    for ab = 1:max(macs_GBS_unexp_vec)
        for s = 1: max(state_vec_new)
            GBS_unexp_state(ab,s) = sum((macs_GBS_unexp_vec' == ab)&(state_vec_new == s));
        end
    end
    perc_state_mut = GBS_unexp_state ./ sum(GBS_unexp_state,2);
   
    bac_names = {'GBS','unexp'}; 
    state_names = {'1','2','3','4','5'};
    
    figure;
    h= bar(perc_state_mut,'stacked');
    h(1).FaceColor = state_colors(1,:);
    h(2).FaceColor = state_colors(2,:);
    h(3).FaceColor = state_colors(3,:);
    h(4).FaceColor = state_colors(4,:);
    h(5).FaceColor = state_colors(5,:);
    legend(state_names);
    set(gca,'xtick',1:length(bac_names)); set(gca,'xticklabel',bac_names);
    ylim([0 1]);  
   
end

if fig_2c
   load GBS_unexp_states.mat
   load ../prepare_the_data/GBS_macs.mat macs_tpm_GBS_unexp macs_ft_GBS_unexp 
   load ../prepare_the_data/general_info.mat ribo_mito gene_names state_colors
   load GO_human_mat_10X.mat GO_mat_human_50_10X GO_names_50_10X
   load cmap_color_blind.mat
   
   names_no_ribo_mito = gene_names(~ribo_mito);
   no_ribo_mito_tpm = macs_tpm_GBS_unexp(~ribo_mito,:);
   GO_no_ribo_mito    = GO_mat_human_50_10X(~ribo_mito,:);
    
   %[p_values, ~, ~, ~, ~, ind_p_values] = diff_exp_genes(names_no_ribo_mito,state_vec_new, no_ribo_mito_tpm);
   %load states_diff_genes.mat p_values ind_p_values
    
    Pval_thersh = 0.005;
    P_genes = p_values < Pval_thersh; sum(P_genes)  
    % take top X genes
    gene_count = 100;
    genes_up   = ind_p_values(1:gene_count,:);    
    genes_up_P = genes_up.*(P_genes(genes_up));    
    diff_genes = genes_up_P(:); diff_genes(diff_genes == 0) = [];
    diff_mat   = no_ribo_mito_tpm(diff_genes,:);
    [~,ind_cells] = sort(state_vec_new);
    
    figure;
    subplot(10,1,1:9);
    imagesc(zscore(diff_mat(:,ind_cells),0,2))
    caxis([-1.5 1.5]);
    cells = 0;
    genes = 0;
    for s = 1:max(state_vec_new)
        rectangle('position',[cells genes sum(state_vec_new ==s) sum(genes_up_P(:,s)>0)],'LineWidth',1);
        cells = cells + sum(state_vec_new == s);
        genes = genes +sum(genes_up_P(:,s)>0);
    end
    colormap(gca,cmap_color_blind);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    subplot(10,1,10);
    imagesc(state_vec_new(ind_cells)');
    colormap(gca,state_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;

    % GO terms per state
    % change state and p thresh
    state = 3;
    p_thresh_go = 0.0005; 
    genes_to_take = genes_up_P(:,state);
    genes_to_take(genes_to_take ==0) = [];
    [~,GO_up] = Enrichment(genes_to_take, p_thresh_go, GO_no_ribo_mito, GO_names_50_10X);
    
    save states_diff_genes.mat p_values ind_p_values genes_up_P diff_genes diff_mat names_no_ribo_mito no_ribo_mito_tpm
    
end

if fig_2def
    unexp = csvread('export_unexp A H_Data Source - 1.csv', 1, 0);
    GBS   = csvread('export_GBS 4 A H_Data Source - 1.csv', 1, 0);
    unstained = csvread('export_unstained change gates_Data Source - 1.csv', 1, 0);
   
    % set values lower than 0 to 0
    unexp(unexp<0) = 0;
    GBS(GBS<0) = 0;
    unstained(unstained<0) = 0;
    
    % Channels
    FSC     = 1;
    PE      = 7; %GBS
    state_1 = 5; %CXCR4
    state_3 = 10;%CD80 
    state_5 = 8; %IL7R
    % set threshold for infection (based on unstained)
    cells_thresh   = 2.7;
    GBS_thresh     = 3.7;
    state_1_thresh = 3.2;
    state_5_thresh = 4;
    state_3_thresh = 3.5;
    PE_thresh      = 3.5;

    bys_cells   = log10(1+GBS(:,PE))> cells_thresh & log10(1+GBS(:,PE))<= GBS_thresh;
    inf_cells   = log10(1+GBS(:,PE))> GBS_thresh;
    GBS_cells   = log10(1+GBS(:,PE)) > cells_thresh;
    unexp_cells = log10(1+unexp(:,PE))> cells_thresh;
    unsta_cells = log10(1+unstained(:,PE))> cells_thresh;
     
    % per state- are you infected or bystendar
    states_sum = NaN(3,2); %three states and infected/bystanders
    %state 1,3,5
    states_sum(1,1) = sum(log10(1+GBS(inf_cells,state_1))>state_1_thresh);
    states_sum(1,2) = sum(log10(1+GBS(bys_cells,state_1))>state_1_thresh);  
    states_sum(2,1) = sum(log10(1+GBS(inf_cells,state_3))>state_3_thresh);
    states_sum(2,2) = sum(log10(1+GBS(bys_cells,state_3))>state_3_thresh);
    states_sum(3,1) = sum(log10(1+GBS(inf_cells,state_5))>state_5_thresh);
    states_sum(3,2) = sum(log10(1+GBS(bys_cells,state_5))>state_5_thresh);
    states_freq = states_sum./sum(states_sum,2);
     
    % state 1 vs 5
    Q_samples = NaN(4,3); % 4 Q's and 3 conditions
    %infected
    Q_samples(1,1) = sum(log10(1+GBS(inf_cells,state_1)) > state_1_thresh & log10(1+GBS(inf_cells,state_5)) <= state_5_thresh); %high state 1 low state 5
    Q_samples(2,1) = sum(log10(1+GBS(inf_cells,state_1)) > state_1_thresh & log10(1+GBS(inf_cells,state_5)) > state_5_thresh); %high state 1 high state 5
    Q_samples(3,1) = sum(log10(1+GBS(inf_cells,state_1)) <= state_1_thresh & log10(1+GBS(inf_cells,state_5)) <= state_5_thresh); %low state 1 low state 5
    Q_samples(4,1) = sum(log10(1+GBS(inf_cells,state_1)) <= state_1_thresh & log10(1+GBS(inf_cells,state_5)) > state_5_thresh); %low state 1 high state 5
    %bys
    Q_samples(1,2) = sum(log10(1+GBS(bys_cells,state_1)) > state_1_thresh & log10(1+GBS(bys_cells,state_5)) <= state_5_thresh); %high state 1 low state 5
    Q_samples(2,2) = sum(log10(1+GBS(bys_cells,state_1)) > state_1_thresh & log10(1+GBS(bys_cells,state_5)) > state_5_thresh); %high state 1 high state 5
    Q_samples(3,2) = sum(log10(1+GBS(bys_cells,state_1)) <= state_1_thresh & log10(1+GBS(bys_cells,state_5)) <= state_5_thresh); %low state 1 low state 5
    Q_samples(4,2) = sum(log10(1+GBS(bys_cells,state_1)) <= state_1_thresh & log10(1+GBS(bys_cells,state_5)) > state_5_thresh); %low state 1 high state 5 
    %unexp
    Q_samples(1,3) = sum(log10(1+unexp(unexp_cells,state_1)) > state_1_thresh & log10(1+unexp(unexp_cells,state_5)) <= state_5_thresh); %high state 1 low state 5
    Q_samples(2,3) = sum(log10(1+unexp(unexp_cells,state_1)) > state_1_thresh & log10(1+unexp(unexp_cells,state_5)) > state_5_thresh); %high state 1 high state 5
    Q_samples(3,3) = sum(log10(1+unexp(unexp_cells,state_1)) <= state_1_thresh & log10(1+unexp(unexp_cells,state_5)) <= state_5_thresh); %low state 1 low state 5
    Q_samples(4,3) = sum(log10(1+unexp(unexp_cells,state_1)) <= state_1_thresh & log10(1+unexp(unexp_cells,state_5)) > state_5_thresh); %low state 1 high state 5
    Q_percent = Q_samples./sum(Q_samples);
    
    state_names = {'state 1','state 3','state 5'}; 
    inf_names   = {'infected','bystanders'};
    inf_colors  = [0.75 0.75 0.75; 0 0 0];
    cond_names = {'infected','bystanders','unexposed'};
    colors = brewermap(12,'Set3');
    Q_colors = colors(1:4,:);
    
    figure;
    h= bar(states_freq,'stacked');
    h(1).FaceColor = inf_colors(1,:);
    h(2).FaceColor = inf_colors(2,:);
    legend(inf_names);
    set(gca,'xtick',1:length(state_names)); set(gca,'xticklabel',state_names);
    ylim([0 1]);  
    figure;
    plot(log10(1+GBS(bys_cells,state_5)),log10(1+GBS(bys_cells,state_1)),'.','color','k'); hold on;
    plot(log10(1+GBS(inf_cells,state_5)),log10(1+GBS(inf_cells,state_1)),'.','color',[0.75 0.75 0.75]); 
    xlim([2.5 6]); ylim([2 5]);
    line([state_5_thresh state_5_thresh],[2 5],'color','k');
    line([2.5 6],[state_1_thresh state_1_thresh],'color','k');
    title('GBS'); xlabel('state 5'); ylabel('state 1');
    figure;
    h= bar(Q_percent','stacked');
    h(1).FaceColor = Q_colors(1,:);
    h(2).FaceColor = Q_colors(2,:);
    h(3).FaceColor = Q_colors(3,:);
    h(4).FaceColor = Q_colors(4,:);
    set(gca,'xtick',1:length(cond_names)); set(gca,'xticklabel',cond_names);
    xtickangle(45);
    legend({'Q1','Q2','Q3','Q4'});
    
end

if fig_2g
    load GBS_unexp_states.mat
    load ../prepare_the_data/macs_mut.mat macs_ft_mut macs_ab_mut
    load ../prepare_the_data/GBS_macs.mat macs_tpm_GBS_unexp macs_ft_GBS_unexp macs_GBS_unexp_vec
    load ../prepare_the_data/general_info.mat ribo_mito gene_names state_colors
  
    [info_macs] = informative_genes(macs_tpm_GBS_unexp,2,2, ribo_mito);
    m = fitcknn(macs_ft_GBS_unexp(info_macs,:)',state_vec_new);
    c = predict(m,macs_ft_mut(info_macs,:)');
     
    mut_state = NaN(max(macs_ab_mut),max(c));
    for ab = 1:max(macs_ab_mut)
        for s = 1: max(c)
            mut_state(ab,s) = sum((macs_ab_mut == ab)&(c == s));
        end
    end
    perc_state_mut = mut_state ./ sum(mut_state,2);
  
    bac_names = {'WT','cpsE','cylE','unexp'}; 
    state_names = {'1','2','3','4','5'};
    
    figure;
    h= bar(perc_state_mut,'stacked');
    h(1).FaceColor = state_colors(1,:);
    h(2).FaceColor = state_colors(2,:);
    h(3).FaceColor = state_colors(3,:);
    h(4).FaceColor = state_colors(4,:);
    h(5).FaceColor = state_colors(5,:);
    legend(state_names);
    set(gca,'xtick',1:length(bac_names)); set(gca,'xticklabel',bac_names);
    ylim([0 1]);  
    
end

%sup figures:
fig_s2ab       = 0; % poodle GBS/unexp and pseudotime 
fig_s2cd       = 0; % tsne states and GBS/unexp
fig_s2ef      = 0; % FACS state 1/5 unstained and unexp
fig_s2ghi     = 0; % inffered poodle mut states and conditions and stack bar mut replicates

if fig_S2ab
    load GBS_unexp_states.mat coord_x coord_y
    load ../prepare_the_data/GBS_macs.mat macs_GBS_unexp_vec
    load ../prepare_the_data/general_info.mat bac_colors
    
    pseudotime_import = readtable('pseudotime_vec_R_GBS_unexp');
    pseudotime_vec    = pseudotime_import{:,1};
    Jet_cmap = jet;
    [~,Xi] = sort(coord_x);
    rand_vec = randperm(RandStream('mt19937ar','Seed',7),length(macs_GBS_unexp_vec));
    figure;
    scatter(coord_x(rand_vec),coord_y(rand_vec),20,pseudotime_vec(rand_vec),'filled')
    colormap(gca,Jet_cmap);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    figure;
    scatter(coord_x(rand_vec),coord_y(rand_vec),20,macs_GBS_unexp_vec(rand_vec),'filled')
    colormap(gca,bac_colors(7:8,:));
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
end

if fig_S2cd
   load ../prepare_the_data/GBS_macs.mat 
   load ../prepare_the_data/general_info.mat ribo_mito gene_names bac_colors state_colors
   load cmap_color_blind.mat
   load ../figure_2/GBS_unexp_states.mat
   
   rand_vec = randperm(RandStream('mt19937ar','Seed',7),length(macs_GBS_unexp_vec));

   info_macs = informative_genes(macs_tpm_GBS_unexp,2,2,ribo_mito);
   rng(7);
   T = tsne(macs_ft_GBS_unexp(info_macs,:)','Perplexity', 15);
   T1 = T(:,1);
   T2 = T(:,2);
   figure;
   scatter(T1(rand_vec),T2(rand_vec),20,state_vec_new(rand_vec),'filled');
   colormap(state_colors);
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
   
   figure;
   scatter(T1(rand_vec),T2(rand_vec),20,macs_GBS_unexp_vec(rand_vec),'filled');
   colormap(bac_colors(7:8,:));
   set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
 
   
end

if fig_S2ef
    unstained = csvread('export_unstained change gates_Data Source - 1.csv', 1, 0);
    unexp = csvread('export_unexp A H_Data Source - 1.csv', 1, 0);

    % set values lower than 0 to 0
    unstained(unstained<0) = 0;
    unexp(unexp<0) = 0;

    % set threshold for infection
    PE  = 7;
    FSC = 1;
    cells_thresh = 2.7;
    unsta_cells = log10(1+unstained(:,PE))> cells_thresh;    
    unexp_cells = log10(1+unexp(:,PE))> cells_thresh;    
    state_1 = 5;
    state_5 = 8;
    
    %set thresh by unstained
    figure;
    plot(log10(1+unstained(unsta_cells,state_5)),log10(1+unstained(unsta_cells,state_1)),'.','color','k');
    line([4 4],[1 4],'color','k');
    line([2.5 5],[3.2 3.2],'color','k');
    title('unstained'); xlabel('state 5'); ylabel('state 1');

    figure;
    plot(log10(1+unexp(unexp_cells,state_5)),log10(1+unexp(unexp_cells,state_1)),'.','color','k');
    line([4 4],[0 6],'color','k');
    line([2.5 5.5],[3.2 3.2],'color','k');
    title('unexp'); xlabel('state 5'); ylabel('state 1');

end

if fig_S2ghi
    load GBS_unexp_states.mat
    load ../prepare_the_data/macs_mut.mat macs_ft_mut macs_ab_mut macs_donor_mut
    load ../prepare_the_data/GBS_macs.mat macs_tpm_GBS_unexp macs_ft_GBS_unexp macs_GBS_unexp_vec
    load ../prepare_the_data/general_info.mat ribo_mito gene_names bac_colors state_colors
    
    [info_macs] = informative_genes(macs_tpm_GBS_unexp,2,2, ribo_mito);
    m = fitcknn(macs_ft_GBS_unexp(info_macs,:)',state_vec_new);
    c = predict(m,macs_ft_mut(info_macs,:)');
     
    mut_state = NaN(max(macs_ab_mut),max(c),max(macs_donor_mut));
    for ab = 1:max(macs_ab_mut)
        for s = 1: max(c)
            for d = 1:max(macs_donor_mut)
                mut_state(ab,s,d) = sum((macs_ab_mut == ab)&(c == s)&(macs_donor_mut' == d));
            end
        end
    end
    perc_state_mut = mut_state ./ sum(mut_state,2);
  
    bac_names = {'WT','cpsE','cylE','unexp'}; 
    state_names = {'1','2','3','4','5'};
    
    figure;
    for d = 1:max(macs_donor_mut)
        subplot(1,2,d)
        h= bar(perc_state_mut(:,:,d),'stacked');
        h(1).FaceColor = state_colors(1,:);
        h(2).FaceColor = state_colors(2,:);
        h(3).FaceColor = state_colors(3,:);
        h(4).FaceColor = state_colors(4,:);
        h(5).FaceColor = state_colors(5,:);
        legend(state_names);
        set(gca,'xtick',1:length(bac_names)); set(gca,'xticklabel',bac_names);
        ylim([0 1]);
        title(['donor  ' num2str(d)])
    end
    
    num_neighbors = 5;
    m_x = fitcknn(macs_ft_GBS_unexp(info_macs,:)',coord_x,'NumNeighbors',num_neighbors,'BreakTies','nearest');
    c_x = predict(m_x,macs_ft_mut(info_macs,:)');
    m_y = fitcknn(macs_ft_GBS_unexp(info_macs,:)',coord_y,'NumNeighbors',num_neighbors,'BreakTies','nearest');
    c_y = predict(m_y,macs_ft_mut(info_macs,:)');
    
    figure;
    scatter(coord_x,coord_y,20,[0.75 0.75 0.75],'filled'); hold on;
    scatter(c_x,c_y,20,c,'filled');
    colormap(gca,state_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
    
    pinks = brewermap(11,'PiYG');
    mut_colors = [bac_colors(7,:);pinks(1:2,:);bac_colors(8,:)];
    figure;
    scatter(coord_x,coord_y,20,[0.75 0.75 0.75],'filled'); hold on;
    scatter(c_x,c_y,20,macs_ab_mut,'filled');
    colormap(gca,mut_colors);
    set(gca,'xtick',[]); set(gca,'ytick',[]); box off;
 
    
end

