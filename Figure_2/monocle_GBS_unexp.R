#### Installing and loading packages ####

setwd("C:/Users/Gal/Dropbox (NYU Langone Health)/Projects/Pap36_Pathcourse/analysis/analysis_by_figure/new_2020_5_16/Figure_2")

options(warn=-1)

#source("http://bioconductor.org/biocLite.R")
#biocLite("monocle")
library(monocle)

#### Loading data ####
## you need raw matrix, gene names as .txt, sample/time point/ treatment IDs as .txt.
expr = as.matrix(read.delim("C:/Users/Gal/Dropbox (NYU Langone Health)/Projects/Pap36_Pathcourse/analysis/analysis_by_figure/new_2020_5_16/Figure_2/GBS_unexp_macs_mat.txt", header=FALSE))
gene_names <- read.table("C:/Users/Gal/Dropbox (NYU Langone Health)/Projects/Pap36_Pathcourse/analysis/analysis_by_figure/new_2020_5_16/Figure_2/gene_names", quote="\"", comment.char="", stringsAsFactors = FALSE, header = FALSE)[,1]
rownames(expr) = gene_names

#### creating it myself
pd = new("AnnotatedDataFrame", data = data.frame('Cell' = colnames(expr), row.names = colnames(expr)))
fd = new("AnnotatedDataFrame", data = data.frame('gene_short_name' = rownames(expr), row.names = rownames(expr)))
cds = newCellDataSet(expr, phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


#cds <- load_lung()
cds
View(pData(cds))

info_gene_names <- read.table("C:/Users/Gal/Dropbox (NYU Langone Health)/Projects/Pap36_Pathcourse/analysis/analysis_by_figure/new_2020_5_16/Figure_2/info_genes", quote="\"", comment.char="", stringsAsFactors = FALSE, header = FALSE)[,1]

#temp_genes = sample(ordering_genes, 1000)

cds <- setOrderingFilter(cds, info_gene_names)
#plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
#cds <- orderCells(cds, root_state = "5")

#plot_cell_trajectory(cds, color_by = "orig.ident")
plot_cell_trajectory(cds, color_by = "State")


#### Building tree ####

#cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
#cds <- orderCells(cds)

#plot_cell_trajectory(cds, color_by = "my_sample")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")

pseudotime_vec = cds$Pseudotime
write(pseudotime_vec,file = "pseudotime_vec_R_GBS_unexp",sep = "\n")

state_vec = cds$State
write(state_vec,file = "state_vec_R_GBS_unexp",sep = "\n")

all_coord = reducedDimS(cds)
write(all_coord[1,],file = "coord_x_GBS_unexp",sep = "\n")
write(all_coord[2,],file = "coord_y_GBS_unexp",sep = "\n")

