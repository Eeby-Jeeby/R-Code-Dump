#how to mitigate the effects of cell cycle heterogeneity in scRNA-seq data by 
#calculating cell cycle phase scores based on canonical markers, and regressing these 
#out of the data during pre-processing

library(Seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "C:/Users/Ryan/Desktop/research_r/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
#using FindVariableFeatures, we see that while most of the variance can be explained
#by lineage, PC8 and PC10 are split on cell-cycle genes including TOP2A and MKI67
marrow <- ScaleData(marrow, features = rownames(marrow))

#attempt to regress this signal from the data, so that cell-cycle heterogeneity does 
#not contribute to PCA or downstream analysis
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
DimHeatmap(marrow, dims = c(8, 10))

##Assign Cell-Cycle Scores##
#1. we assign each cell a score, based on its expression of G2/M and S phase markers
#assign scores in the CellCycleScoring() function, which stores S and G2/M scores in object meta data
#CellCycleScoring() also set the identity of the Seurat object to the cell-cycle 
#phase by passing set.ident = TRUE
marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

##Regress out cell cycle scores during data scaling##
#attempt to subtract (‘regress out’) this source of heterogeneity from the data
marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

##ALTERNATE WORKFLOW##
#procedure above removes all signal associated with cell cycle
#in some cases, it can negatively impact downstream analysis, particularly in 
#differentiating processes

#As an alternative, we suggest regressing out the difference between the G2M and S 
#phase scores
#means that signals separating non-cycling cells and cycling cells will be maintained,
#but differences in cell cycle phase among proliferating cell will be regressed out 
#of the data

marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))

# cell cycle effects strongly mitigated in PCA
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)

# when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# cells however, within actively proliferating cells, G2M and S phase cells group together
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)
