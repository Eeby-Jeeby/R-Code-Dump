library(Seurat)
library(SeuratData)

InstallData("panc8")
data("panc8")

pancreas_list <- SplitObject(panc8, split.by = "tech")
pancreas_list <- pancreas_list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

#perform standard preprocessing (log-normalization) 
#and identify variable features individually for each
for (i in 1:length(pancreas_list)) {
  pancreas_list[[i]] <- NormalizeData(pancreas_list[[i]], verbose = FALSE)
  pancreas_list[[i]] <- FindVariableFeatures(pancreas_list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}

##Integration of 3 pancreatic islet cell datasets##
#identify anchors using FindIntegrationAnchors() function, which
#takes a list of Seurat objects as input
reference_list <- pancreas_list[c("celseq", "celseq2", "smartseq2")]
pancreas_anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:30)

#pass these anchors to the IntegrateData() function, which returns a Seurat object
pancreas_integrated <- IntegrateData(anchorset = pancreas_anchors, dims = 1:30)

#After running IntegrateData(), the Seurat object will contain a 
#new Assay with the integrated expression matrix

#We can then use this new integrated matrix for downstream analysis
#and visualization. Here we scale the integrated data, run PCA, and
#visualize the results with UMAP
library(tidyverse)
library(lubridate)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas_integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas_integrated <- ScaleData(pancreas_integrated, verbose = FALSE)
pancreas_integrated <- RunPCA(pancreas_integrated, npcs = 30, verbose = FALSE)
pancreas_integrated <- RunUMAP(pancreas_integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas_integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2

##Cell type classification using an integrated reference##
#After finding anchors, we use TransferData() function to classify
#the query cells based on reference data
#TransferData() returns a matrix with predicted IDs and prediction
#scores, which we can add to the query metadata
pancreas_query <- pancreas_list[["fluidigmc1"]]
pancreas_anchors <- FindTransferAnchors(reference = pancreas_integrated, query = pancreas_query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas_anchors, refdata = pancreas_integrated$celltype,
                            dims = 1:30)
pancreas_query <- AddMetaData(pancreas_query, metadata = predictions)

#we find that there is a high agreement in cell type classification,
#with over 96% of cells being labeled correctly
pancreas_query$prediction_match <- pancreas_query$predicted.id == pancreas_query$celltype
table(pancreas_query$prediction_match)

#To verify this further, we can examine some canonical cell type markers for specific pancreatic islet cell populations
table(pancreas_query$predicted.id)
VlnPlot(pancreas_query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

##Unimodal UMAP Projection##
#projection of a query onto the reference UMAP structure
#computing the reference UMAP model and then call MapQuery() instead of TransferData()

pancreas_integrated <- RunUMAP(pancreas_integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas_query <- MapQuery(anchorset = pancreas_anchors, reference = pancreas_integrated, query = pancreas_query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

#MapQuery() is a wrapper around 3 functions: TransferData(), IntegrateEmbeddings(), and ProjectUMAP()
#TransferData() is used to transfer cell type labels and impute the ADT values
#IntegrateEmbeddings() is used to integrate reference with query by correcting the 
#query’s projected low-dimensional embeddings
#ProjectUMAP() is used to project the query data onto the UMAP structure of the 
#reference

#Visualize the query cells along our reference
p1 <- DimPlot(pancreas_integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2