# Set the correct default repository
# So you don't have to choose one every time you install a package
r = getOption("repos")
r["CRAN"] = "http://cran.rstudio.com"
options(repos = r)

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("knitr")) {
  install.packages("knitr")
  library(knitr)
}

if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require("tibble")) {
  install.packages("tibble")
  library(tibble)
}

if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)
}

if (!require("leiden")) {
  install.packages("leiden")
  library(leiden)
}

if (!require("igraph")) {
  install.packages("igraph")
  library(igraph)
}

if (!require("MAST")) {
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("MAST")
  library(MAST)
}

if (!require("Seurat")) {
  install.packages("Seurat")
  library(Seurat)
}

if (!require("scatterplot3d")) {
  install.packages("scatterplot3d")
  library(scatterplot3d)
}

if (!require("gplots")) {
  install.packages("gplots")
  library(gplots)
}

if (!require("fields")) {
  install.packages("fields")
  library(fields)
}

if (!require("hablar")) {
  install.packages("hablar")
  library(hablar)
}

if (!require('gprofiler2')) {
  install.packages("gprofiler2")
  library(gprofiler2)
}

if (!require("devtools")) {
  install.packages("devtools")
  library(devtools)
}

if (!require("ggdendro")) {
  install.packages("ggdendro")
  library(ggdendro)
}

if (!require("mosaic")) {
  install.packages("mosaic")
  library(mosaic)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

if (!require("patchwork")) {
  install.packages("patchwork")
  library(patchwork)
}

if (!require('BiocNeighbors')) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("BiocNeighbors")
  library(BiocNeighbors)
}

if (!require("ggridges")) {
  install.packages("ggridges")
  library(ggridges)
}

if (!require('viridis')) {
  install.packages("viridis")
  library(viridis)
}

if (!require('rrvgo')) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rrvgo")
  library(rrvgo)
}

if (!require('org.Hs.eg.db')) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
}

if (!require('Nebulosa')) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Nebulosa")
  library(Nebulosa)
}

if (!require('Nebulosa')) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Nebulosa")
  library(Nebulosa)
}

if (!require("NMF")) {
  install.packages('NMF')
  library(NMF)
}

if (!require("circlize")) {
  install.packages("circlize")
  library(circlize)
}

if (!require('ComplexHeatmap')) {
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap)
}

if (!require(CellChat)) {
  if (!require(devtools)) {
    install.packages("devtools")}
  install_github("sqjin/CellChat")
  library(CellChat)
}

if (!require(NeuronChat)) {
  if (!require(devtools)) {
    install.packages("devtools")}
  install_github("Wei-BioMath/NeuronChat")
  library(NeuronChat)
}
