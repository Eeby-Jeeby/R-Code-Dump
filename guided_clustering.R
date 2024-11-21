#set up seurat data#
library(tidyverse)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)

LoadData("ifnb")
pbmc.data <- Read10X(data.dir = ifnb)
