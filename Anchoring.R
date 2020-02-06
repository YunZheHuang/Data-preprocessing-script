library(dplyr)
library(Seurat)
library(Matrix)
library(reticulate)
library(devtools)
library(ggplot2)
py_config()
py_install("umap-learn")
options(future.globals.maxSize = 4000 * 1024^2)


sample5.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019936/filtered_feature_bc_matrix/")
sample5 <- CreateSeuratObject(counts = sample5.rna, project = "sample5")

sample6.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019938/filtered_feature_bc_matrix")
sample6 <- CreateSeuratObject(counts = sample6.rna, project = "sample6")

sample7.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019939/filtered_feature_bc_matrix")
sample7 <- CreateSeuratObject(counts = sample7.rna, project = "sample7")

Blood.B <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Bloodsample2/",gene.column = 2)
Blood2 <- CreateSeuratObject(counts = Blood.B, project = "Blood2")


BloodMerged <- merge(sample5, y = c(sample6, sample7, Blood2), add.cell.ids = c("S5", "S6", "S7", "B2"), project = "BloodMerged", merge.data = FALSE)


Blood.B <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Bloodsample2/",gene.column = 2)
Blood1 <- CreateSeuratObject(counts = Blood.B, project = "Blood1")


BEI.A<- Read10X(data.dir = "/Users/yunzhehuang/Downloads/BEI/",gene.column = 2)
BEI.RNA <- BEI.A[["Gene Expression"]]
Blood2 <- CreateSeuratObject(counts = BEI.RNA, project = "Blood2")


BloodMerged <- merge(sample5, y = sample6, add.cell.ids = c("S5", "S6"), project = "BloodMerged", 
                                        merge.data = FALSE)

BloodList <- SplitObject(BloodMerged)

for (i in 1:length(x = BloodList)) {
  BloodList[[i]] <- SCTransform(object = BloodList[[i]], verbose = FALSE)
    
}


Blood.features <- SelectIntegrationFeatures(object.list = BloodList, nfeatures = 4000)

BloodList <- PrepSCTIntegration(object.list = BloodList, anchor.features = Blood.features, 
                                    verbose = FALSE)





BloodAnchors <- FindIntegrationAnchors(object.list = BloodList, dims = 1:20, normalization.method = "SCT", anchor.features = Blood.features)

BloodIntegrated <- IntegrateData(anchorset = BloodAnchors, normalization.method = "SCT",dims = 1:20, verbose = FALSE)

DefaultAssay(object = BloodIntegrated) <- "integrated"


BloodIntegrated <- RunPCA(object = BloodIntegrated, npcs = 20, verbose = FALSE)
BloodIntegrated <- RunUMAP(object = BloodIntegrated, reduction = "pca", 
                               dims = 1:20)

DimPlot(object = BloodIntegrated, reduction = "umap")


BloodIntegrated <- FindNeighbors(BloodIntegrated, reduction = "pca", dims = 1:20)
BloodIntegrated <- FindClusters(BloodIntegrated, resolution = 0.5)
BloodIntegrated <- RunTSNE(BloodIntegrated, dims = 1:20, method = "FIt-SNE")

TSNEPlot(BloodIntegrated, methods = "FIt-SNE")


BloodMerged <- NormalizeData(BloodMerged, assay = "RNA")

BloodMerged <- FindVariableFeatures(BloodMerged, selection.method = "vst")
BloodMerged <- ScaleData(BloodMerged)
BloodMerged <- RunPCA(BloodMerged, features = VariableFeatures(object = BloodMerged))

BloodMerged <- FindNeighbors(BloodMerged, reduction = "pca", dims = 1:20)
BloodMerged <- FindClusters(BloodMerged, resolution = 0.5)
BloodMerged <- RunTSNE(BloodMerged, dims = 1:20, method = "FIt-SNE")

TSNEPlot(BloodMerged, methods = "FIt-SNE", label = TRUE)

WhichCells(BloodMerged, idents = '8')



--------------------------------
  
RPS.genes <- grep(pattern = "^RPS", x = rownames(x = sample6@assays[["RNA"]]@data), value = TRUE)

sample6 <- NormalizeData(sample6, assay = "RNA")

sample6 <- FindVariableFeatures(sample6, selection.method = "vst", nfeatures = 4000)

var.genes.no.RPS <- dplyr::setdiff(sample6@assays[["RNA"]]@var.features, RPS.genes)

IGHV.genes <- grep(pattern = "^IGHV", x = rownames(x = sample6@assays[["RNA"]]@data), value = TRUE)

sample6 <- RunPCA(object = sample6, pc.genes = var.genes.no.RPS, do.print = TRUE, pcs.print = 1:12, genes.print = 10)

VariableFeatures(sample6) <- c("FOXP3","CTLA4","CD19","S100A9","IGHM","IGHD","CD45","CD34","XBP1","IRF4","CD9","CD38","ITGB7","VCAM1",
                               "CD138","FCRL5","FCLR4","S100A6","CD29","CD30","FBXW12","ZBTB32","TBX21","NR4A1","NR4A3","PDCD1","IFI44",
                               "SPRY2","CD79A","MRC1","CD163","CCR4","CCR6","CCR7","IL2RA","FCGR3A","CD14","ITGAX")
sample6@assays[["RNA"]]@var.features <- c(sample6@assays[["RNA"]]@var.features, IGHV.genes)
            
sample6@assays[["RNA"]]@var.features > "FOXP3"
sample6 <- ScaleData(sample6)
sample6 <- FindVariableFeatures(sample6, selection.method = "vst")
sample6 <- RunPCA(sample6, features = VariableFeatures(object = sample6))

sample6 <- FindNeighbors(sample6, reduction = "pca", dims = 1:20)
sample6 <- FindClusters(sample6, resolution = 0.5)
sample6 <- RunTSNE(sample6, dims = 1:20, method = "FIt-SNE")

TSNEPlot(sample6, methods = "FIt-SNE", label = TRUE)


sample6 <- RunUMAP(object = sample6, reduction = "pca", 
                           dims = 1:20)

DimPlot(object = sample6, reduction = "umap", label = TRUE)


FeaturePlot(sample6, features = c("CCR4"), ncol = 2)
















