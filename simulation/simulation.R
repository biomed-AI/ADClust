rm(list = ls())
library(splatter)
library(Seurat)
library(SeuratDisk)

dropout.rate <- c()
facScale_list=c(0.4, 0.35,0.3,0.25,0.2)
batchCells_input=2000
nGroups_input=4
for (facScale in facScale_list){
for(i in 1:5) {
  simulate <- function(nGroups=nGroups_input, nGenes=2000, batchCells=batchCells_input)
  {
    if (nGroups > 1) method <- 'groups'
    else             method <- 'single'

   group.prob <-c(0.1,0.15, 0.25,0.5)

    sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                         dropout.type="experiment", method=method, seed=1000+i,
    de.facScale=facScale, dropout.mid=2)

    print(str(sim))

    counts     <- as.data.frame(t(counts(sim)))
    truecounts <- as.data.frame(t(assays(sim)$TrueCounts))

    dropout    <- assays(sim)$Dropout
    mode(dropout) <- 'integer'

    cellinfo   <- as.data.frame(colData(sim))
    geneinfo   <- as.data.frame(rowData(sim))

    list(sim=sim,
         counts=counts,
         cellinfo=cellinfo,
         geneinfo=geneinfo,
         truecounts=truecounts)
  }

  sim <- simulate()

  simulation <- sim$sim
  counts <- sim$counts
  geneinfo <- sim$geneinfo
  cellinfo <- sim$cellinfo
  truecounts <- sim$truecounts

  dropout.rate <- c(dropout.rate, (sum(counts==0)-sum(truecounts==0))/sum(truecounts>0))
    print(paste0(" dropout.rate", dropout.rate))
    print(paste0("zero ratio: ", sum(counts==0)/(sum(counts==0)+sum(counts>0))))

  X <- t(counts)
  Y <- as.integer(substring(cellinfo$Group,6))
  Y <- Y-1

     reference_sample = colnames(X)
      reference_frame=data.frame(as.character(Y), row.names=reference_sample)
      reference_frame[[1]]=as.character(reference_frame[[1]])
      names(reference_frame) <-NULL


    names(reference_frame) <-"CellType"
    data <- CreateSeuratObject(X, meta.data = reference_frame, project = "Seurat3_benchmark")


    # # save seurat object data
    root_path="./"
    save(data,file=paste0(root_path, "splatter_facScale_",facScale,"_",i , ".RData"))



    # # save hd5a count data
    namestr=paste0("splatter_facScale_",facScale,"_",i)
    setwd("./")
    filename=paste0(root_path,namestr,'_count_data.h5Seurat')
    if (file.exists(filename)){        file.remove(filename)}
    SaveH5Seurat(data, filename=filename)

    if (file.exists(paste0(root_path,namestr,'_count_data.h5ad'))){
        file.remove(paste0(root_path,namestr,'_count_data.h5ad'))}
    Convert( filename, dest = "h5ad")



    # save top2000 norm data
    save_base="./"
    name=paste0("splatter_facScale_",facScale,"_",i)
    data <- NormalizeData(data, normalization.method='RC')
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

    count=data@assays$RNA@data

    meta.data=data@meta.data
    batches <- CreateSeuratObject(count, meta.data = meta.data, project = "Seurat3_benchmark")

    filename=paste0(save_base,name,'_seurat_norm.h5Seurat')
    if (file.exists(filename)){        file.remove(filename)}
    SaveH5Seurat(batches, filename=filename)
    print(paste0('converting ',name))

    if (file.exists(paste0(root_path,namestr,'_seurat_norm.h5ad'))){
    file.remove(paste0(root_path,namestr,'_seurat_norm.h5ad'))}
    Convert( filename, dest = "h5ad")

    print(namestr)

}
}

