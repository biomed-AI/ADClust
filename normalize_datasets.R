library(SeuratDisk)
library(Seurat)


norm_data<-function(name){
    # base path is the raw datasets with seurat object format
	base="./seurat_object/"
    # save the data to the dir data of ADClust
	save_base="./data/"
	
	print(paste0('loading ',name))
	data=get(load(paste0(base, name, ".RData")))
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

    count=data@assays$RNA@data[data@assays$RNA@var.features,]
    meta.data=data@meta.data
    batches <- CreateSeuratObject(count, meta.data = meta.data, project = "Seurat3_benchmark")
    filename=paste0(save_base,name,'_seurat_norm.h5Seurat')
    SaveH5Seurat(batches, filename=filename)
    print(paste0('converting ',name))
    Convert( filename, dest = "h5ad")


}

filename=c(
   "Baron_Human"
)

for(i in filename){
  norm_data(i)
}
