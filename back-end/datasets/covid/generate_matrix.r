library(Seurat)

objects <- c("Covid_555_1", "HC_HIP002")
for (object in objects) {
  print(object)
  print("Loading Object")
  seurat_obj <- readRDS(paste0(object,".rds"))
  raw_counts <- GetAssayData(seurat_obj, layer='counts', assay='RNA')

  dir.create("./temp")

  print("Saving metadata")
  path <- paste0('temp/', object, '_meta.csv')
  write.table(seurat_obj@meta.data[,c("Donor", "Status")], file=path, sep=",") # keeps the rownames
  print(paste0("Saved metadata to", path))

  print("Saving Matrix")
  path <- paste0('temp/', object, '_matrix.csv')
  write.table(t(raw_counts), file=path, sep=",") # keeps the rownames
  print(paste0("Saved matrix to", path))
}

