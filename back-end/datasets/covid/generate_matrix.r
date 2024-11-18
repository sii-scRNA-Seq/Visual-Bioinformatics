library(Seurat)

print("Loading Object")
seurat_obj <- readRDS('Covid_555_1.rds')

print(seurat_obj)
print(colnames(seurat_obj@meta.data))

raw_counts <- GetAssayData(seurat_obj, layer='counts', assay='RNA')

dir.create("./temp")

print("Saving metadata")
write.table(seurat_obj@meta.data[,c("Donor", "Status")], file="temp/covid_meta.csv", sep=",") # keeps the rownames
print("Saved metadata to temp/covid_meta.csv")

print("Saving Matrix")
write.table(t(raw_counts), file="temp/covid_matrix.csv", sep=",") # keeps the rownames
print("Saved matrix to temp/covid_matrix.csv")

