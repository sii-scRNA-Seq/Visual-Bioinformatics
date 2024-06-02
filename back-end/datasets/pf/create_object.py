import scanpy as sc
import anndata
from scipy import sparse
import pandas 

print("reading day0")
day0 = sc.read_10x_mtx("./day0_filtered_feature_bc_matrix")
print("reading day3")
day3 = sc.read_10x_mtx("./day3_filtered_feature_bc_matrix")
print("reading day5")
day5 = sc.read_10x_mtx("./day5_filtered_feature_bc_matrix")
print("reading day10a")
day10a = sc.read_10x_mtx("./day10a_filtered_feature_bc_matrix")
print("reading day10b")
day10b = sc.read_10x_mtx("./day10b_filtered_feature_bc_matrix")

day0.obs["day"] = "0"
day3.obs["day"] = "3"
day5.obs["day"] = "5"
day10a.obs["day"] = "10a"
day10b.obs["day"] = "10b"

print("concatenating")
adata = anndata.concat([day0, day3, day5, day10a, day10b])
adata.obs_names_make_unique()

print("convert day metadata to category")
cat_type = pandas.CategoricalDtype(categories=['0', '3', '5', '10a', '10b'], ordered=True)
adata.obs["day"] = adata.obs["day"].astype(cat_type)

print("sparsify matrix for save")
adata.X = sparse.csr_matrix(adata.X) 

sc.pp.subsample(adata, n_obs=10000, random_state=7, copy=False)

print("saving")
adata.write_h5ad('MCA_PF_DOGGA.h5ad')
