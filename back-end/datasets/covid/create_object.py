from pathlib import Path
import scanpy as sc
import pandas as pd
from scipy import sparse

print("reading data")
meta = pd.read_csv("./temp/Covid_555_1_meta.csv")
adata = sc.read_csv(Path("./temp/Covid_555_1_matrix.csv"))
adata.obs = meta

meta_hc = pd.read_csv("./temp/HC_HIP002_meta.csv")
adata_hc = sc.read_csv(Path("./temp/HC_HIP002_matrix.csv"))
adata_hc.obs = meta_hc

adata = adata.concatenate(adata_hc)

# Convert strings to categories to reduce file size
for col in ['Donor', 'Status']:
    adata.obs[col] = adata.obs[col].astype('category')

print(adata)

print("Compressing matrix")
sparse_X = sparse.csr_matrix(adata.X)
adata.X = sparse_X

print("saving")
adata.write_h5ad(Path('COVID.h5ad'), compression="gzip")
