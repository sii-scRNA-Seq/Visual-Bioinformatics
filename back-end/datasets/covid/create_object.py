from pathlib import Path
import scanpy as sc
import pandas as pd
from scipy import sparse

print("reading data")
meta = pd.read_csv("./temp/covid_meta.csv")

# Convert strings to categories to reduce file size
for col in ['Donor', 'Status']:
    meta[col] = meta[col].astype('category')

adata = sc.read_csv(Path("./temp/covid_matrix.csv"))
adata.obs = meta
print(adata)

print("Compressing matrix")
sparse_X = sparse.csr_matrix(adata.X)
adata.X = sparse_X

print("saving")
adata.write_h5ad(Path('COVID.h5ad'), compression="gzip")
