from pathlib import Path
import scanpy as sc

print("reading data")
adata = sc.datasets.pbmc3k()
adata.obs["sample"] = "1"

print("saving")
adata.write_h5ad(Path('PBMC3K.h5ad'), compression="gzip")
