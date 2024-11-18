# Covid dataset
This dataset taken from the lab reference dataset generated here: https://www.science.org/doi/10.1126/science.adj4088

## To recreate:

### Requirements
1. Conda environment (use ./environment.yml) with:
   * Python3 with pandas and scanpy
   * R with seurat

### Steps
* Create conda env: `conda env create -f ./environment.yml`
* Activate conda env `conda activate scampi-generate-covid-dataset`
* Use `RHOME="" Rscript generate_matrix.r` to generate a count matrix
* Use `python create_object.py` to create the h5ad object.
