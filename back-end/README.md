From this directory, in a Conda command line, run "conda env -f create environment.yml" to create an evironment called "visual-bioinformatics".
Activate this environment using "conda activate visual-bioinformatics"
Run "python back-end.py" to launch the server on http://127.0.0.1:5000/

Use this code from https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html to download the data:
# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write