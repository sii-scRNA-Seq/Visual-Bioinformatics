#! /usr/bin/bash
token=$1
version=$2
git fetch
git checkout tags/$version
cd front-end
source ~/.bashrc
npm install
npm run build
mv dist ..; cd ..
mv dist back-end; cd back-end
mkdir data
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
mkdir write
conda remove --name visual-bioinformatics --all
conda env create -f environment.yml
conda activate visual-bioinformatics
