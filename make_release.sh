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
conda remove --name visual-bioinformatics --all
conda env create -f environment.yml
conda activate visual-bioinformatics
