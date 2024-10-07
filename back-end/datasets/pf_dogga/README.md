# PF dataset
This dataset taken from the lab reference dataset generated here: https://www.science.org/doi/10.1126/science.adj4088

## To recreate:

### Requirements
1. Python3 with pandas and scanpy
2. Cellranger
3. AGAT for conversion between gff and gtf
4. 300-400GB spare hdd space

### Steps
* Use `download_fastqs.py` to download fastq files
    * You will need to rename the fastq files to match the cellranger specification, for example: 
    * `ERR11471992_1.fastq.gz` -> `ERR11471992_S1_L001_R1_001.fastq.gz`
    * `ERR11471992_2.fastq.gz` -> `ERR11471992_S1_L001_R2_001.fastq.gz`
* Use `make_reference.sh` to make a cellranger reference from the reference we have from PlasmoDb.
* Use `1_map_DayX.sh` to map each day
* Use `create_object.py` to create the h5ad object.
