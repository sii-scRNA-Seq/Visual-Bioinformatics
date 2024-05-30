# download_pk_mca_fastqs.py
# Run me from project root (i.e. /users/2117532m/home/SCAMPI_PF_dataset)

# $ module load apps/python3
# $ pip3 install pandas

import requests
import json
from pathlib import Path
import os
import urllib.request
import pandas as pd

extended_metadata = pd.read_csv('PRJEB58790.tsv', header=0, sep="\t")

extended_metadata = extended_metadata[extended_metadata["Sanger_SampleName"].str.contains("Day")]
extended_metadata = extended_metadata[extended_metadata["Short_or_Long_Read"].str.contains("short read")]

extended_metadata_obj_list = extended_metadata.to_dict(orient='records')

root_dir = "MCA_data"
Path(root_dir).mkdir(parents=True, exist_ok=True)

for assay in extended_metadata_obj_list:
    print(f"Found {assay['Sanger_SampleName']} - {assay['run_accession']}")
    assay_path = f"{root_dir}/{assay['run_accession']}"
    Path(assay_path).mkdir(parents=True, exist_ok=True)

    # Clear old downloads which may have been interrupted/corrupt
    for filename in os.listdir(assay_path):
        file_path = os.path.join(assay_path, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.remove(file_path)

    fastq_urls = assay['fastq_ftp'].split(";")
    assert len(fastq_urls) == 2
    print(fastq_urls)

    for fastq in fastq_urls:
        fastq_filepath = os.path.join(assay_path, fastq.split("/")[-1])
        print(f"saving {fastq_filepath}")
        urllib.request.urlretrieve(f"ftp://{fastq}", fastq_filepath)

print("DONE")

