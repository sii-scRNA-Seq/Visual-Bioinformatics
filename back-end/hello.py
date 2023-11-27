from flask import Flask, jsonify
from flask_cors import CORS
import numpy as np
import pandas as pd
import scanpy as sc

app = Flask(__name__)
CORS(app)

@app.route('/')
def get_data():
    data = {
        'text': 'Hello, World!',
        'other': "You can't see this bit"
    }
    return jsonify(data)

@app.route('/loaddata/')
def loaddata():
    adata = sc.read_10x_mtx(
    'C:/Users/rober/OneDrive/Desktop/Bioinformatics Project/scanpy-tutorials-master/data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols', cache=True) 
    adata.var_names_make_unique()
    message = {
        'text': str(adata),
        'other': ''
    }
    return jsonify(message)

if __name__ == '__main__':
    app.run()