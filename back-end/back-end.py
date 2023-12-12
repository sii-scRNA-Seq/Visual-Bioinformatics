from flask import Flask, jsonify
from flask_cors import CORS
import numpy as np
import pandas as pd
import scanpy as sc
import copy
sc.settings.verbosity = 3

app = Flask(__name__)
CORS(app)

data = {
    'pbmc3k': None,
    'filtered' : None
}

@app.route('/loaddata/')
def loaddata():
    #if data['pbmc3k'] is None:
        #data['pbmc3k'] = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True) 
        #data['pbmc3k'].var_names_make_unique()
    message = {
        'text': str(data['pbmc3k']),
        'other': ''
    }
    return jsonify(message)

@app.route('/basicfiltering/<min_genes>/<min_cells>/')
def basicfiltering(min_genes, min_cells):
    if data['pbmc3k'] is None:
        raise Exception('ERROR: Cannot complete request')
    else:
        data['filtered'] = copy.copy(data['pbmc3k'])
        sc.pp.filter_cells(data['filtered'], min_genes = int(min_genes))
        sc.pp.filter_genes(data['filtered'], min_cells = int(min_cells))
        message = {
            'text': str(data['filtered']),
            'other': ''
        }
        return jsonify(message)

if __name__ == '__main__':
    app.run()
