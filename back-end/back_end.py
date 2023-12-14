import json
from flask import Flask, jsonify, request
from flask_cors import CORS
import werkzeug.exceptions
import numpy as np
import pandas as pd
import scanpy as sc
import copy
sc.settings.verbosity = 3


def create_app():
    
    app = Flask(__name__)
    CORS(app)

    data = {
        'pbmc3k': None,
        'filtered' : None
    }

    class IncorrectOrderException(werkzeug.exceptions.HTTPException):
        code = 406
        description = 'The blocks you have executed are not a valid order. Please check the order and try again.'

    def handle_exception(e):
        """Return JSON instead of HTML for HTTP errors."""
        response = e.get_response()
        response.data = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.content_type = "application/json"
        return response

    app.register_error_handler(IncorrectOrderException, handle_exception)
    app.register_error_handler(werkzeug.exceptions.BadRequest, handle_exception)

    @app.route('/loaddata')
    def loaddata():
        if data['pbmc3k'] is None:
            data['pbmc3k'] = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True) 
            data['pbmc3k'].var_names_make_unique()
        message = {
            'text': str(data['pbmc3k']),
            'other': ''
        }
        return jsonify(message)

    @app.route('/basicfiltering')
    def basicfiltering():
        invalid_params = getInvalidParameters(['min_genes','min_cells'])
        if invalid_params != []:
            raise werkzeug.exceptions.BadRequest('Missing parameters: ' + str(invalid_params))
        elif data['pbmc3k'] is None:
            raise IncorrectOrderException()
        else:
            min_genes = request.args.get('min_genes')
            min_cells = request.args.get('min_cells')
            data['filtered'] = copy.copy(data['pbmc3k'])
            sc.pp.filter_cells(data['filtered'], min_genes = int(min_genes))
            sc.pp.filter_genes(data['filtered'], min_cells = int(min_cells))
            message = {
                'text': str(data['filtered']),
                'other': ''
            }
            return jsonify(message)
    
    def getInvalidParameters(params):
        invalid_params = []
        for param in params:
            if request.args.get(param) is None:
                invalid_params.append(param)
        return invalid_params

    return app


app = create_app()

if __name__ == '__main__':

    app.run()
