from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
import copy
import json
import numpy as np
import pandas as pd
import scanpy as sc
import uuid
import werkzeug.exceptions as we
sc.settings.verbosity = 3


def create_app():

    config = {
        "CACHE_TYPE": "SimpleCache",
        "CACHE_DEFAULT_TIMEOUT": 3*24*60*60
    }
    app = Flask(__name__)
    app.config.from_mapping(config)
    cache = Cache(app)
    CORS(app)

    data = {
        'pbmc3k': None,
        'filtered': None
    }

    class IncorrectOrderException(we.HTTPException):
        code = 406
        description = ('The blocks you have executed are not a valid order. '
                       'Please check the order and try again.')

    def handle_exception(e):
        response = e.get_response()
        response.data = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.content_type = "application/json"
        return response

    app.register_error_handler(IncorrectOrderException, handle_exception)
    app.register_error_handler(we.BadRequest, handle_exception)

    @app.route('/useridisvalid')
    def user_id_is_valid():
        user_id = request.args.get('user_id')
        message = {
            'text': str(cache.get(user_id) is not None)
        }
        return jsonify(message)

    @app.route('/createuserid')
    def create_user_id():
        user_id = str(uuid.uuid4())
        cache.set(user_id, {})
        print(cache.get(user_id))
        message = {
            'text': user_id
        }
        return jsonify(message)

    @app.route('/loaddata')
    def load_data():
        if data['pbmc3k'] is None:
            data['pbmc3k'] = sc.read_10x_mtx(
                'data/filtered_gene_bc_matrices/hg19/',
                var_names='gene_symbols',
                cache=True)
            data['pbmc3k'].var_names_make_unique()
        message = {
            'text': str(data['pbmc3k']),
        }
        return jsonify(message)

    @app.route('/basicfiltering')
    def basic_filtering():
        invalid_params = get_invalid_parameters(['min_genes', 'min_cells'])
        if invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif data['pbmc3k'] is None:
            raise IncorrectOrderException()
        else:
            min_genes = request.args.get('min_genes')
            min_cells = request.args.get('min_cells')
            data['filtered'] = copy.copy(data['pbmc3k'])
            sc.pp.filter_cells(data['filtered'], min_genes=int(min_genes))
            sc.pp.filter_genes(data['filtered'], min_cells=int(min_cells))
            message = {
                'text': str(data['filtered']),
            }
            return jsonify(message)

    def get_invalid_parameters(params):
        invalid_params = []
        for param in params:
            if request.args.get(param) is None:
                invalid_params.append(param)
        return invalid_params

    return app


app = create_app()

if __name__ == '__main__':

    app.run()
