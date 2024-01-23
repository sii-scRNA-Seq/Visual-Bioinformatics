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

    @app.route('/getuserid')
    def get_user_id():
        user_id = request.args.get('user_id')
        if cache.get(user_id) is None:
            user_id = str(uuid.uuid4())
            cache.set(user_id, {
                'basic_filtering': (None, None)
                # TODO: Reset cache
            })
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
            user_id = request.args.get('user_id')
            min_genes = int(request.args.get('min_genes'))
            min_cells = int(request.args.get('min_cells'))
            if cache.get(user_id)['basic_filtering'][0] != (min_genes,
                                                            min_cells):
                filtered_data = copy.copy(data['pbmc3k'])
                sc.pp.filter_cells(filtered_data, min_genes=min_genes)
                sc.pp.filter_genes(filtered_data, min_cells=min_cells)
                cache.set(user_id, {
                    'basic_filtering': ((min_genes, min_cells), filtered_data)
                    # TODO: Reset cache
                })
            message = {
                'text': str(cache.get(user_id)['basic_filtering'][1]),
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
