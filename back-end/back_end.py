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
THREE_DAYS = 3 * 24 * 60 * 60


def create_app(test_mode=False):

    if test_mode:
        config = {
            "CACHE_TYPE": "SimpleCache",
            "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
        }
    else:
        config = {
            "CACHE_TYPE": "FileSystemCache",
            "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
            "CACHE_IGNORE_ERRORS": False,  # Default
            "CACHE_DIR": 'back-end-cache',
            "CACHE_THRESHOLD": 500,        # Default
        }
    app = Flask(__name__)
    app.config.from_mapping(config)
    user_cache = Cache(app)
    raw_data_cache = Cache(app)
    CORS(app)

    class IncorrectOrderException(we.HTTPException):
        code = 406
        description = ('The blocks you have executed are not a valid order. Please check the order and try again.')

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
        if user_id is None:
            raise we.BadRequest('Not a valid user_id')
        else:
            if user_id == '':
                user_id = str(uuid.uuid4())
            if user_cache.get(user_id) is None:
                user_cache.set(user_id, {
                    'basic_filtering': (None, None),
                    'qc_plots': (None, None),
                    # Reset cache
                })
            message = {
                'text': user_id
            }
            return jsonify(message)

    @app.route('/loaddata')
    def load_data():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '':
            raise we.BadRequest('Not a valid user_id')
        else:
            if user_cache.get(user_id) is None:
                user_cache.set(user_id, {
                    'basic_filtering': (None, None),
                    'qc_plots': (None, None),
                    # Reset cache
                })
            if raw_data_cache.get('pbmc3k') is None:
                data = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
                data.var_names_make_unique()
                raw_data_cache.set('pbmc3k', data)
            message = {
                'text': str(raw_data_cache.get('pbmc3k')),
            }
            return jsonify(message)

    @app.route('/basicfiltering')
    def basic_filtering():
        user_id = request.args.get('user_id')
        invalid_params = get_invalid_parameters(['min_genes', 'min_cells'])
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif raw_data_cache.get('pbmc3k') is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            min_genes = int(request.args.get('min_genes'))
            min_cells = int(request.args.get('min_cells'))
            if user_cache.get(user_id)['basic_filtering'][0] != (min_genes, min_cells):
                new_adata = copy.copy(raw_data_cache.get('pbmc3k'))
                sc.pp.filter_cells(new_adata, min_genes=min_genes)
                sc.pp.filter_genes(new_adata, min_cells=min_cells)
                user_cache.set(user_id, {
                    'basic_filtering': ((min_genes, min_cells), new_adata),
                    'qc_plots': (None, None),
                    # Reset cache
                })
            message = {
                'text': str(user_cache.get(user_id)['basic_filtering'][1]),
            }
            return jsonify(message)

    @app.route('/qcplots')
    def qc_plots():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif user_cache.get(user_id)['basic_filtering'][1] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            if user_cache.get(user_id)['qc_plots'][0] != ():
                new_adata = copy.copy(user_cache.get(user_id)['basic_filtering'][1])
                new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
                sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                # sc.pl.violin(new_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
                user_cache.set(user_id, {
                    'basic_filtering': user_cache.get(user_id)['basic_filtering'],
                    'qc_plots': ((), new_adata),
                    # Reset cache
                })
            message = {
                'text': str(user_cache.get(user_id)['qc_plots'][1]),
            }
            return jsonify(message)

    @app.route('/qcfiltering')
    def qc_filtering():
        # TODO: Implement this
        message = {
            'text': 'Implement this',
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
