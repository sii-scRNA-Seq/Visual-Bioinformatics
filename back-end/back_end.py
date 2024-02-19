import base64
from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
import copy
import json
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
import uuid
import werkzeug.exceptions as we
import codecs
import io

sc.settings.verbosity = 3
plt.switch_backend('agg')
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
                    'working_data': None
                })
            message = {
                'user_id': user_id
            }
            return jsonify(message)

    @app.route('/loaddata')
    def load_data():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '':
            raise we.BadRequest('Not a valid user_id')
        else:
            # if user_cache.get(user_id) is None:
            #     user_cache.set(user_id, {
            #         'basic_filtering': (None, None),
            #         'qc_plots': (None, None),
            #         # Reset cache
            #     })
            if raw_data_cache.get('pbmc3k') is None:
                data = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
                data.var_names_make_unique()
                raw_data_cache.set('pbmc3k', data)
            user_cache.set(user_id, {
                'working_data': copy.copy(raw_data_cache.get('pbmc3k')),
            })
            message = {
                'text': str(user_cache.get(user_id)['working_data']),
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
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            min_genes = int(request.args.get('min_genes'))
            min_cells = int(request.args.get('min_cells'))
            new_adata = copy.copy(user_cache.get(user_id)['working_data'])
            sc.pp.filter_cells(new_adata, min_genes=min_genes)
            sc.pp.filter_genes(new_adata, min_cells=min_cells)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            message = {
                'text': str(user_cache.get(user_id)['working_data']),
            }
            return jsonify(message)

    @app.route('/qcplots')
    def qc_plots():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            new_adata = copy.copy(user_cache.get(user_id)['working_data'])
            new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            sc.pl.violin(new_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
            my_stringIObytes = io.BytesIO()
            plt.savefig(my_stringIObytes, format='png')
            my_stringIObytes.seek(0)
            message = {
                'img': str(codecs.encode(my_stringIObytes.read(), 'base64'))
            }
            return jsonify(message)

    @app.route('/qcfiltering')
    def qc_filtering():
        user_id = request.args.get('user_id')
        invalid_params = get_invalid_parameters(['n_genes_by_counts', 'pct_counts_mt'])
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            n_genes_by_counts = int(request.args.get('n_genes_by_counts'))
            pct_counts_mt = int(request.args.get('pct_counts_mt'))
            new_adata = copy.copy(user_cache.get(user_id)['working_data'])
            new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            new_adata = new_adata[new_adata.obs.n_genes_by_counts < n_genes_by_counts, :]
            new_adata = new_adata[new_adata.obs.pct_counts_mt < pct_counts_mt, :]
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            message = {
                'text': str(user_cache.get(user_id)['working_data']),
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
