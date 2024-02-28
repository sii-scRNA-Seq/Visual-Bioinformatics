from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
from matplotlib import pyplot as plt
import codecs
import copy
import io
import json
# import numpy as np
# import pandas as pd
import scanpy as sc
import uuid
import werkzeug.exceptions as we
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
        description = ('The blocks you have executed are not a valid order. Please check the blocks and try again.')

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
            # To be added once user cache is fully implemented
            # if user_cache.get(user_id) is None:
            #     user_cache.set(user_id, {
            #         'basic_filtering': (None, None),
            #         'qc_plots': (None, None),
            #         'qc_filtering': (None, None),
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
            min_genes = float(request.args.get('min_genes'))
            min_cells = float(request.args.get('min_cells'))
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
            plt.rcParams['font.size'] = 18
            sc.pl.violin(new_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A violin plot displaying quality control metrics generated by a QC Plots block',
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
            n_genes_by_counts = float(request.args.get('n_genes_by_counts'))
            pct_counts_mt = float(request.args.get('pct_counts_mt'))
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

    @app.route('/variablegenes')
    def variable_genes():
        user_id = request.args.get('user_id')
        invalid_params = get_invalid_parameters(['min_mean', 'max_mean', 'min_disp'])
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            min_mean = float(request.args.get('min_mean'))
            max_mean = float(request.args.get('max_mean'))
            min_disp = float(request.args.get('min_disp'))
            new_adata = copy.copy(user_cache.get(user_id)['working_data'])
            sc.pp.normalize_total(new_adata, target_sum=1e4)
            sc.pp.log1p(new_adata)
            sc.pp.highly_variable_genes(new_adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            plt.rcParams['font.size'] = 14
            sc.pl.highly_variable_genes(new_adata)
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A scatter plot displaying dispersions of genes generated by an Identify Highly Variable Genes block',
            }
            return jsonify(message)

    @app.route('/pca')
    def pca():
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
            sc.pp.regress_out(new_adata, ['total_counts', 'pct_counts_mt'])
            sc.pp.scale(new_adata, max_value=10)
            sc.tl.pca(new_adata, svd_solver='arpack')
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            sc.pl.pca_variance_ratio(new_adata, log=True)
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A scatter plot displaying the contribution of each PC to the total variance in the data, generated by a Principle Component Analysis block',
            }
            return jsonify(message)

    @app.route('/runumap')
    def run_umap():
        user_id = request.args.get('user_id')
        invalid_params = get_invalid_parameters(['n_neighbors', 'n_pcs'])
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif not (request.args.get('n_neighbors').isdecimal() and request.args.get('n_pcs').isdecimal()):
            raise we.BadRequest('Number of Neighbors and Number of Principle Components must both be integers')
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            n_neighbors = int(request.args.get('n_neighbors'))
            n_pcs = int(request.args.get('n_pcs'))
            new_adata = copy.copy(user_cache.get(user_id)['working_data'])
            sc.pp.neighbors(new_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.umap(new_adata)
            sc.tl.leiden(new_adata)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            sc.pl.umap(new_adata, color=['leiden'])
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A UMAP of Leiden clusters using the principle components generated by a Run UMAP block',
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
