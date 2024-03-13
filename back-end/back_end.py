from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
from matplotlib import pyplot as plt
import codecs
import io
import json
from anndata import AnnData
import scanpy as sc
import uuid
import werkzeug.exceptions as we
from threadpoolctl import threadpool_limits
import logging
import logging.config
import yaml

sc.settings.verbosity = 0
plt.switch_backend('agg')

THREE_DAYS = 3 * 24 * 60 * 60

with open('logging.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)


def create_app(test_mode=False):

    if test_mode:
        logger = logging.getLogger('scampi-test')
        config = {
            "CACHE_TYPE": "SimpleCache",
            "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
        }
    else:
        logger = logging.getLogger('scampi')
        config = {
            "CACHE_TYPE": "FileSystemCache",
            "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
            "CACHE_IGNORE_ERRORS": False,  # Default
            "CACHE_DIR": 'back-end-cache',
            "CACHE_THRESHOLD": 500,        # Default
        }

    external_loggers = [
        logging.getLogger('werkzeug'),
        logging.getLogger('waitress'),
        logging.getLogger('scanpy')
    ]
    for external_logger in external_loggers:
        external_logger.handlers.clear()
        external_logger.handlers.extend(logger.handlers)
        external_logger.setLevel(logging.WARNING)

    logger.info("Starting App")
    app = Flask(__name__, static_folder='dist/visual-bioinformatics', static_url_path='/dist/visual-bioinformatics')
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

    @app.route('/api/getuserid')
    def get_user_id():
        user_id = request.args.get('user_id')
        if user_id is None:
            raise we.BadRequest('Not a valid user_id')
        else:
            if user_id == '':
                user_id = str(uuid.uuid4())
                logger.info(f"created user_id={user_id}")
            if user_cache.get(user_id) is None:
                user_cache.set(user_id, {
                    'working_data': None
                })
            message = {
                'user_id': user_id
            }
            return jsonify(message)

    @app.route('/api/loaddata')
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
            logger.info(f"user_id={user_id}")

            if raw_data_cache.get('pbmc3k') is None:
                data = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)
                data.var_names_make_unique()
                raw_data_cache.set('pbmc3k', data)
            user_cache.set(user_id, {
                'working_data': raw_data_cache.get('pbmc3k').copy(),
            })
            message = {
                'text': adata_text(user_cache.get(user_id)['working_data']),
            }
            return jsonify(message)

    @app.route('/api/basicfiltering')
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
            logger.info(f"min_genes={min_genes} min_cells={min_cells} user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()
            sc.pp.filter_cells(new_adata, min_genes=min_genes)
            sc.pp.filter_genes(new_adata, min_cells=min_cells)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            message = {
                'text': adata_text(user_cache.get(user_id)['working_data']),
            }
            return jsonify(message)

    @app.route('/api/qcplots')
    def qc_plots():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            logger.info(f"user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()

            new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

            # Rename total_counts to total_UMIs
            new_adata.obs['total_UMIs'] = new_adata.obs['total_counts']
            new_adata.obs = new_adata.obs.drop('total_counts', axis=1)

            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            plt.rcParams['font.size'] = 18

            sc.pl.violin(new_adata, ['n_genes_by_counts', 'total_UMIs', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A violin plot displaying quality control metrics generated by a QC Plots block',
            }
            return jsonify(message)

    @app.route('/api/qcfiltering')
    def qc_filtering():
        user_id = request.args.get('user_id')
        invalid_params = get_invalid_parameters(['min_n_genes_by_counts', 'max_n_genes_by_counts', 'pct_counts_mt'])
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif invalid_params != []:
            raise we.BadRequest('Missing parameters: ' + str(invalid_params))
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')
            min_n_genes_by_counts = float(request.args.get('min_n_genes_by_counts'))
            max_n_genes_by_counts = float(request.args.get('max_n_genes_by_counts'))
            pct_counts_mt = float(request.args.get('pct_counts_mt'))
            logger.info(f"min_n_genes_by_counts={min_n_genes_by_counts} max_n_genes_by_counts={max_n_genes_by_counts} pct_counts_mt={pct_counts_mt} user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()
            new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

            # Rename total_counts to total_UMIs
            new_adata.obs['total_UMIs'] = new_adata.obs['total_counts']
            new_adata.obs = new_adata.obs.drop('total_counts', axis=1)

            new_adata = new_adata[new_adata.obs.n_genes_by_counts < max_n_genes_by_counts, :]
            new_adata = new_adata[new_adata.obs.n_genes_by_counts > min_n_genes_by_counts, :]
            new_adata = new_adata[new_adata.obs.pct_counts_mt < pct_counts_mt, :]
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            message = {
                'text': adata_text(user_cache.get(user_id)['working_data']),
            }
            return jsonify(message)

    @app.route('/api/variablegenes')
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
            logger.info(f"min_mean={min_mean} max_mean={max_mean} min_disp={min_disp} user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()
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

    @app.route('/api/pca')
    def pca():
        user_id = request.args.get('user_id')
        if user_id is None or user_id == '' or user_cache.get(user_id) is None:
            raise we.BadRequest('Not a valid user_id')
        elif user_cache.get(user_id)["working_data"] is None:
            raise IncorrectOrderException()
        else:
            user_id = request.args.get('user_id')

            logger.info(f"user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()
            new_adata.var['mt'] = new_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

            # Rename total_counts to total_UMIs
            new_adata.obs['total_UMIs'] = new_adata.obs['total_counts']
            new_adata.obs = new_adata.obs.drop('total_counts', axis=1)

            sc.pp.regress_out(new_adata, ['total_UMIs', 'pct_counts_mt'], n_jobs=1)
            sc.pp.scale(new_adata, max_value=10)

            with threadpool_limits(limits=1, user_api='blas'):
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

    @app.route('/api/runumap')
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
            logger.info(f"n_neighbors={n_neighbors} n_pcs={n_pcs} user_id={user_id}")

            new_adata = user_cache.get(user_id)['working_data'].copy()
            sc.pp.neighbors(new_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.umap(new_adata)
            sc.tl.leiden(new_adata)
            user_cache.set(user_id, {
                'working_data': new_adata,
            })
            sc.pl.umap(new_adata, color=['leiden'], legend_loc='on data')
            image_stream = io.BytesIO()
            plt.savefig(image_stream, format='png')
            image_stream.seek(0)
            message = {
                'img': str(codecs.encode(image_stream.read(), 'base64')),
                'alttext': 'A UMAP of Leiden clusters using the principle components generated by a Run UMAP block',
            }
            return jsonify(message)

    @app.route('/<path:path>')
    def serve_spa_files(path):
        print("Serve file: ", path)
        return app.send_static_file(path)

    @app.route('/')
    def serve_spa_default():
        return app.send_static_file('index.html')

    def get_invalid_parameters(params):
        invalid_params = []
        for param in params:
            if request.args.get(param) is None:
                invalid_params.append(param)
        return invalid_params

    def adata_text(adata: AnnData) -> str:
        return f'Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes'

    return app


app = create_app()

if __name__ == '__main__':

    app.run()
