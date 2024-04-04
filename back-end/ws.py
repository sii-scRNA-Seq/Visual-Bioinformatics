from gevent import monkey
import gevent
from joblib import parallel_backend
monkey.patch_all()

from anndata import AnnData
from flask import Flask, jsonify, request, render_template
from flask_caching import Cache
from flask_cors import CORS
from flask_socketio import SocketIO
from matplotlib import pyplot as plt
from threadpoolctl import threadpool_limits
import codecs
import io
import json
import logging
import logging.config
import scanpy as sc
import uuid
import yaml

sc.settings.verbosity = 0
plt.switch_backend('agg')

THREE_DAYS = 3 * 24 * 60 * 60

with open('logging.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)

simple_cache_config = {
    "CACHE_TYPE": "SimpleCache",
    "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
}

file_cache_config = {
    "CACHE_TYPE": "FileSystemCache",
    "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
    "CACHE_IGNORE_ERRORS": False,  # Default
    "CACHE_DIR": 'back-end-cache',
    "CACHE_THRESHOLD": 500,        # Default
}


def create_app(test_mode=False):

    if test_mode:
        logger = logging.getLogger('scampi-test')
        user_cache_config = simple_cache_config
    else:
        logger = logging.getLogger('scampi')
        user_cache_config = file_cache_config

    logger.debug("DEBUG logs enabled")

    external_loggers = [
        logging.getLogger('werkzeug'),
        logging.getLogger('waitress'),
        logging.getLogger('scanpy'),
        logging.getLogger('geventwebsocket.handler')
    ]
    for external_logger in external_loggers:
        external_logger.handlers.clear()
        external_logger.handlers.extend(logger.handlers)
        external_logger.setLevel(logging.WARNING)

    logger.info("Starting App")

    logger.info("Loading raw data...")
    raw_data_cache = {
        "pbmc3k": sc.datasets.pbmc3k()
    }
    logger.info("Finished loading raw data")

    app = Flask(__name__, static_folder='dist/visual-bioinformatics', static_url_path='/dist/visual-bioinformatics')
    # TODO: I don't know if this is required but it was in all the flask_socketio examples
    app.config['SECRET_KEY'] = 'secret!'

    accepting_user_requests = Cache(config=user_cache_config)
    accepting_user_requests.init_app(app)

    # TODO, should we really be accepting CORS requests *all* the time?
    CORS(app)
    socketio = SocketIO(app, cors_allowed_origins="*")

    logger.info(f"Server async mode: {socketio.server.eio.async_mode}")

    @app.route('/')
    def index():
        return render_template('index.html')

    @app.route('/api/getuserid')
    def get_user_id():
        user_id = request.args.get('user_id')
        if user_id is None:
            raise ValueError
        else:
            if user_id == '':
                user_id = str(uuid.uuid4())
                logger.info(f"created user_id={user_id}")
            else:
                logger.info(f"User ID kept user_id={user_id}")
            message = {
                'user_id': user_id
            }
            return jsonify(message)

    class UserIDException(Exception):
        pass

    class NotAcceptingRequestException(Exception):
        pass

    @socketio.on('json')
    def execute_blocks(message):
        logger.info(f"Executing blocks, json={message}")
        try:
            user_id = message['user_id']
            if user_id is None or user_id == '':
                logger.error("User id not in message")
                raise UserIDException
            elif accepting_user_requests.get(user_id) is False:
                logger.error("User already processing blocks")
                raise NotAcceptingRequestException
            logger.info(f"user_id={user_id}")
            user_data = None
            accepting_user_requests.set(user_id, False)
            blocks = message['blocks']
            for block in blocks:
                logger.debug(f"Executing block={block} user={user_id}")
                if block['block_id'] == 'loaddata':
                    user_data, output_message = load_data(user_data, block)
                elif block['block_id'] == 'basicfiltering':
                    user_data, output_message = basic_filtering(user_data, block)
                elif block['block_id'] == 'qcplots':
                    user_data, output_message = qc_plots(user_data, block)
                elif block['block_id'] == 'qcfiltering':
                    user_data, output_message = qc_filtering(user_data, block)
                elif block['block_id'] == 'variablegenes':
                    user_data, output_message = variable_genes(user_data, block)
                elif block['block_id'] == 'pca':
                    user_data, output_message = pca(user_data, block)
                elif block['block_id'] == 'runumap':
                    user_data, output_message = run_umap(user_data, block)
                else:
                    raise ValueError
                socketio.emit("json", json.dumps(output_message))
                logger.debug('emitted:' + json.dumps(output_message))
                # Allow other threads to execute
                #https://stackoverflow.com/questions/30901998/threading-true-with-flask-socketio
                gevent.sleep()
            logger.debug(f"Finished processing blocks for user={user_id}")
            end_connection = json.dumps({'end_connection': 'end_connection'})
        except UserIDException as e:
            logger.error(e, exc_info=True)
            end_connection = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
        except NotAcceptingRequestException as e:
            logger.error(e, exc_info=True)
            end_connection = json.dumps({'error': 'You have another request in progress, please wait and try again'})
        except KeyError as e:
            logger.error(e, exc_info=True)
            end_connection = json.dumps({'error': 'Received an incomplete request, please refresh the page and try again'})
        except ValueError as e:
            logger.error(e, exc_info=True)
            end_connection = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
        except Exception as e:
            logger.error(e, exc_info=True)
            end_connection = json.dumps({'error': 'Unknown error, please refresh the page and try again'})
        finally:
            socketio.emit("json", end_connection)
            accepting_user_requests.set(user_id, True)

    def load_data(user_data, block):
        logger.info("")
        user_data = raw_data_cache.get('pbmc3k').copy()
        logger.debug(user_data)
        message = {
            'text': adata_text(user_data)
        }
        return user_data, message

    def basic_filtering(user_data, block):
        min_genes = float(block['min_genes'])
        min_cells = float(block['min_cells'])
        logger.info(f"min_genes={min_genes} min_cells={min_cells}")

        sc.pp.filter_cells(user_data, min_genes=min_genes)
        sc.pp.filter_genes(user_data, min_cells=min_cells)

        message = {
            'text': adata_text(user_data)
        }
        return user_data, message

    def qc_plots(user_data, block):
        logger.info("")

        user_data.var['mt'] = user_data.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(user_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        user_data.obs['total_UMIs'] = user_data.obs['total_counts']
        user_data.obs = user_data.obs.drop('total_counts', axis=1)
        plt.rcParams['font.size'] = 18
        sc.pl.violin(user_data, ['n_genes_by_counts', 'total_UMIs', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
        
        image_stream = io.BytesIO()
        plt.savefig(image_stream, format='png')
        image_stream.seek(0)
        message = {
            'img': str(codecs.encode(image_stream.read(), 'base64')),
            'alttext': 'A violin plot displaying quality control metrics generated by a QC Plots block',
        }
        return user_data, message

    def qc_filtering(user_data, block):
        min_n_genes_by_counts = float(block['min_n_genes_by_counts'])
        max_n_genes_by_counts = float(block['max_n_genes_by_counts'])
        pct_counts_mt = float(block['pct_counts_mt'])
        logger.info(f"min_n_genes_by_counts={min_n_genes_by_counts} max_n_genes_by_counts={max_n_genes_by_counts} pct_counts_mt={pct_counts_mt}")  
        
        user_data.var['mt'] = user_data.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(user_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        user_data.obs['total_UMIs'] = user_data.obs['total_counts']
        user_data.obs = user_data.obs.drop('total_counts', axis=1)
        user_data = user_data[user_data.obs.n_genes_by_counts < max_n_genes_by_counts, :]
        user_data = user_data[user_data.obs.n_genes_by_counts > min_n_genes_by_counts, :]
        user_data = user_data[user_data.obs.pct_counts_mt < pct_counts_mt, :]
        
        message = {
            'text': adata_text(user_data)
        }
        return user_data, message
    
    def variable_genes(user_data, block):
        min_mean = float(block['min_mean'])
        max_mean = float(block['max_mean'])
        min_disp = float(block['min_disp'])
        logger.info(f"min_mean={min_mean} max_mean={max_mean} min_disp={min_disp}")  
        
        sc.pp.normalize_total(user_data, target_sum=1e4)
        sc.pp.log1p(user_data)
        sc.pp.highly_variable_genes(user_data, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
        plt.rcParams['font.size'] = 14
        sc.pl.highly_variable_genes(user_data)
        
        image_stream = io.BytesIO()
        plt.savefig(image_stream, format='png')
        image_stream.seek(0)
        message = {
            'img': str(codecs.encode(image_stream.read(), 'base64')),
            'alttext': 'A scatter plot displaying dispersions of genes generated by an Identify Highly Variable Genes block',
        }
        return user_data, message
    
    def pca(user_data, block):
        logger.info("")
        
        user_data.var['mt'] = user_data.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(user_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        user_data.obs['total_UMIs'] = user_data.obs['total_counts']
        user_data.obs = user_data.obs.drop('total_counts', axis=1)

        # TODO: Why is this breaking the websocket?
        # with parallel_backend('threading', n_jobs=1):
        sc.pp.regress_out(user_data, ['total_UMIs', 'pct_counts_mt'], n_jobs=1)
        sc.pp.scale(user_data, max_value=10)

        with threadpool_limits(limits=1, user_api='blas'):
            sc.tl.pca(user_data, svd_solver='arpack')
        sc.pl.pca_variance_ratio(user_data, log=True)
        
        image_stream = io.BytesIO()
        plt.savefig(image_stream, format='png')
        image_stream.seek(0)
        message = {
            'img': str(codecs.encode(image_stream.read(), 'base64')),
            'alttext': 'A scatter plot displaying the contribution of each PC to the total variance in the data, generated by a Principle Component Analysis block',
        }
        return user_data, message
    
    def run_umap(user_data, block):
        n_neighbors = int(block['n_neighbors'])
        n_pcs = int(block['n_pcs'])
        logger.info(f"n_neighbors={n_neighbors} n_pcs={n_pcs}")
        
        sc.pp.neighbors(user_data, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(user_data)
        sc.tl.leiden(user_data)
        sc.pl.umap(user_data, color=['leiden'], legend_loc='on data')
        
        image_stream = io.BytesIO()
        plt.savefig(image_stream, format='png')
        image_stream.seek(0)
        message = {
            'img': str(codecs.encode(image_stream.read(), 'base64')),
            'alttext': 'A UMAP of Leiden clusters using the principle components generated by a Run UMAP block',
        }
        return user_data, message

    def adata_text(adata: AnnData) -> str:
        return f'Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes'

    @app.route('/<path:path>')
    def serve_spa_files(path):
        return app.send_static_file(path)

    @app.route('/')
    def serve_spa_default():
        return app.send_static_file('index.html')

    return socketio, app


if __name__ == '__main__':
    socketio, app = create_app(True)

    # By default, threading is handled by gevent.spawn:
    # https://github.com/miguelgrinberg/Flask-SocketIO/blob/40007fded6228013ce7e408ea1d9628da8b125fa/src/flask_socketio/__init__.py#L700C36-L700C42
    # https://www.gevent.org/api/gevent.baseserver.html#gevent.baseserver.BaseServer

    socketio.run(app)
    # from gevent import pywsgi
    # from geventwebsocket.handler import WebSocketHandler
    # server = pywsgi.WSGIServer(('127.0.0.1', 5000), app, handler_class=WebSocketHandler)
    # server.serve_forever()
