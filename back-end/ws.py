from anndata import AnnData
from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
from flask_websockets import WebSockets, ws
from matplotlib import pyplot as plt
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

    logger.info("Loading raw data...")
    raw_data_cache = {
        "pbmc3k": sc.datasets.pbmc3k()
    }
    logger.info("Finished loading raw data")

    app = Flask(__name__, static_folder='dist/visual-bioinformatics', static_url_path='/dist/visual-bioinformatics')

    accepting_user_requests = Cache(config=user_cache_config)
    accepting_user_requests.init_app(app)

    CORS(app)
    sockets = WebSockets(app)

    @app.route('/api/getuserid')
    def get_user_id():
        user_id = request.args.get('user_id')
        if user_id is None:
            raise ValueError
        else:
            if user_id == '':
                user_id = str(uuid.uuid4())
                logger.info(f"created user_id={user_id}")
            message = {
                'user_id': user_id
            }
            return jsonify(message)

    class UserIDException(Exception):
        pass

    class NotAcceptingRequestException(Exception):
        pass

    @sockets.on_message
    def execute_blocks(message_json):
        try:
            message = json.loads(message_json)
            user_id = message['user_id']
            logger.info(f"user_id={user_id}")
            if user_id is None or user_id == '':
                raise UserIDException
            elif accepting_user_requests.get(user_id) is False:
                raise NotAcceptingRequestException
            user_data = None
            accepting_user_requests.set(user_id, False)
            blocks = message['blocks']
            for block in blocks:
                if block['block_id'] == 'loaddata':
                    user_data, output_message = load_data(user_data, block)
                elif block['block_id'] == 'basicfiltering':
                    user_data, output_message = basic_filtering(user_data, block)
                elif block['block_id'] == 'qcplots':
                    user_data, output_message = qc_plots(user_data, block)
                else:
                    raise ValueError
                ws.send(json.dumps(output_message))
            end_connection = json.dumps({'end_connection': 'end_connection'})
        except UserIDException:
            end_connection = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
        except NotAcceptingRequestException:
            end_connection = json.dumps({'error': 'You have another request in progress, please wait and try again'})
        except KeyError:
            end_connection = json.dumps({'error': 'Received an incomplete request, please refresh the page and try again'})
        except ValueError:
            end_connection = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
        except Exception:
            end_connection = json.dumps({'error': 'Unknown error, please refresh the page and try again'})
        finally:
            ws.send(end_connection)
            accepting_user_requests.set(user_id, True)

    def load_data(user_data, block):
        logger.info("")
        new_adata = raw_data_cache.get('pbmc3k').copy()
        message = {
            'text': adata_text(new_adata)
        }
        return new_adata, message

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

    def adata_text(adata: AnnData) -> str:
        return f'Object with: {adata.n_obs:,} cells and {adata.n_vars:,} genes'

    @app.route('/<path:path>')
    def serve_spa_files(path):
        return app.send_static_file(path)

    @app.route('/')
    def serve_spa_default():
        return app.send_static_file('index.html')

    return app


app = create_app()

if __name__ == '__main__':
    from gevent import pywsgi
    from geventwebsocket.handler import WebSocketHandler
    server = pywsgi.WSGIServer(('127.0.0.1', 5000), app, handler_class=WebSocketHandler)
    server.serve_forever()
