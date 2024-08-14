from gevent import monkey
monkey.patch_all()

from flask import Flask, jsonify, request
from flask_caching import Cache
from flask_cors import CORS
from flask_socketio import SocketIO

from matplotlib import pyplot as plt
import argparse
import gevent
import json
import logging
import logging.config
import scanpy as sc
import uuid
import werkzeug.exceptions as we
import yaml

from exception.bad_request_exception import BadRequestException
from block.exception.missing_param_exception import MissingParametersException
from exception.not_accepting_request_exception import NotAcceptingRequestException
from exception.user_id_exception import UserIDException

from block.basic_filtering import BasicFiltering
from block.qc_plots import QCPlots
from block.qc_filtering import QCFiltering
from block.variable_genes import VariableGenes
from block.pca import PCA
from block.integration import Integration
from block.run_umap import RunUMAP

from block.block_interface import adata_text

from dataset_info import dataset_info

sc.settings.verbosity = 0
plt.switch_backend("agg")

with open("logging.yml", "rt") as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)

THREE_DAYS = 3 * 24 * 60 * 60

simple_cache_config = {
    "CACHE_TYPE": "SimpleCache",
    "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
}

file_cache_config = {
    "CACHE_TYPE": "FileSystemCache",
    "CACHE_DEFAULT_TIMEOUT": THREE_DAYS,
    "CACHE_IGNORE_ERRORS": False,  # Default
    "CACHE_DIR": "back-end-cache",
    "CACHE_THRESHOLD": 500,        # Default
}


def create_app(test_mode=False):

    if test_mode:
        logger = logging.getLogger("scampi-test")
        user_cache_config = simple_cache_config
        logger.info("Starting in test mode")
    else:
        logger = logging.getLogger("scampi")
        user_cache_config = file_cache_config
        logger.info("Starting in production mode")

    logger.debug("DEBUG logs enabled")

    external_loggers = [
        logging.getLogger("werkzeug"),
        logging.getLogger("waitress"),
        logging.getLogger("scanpy"),
        logging.getLogger("geventwebsocket.handler"),
        logging.getLogger("ws")
    ]
    for external_logger in external_loggers:
        external_logger.handlers.clear()
        external_logger.handlers.extend(logger.handlers)
        external_logger.setLevel(logging.WARNING)

    logger.info("Starting App")

    logger.info("Loading raw data...")
    raw_data_cache = {
        "pbmc3k": sc.read_h5ad("../datasets/pbmc3k/PBMC3K.h5ad"),
        "pf_dogga": sc.read_h5ad("../datasets/pf_dogga/MCA_PF_DOGGA.h5ad")
    }
    logger.info("Finished loading raw data")

    app = Flask(__name__, static_folder="../dist/scampi/browser", static_url_path="/dist/scampi/browser")

    accepting_user_requests = Cache(config=user_cache_config)
    accepting_user_requests.init_app(app)

    # Consider permitted sources for CORS requests (see Trello)
    CORS(app)
    socketio = SocketIO(app, cors_allowed_origins="*")

    def handle_exception(e):
        response = e.get_response()
        response.data = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.content_type = "application/json"
        return response

    app.register_error_handler(we.BadRequest, handle_exception)

    logger.info(f"Server async mode: {socketio.server.eio.async_mode}")

    @app.route("/")
    def serve_spa_default():
        logger.warning("Serving default")
        return app.send_static_file("index.html")

    @app.route("/<path:path>")
    def serve_spa_files(path):
        logger.warning("Serving " + path)
        return app.send_static_file(path)

    @app.route("/api/getuserid")
    def get_user_id():
        user_id = request.args.get("user_id")
        if not isinstance(user_id, str):
            raise we.BadRequest("Not a valid user_id")
        else:
            if user_id == "":
                user_id = str(uuid.uuid4())
                logger.info(f"created user_id={user_id}")
            else:
                logger.info(f"User ID kept user_id={user_id}")
            message = {
                "user_id": user_id
            }
            return jsonify(message)

    @app.route("/api/getdatasetinfo")
    def get_dataset_info():
        logger.info("Sending dataset info")
        message = {
            "datasets": dataset_info
        }
        return jsonify(message)

    @socketio.on("json")
    def execute_blocks(message):
        logger.info(f"Executing blocks, json={message}")

        try:

            client = request.sid

            if "user_id" not in message:
                raise UserIDException("User ID is missing")
            elif not isinstance(message["user_id"], str):
                raise UserIDException("User ID is not a string")
            elif message["user_id"] == "":
                raise UserIDException("User ID is empty")
            user_id = message["user_id"]
            logger.info(f"user_id={user_id}")

            with app.app_context():
                if accepting_user_requests.get(user_id) is False:
                    raise NotAcceptingRequestException("Not accepting requests from user")
                accepting_user_requests.set(user_id, False)
            user_data = None

            if "blocks" not in message:
                raise BadRequestException("Blocks is missing")
            elif not isinstance(message["blocks"], list):
                raise BadRequestException("Blocks is not a list")
            blocks = message["blocks"]
            dataset = None
            executed_blocks = []

            # TODO consider dispatching this whole pipeline to another thread, which seems to be best way to approach
            # TODO cpu intense tasks with gevents. Probably best to do this after we split the code out better.
            # https://github.com/miguelgrinberg/Flask-SocketIO/issues/1473
            for current_block_info in blocks:
                logger.info(f"Executing block={current_block_info} user={user_id}")

                # Allow other threads to execute
                # https://stackoverflow.com/questions/30901998/threading-true-with-flask-socketio
                gevent.sleep(0.1)

                if "block_id" not in current_block_info:
                    raise BadRequestException("Block ID is missing")
                elif current_block_info["block_id"] == "loaddata" and not executed_blocks:
                    user_data, dataset, output_message = load_data(user_data, current_block_info)
                elif current_block_info["block_id"] == "basicfiltering" and "loaddata" in executed_blocks:
                    block = BasicFiltering()
                    user_data, output_message = block.run(user_data, current_block_info)
                elif current_block_info["block_id"] == "qcplots" and "loaddata" in executed_blocks:
                    block = QCPlots()
                    user_data, output_message = block.run(user_data, dataset, current_block_info)
                elif current_block_info["block_id"] == "qcfiltering" and "loaddata" in executed_blocks:
                    block = QCFiltering()
                    user_data, output_message = block.run(user_data, dataset, current_block_info)
                elif current_block_info["block_id"] == "variablegenes" and "loaddata" in executed_blocks:
                    block = VariableGenes()
                    user_data, output_message = block.run(user_data, current_block_info)
                elif current_block_info["block_id"] == "pca" and "variablegenes" in executed_blocks:
                    block = PCA()
                    user_data, output_message = block.run(user_data, current_block_info)
                elif current_block_info["block_id"] == "integration" and "pca" in executed_blocks:
                    block = Integration()
                    user_data, output_message = block.run(user_data, dataset, current_block_info)
                elif current_block_info["block_id"] == "runumap" and "pca" in executed_blocks:
                    block = RunUMAP()
                    user_data, output_message = block.run(user_data, current_block_info)
                elif current_block_info["block_id"] in ["loaddata", "basicfiltering", "qcplots", "qcfiltering", "variablegenes", "pca", "integration", "runumap"]:
                    raise BadRequestException("Blocks have an invalid order")
                else:
                    raise BadRequestException("Block ID does not match expected values")

                executed_blocks.append(current_block_info["block_id"])

                output_message["blockId"] = current_block_info["block_id"]
                socketio.emit("json", json.dumps(output_message), room=client)
                logger.debug("emitted:" + json.dumps(output_message))

                # Allow other threads to execute
                # https://stackoverflow.com/questions/30901998/threading-true-with-flask-socketio
                gevent.sleep(0.1)

            logger.debug(f"Finished processing blocks for user={user_id}")
            end_connection = {"end_connection": "end_connection"}
            with app.app_context():
                accepting_user_requests.set(user_id, True)

        except UserIDException as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": "Your UserID is invalid: " + e.message + ". Please refresh the page and try again."}
        except NotAcceptingRequestException as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": "You have another request in progress. Please wait and try again."}
        except BadRequestException as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": "Received a bad request: " + e.message + ". Please refresh the page and try again."}
            with app.app_context():
                accepting_user_requests.set(user_id, True)
        except MissingParametersException as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": e.message + ". Please refresh the page and try again."}
            with app.app_context():
                accepting_user_requests.set(user_id, True)
        except Exception as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": "Unknown error. Please refresh the page and try again."}
            with app.app_context():
                accepting_user_requests.set(user_id, True)
        finally:
            socketio.emit("json", json.dumps(end_connection), room=client)

    def load_data(user_data, block):
        missing_parameters = get_missing_parameters(["dataset"], block)
        if missing_parameters:
            logger.error("Parameters missing from Block: " + json.dumps(missing_parameters))
            raise MissingParametersException("Missing parameters: " + json.dumps(missing_parameters))

        dataset = block["dataset"]
        if dataset in raw_data_cache:
            user_data = raw_data_cache.get(dataset).copy()
        else:
            raise Exception("Selected dataset does not exist.")

        message = {
            "text": adata_text(user_data)
        }
        return user_data, dataset, message

    def get_missing_parameters(params, block):
        missing_params = []
        for param in params:
            if param not in block:
                missing_params.append(param)
        return missing_params

    return socketio, app


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--production-mode", default=False, action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    production_mode = args.production_mode

    if production_mode:
        # Production configuration
        socketio, app = create_app(False)
        socketio.run(app, port=8080)
    else:
        # Testing configuration
        socketio, app = create_app(True)
        socketio.run(app)

    # By default, threading is handled by gevent.spawn:
    # https://github.com/miguelgrinberg/Flask-SocketIO/blob/40007fded6228013ce7e408ea1d9628da8b125fa/src/flask_socketio/__init__.py#L700C36-L700C42
    # https://www.gevent.org/api/gevent.baseserver.html#gevent.baseserver.BaseServer
