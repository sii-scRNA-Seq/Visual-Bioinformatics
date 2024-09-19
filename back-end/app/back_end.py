from gevent import monkey

from exception.not_accepting_request_exception import NotAcceptingRequestException
from exception.user_id_exception import UserIDException

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

from dataset_info import dataset_info

from tasks import execute_blocks, execute_blocks_celery

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
    """
    Create the SocketIO and Flask apps.

    Parameters:

        - `test_mode`: Will be True if creating an app for testing, and False if creating an app for production.

    Return:

        - The SocketIO app.
        - The Flask app.
    """

    if test_mode:
        logger = logging.getLogger("scampi-test")
        user_cache_config = file_cache_config
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

    app = Flask(__name__, static_folder="../dist/scampi/browser", static_url_path="/dist/scampi/browser")

    accepting_user_requests = Cache(config=user_cache_config)
    accepting_user_requests.init_app(app)

    socketio = SocketIO(app, cors_allowed_origins="*")

    if test_mode:
        # Only required in local mode as in production we serve the frontend from flask.
        CORS(app)

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
        """
        A function called by a HTTP request to /api/getuserid.

        Receive a `user_id` argument. If it is not a string, raise an error; if it is empty, create a new user id; otherwise, use the existing user id.

        Return:

            - A JSON response object containing a user id.

        Exceptions: `werkzeug.exceptions.BadRequest`.
        """
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
        """
        A function called by a HTTP request to /api/getdatasetinfo.

        Send a list of available datasets, with information such as names and samples, to the frontend for configuration.

        Return:

            - A JSON response object containing the dataset info.
        """
        logger.info("Sending dataset info")
        message = {
            "datasets": dataset_info
        }
        return jsonify(message)

    @socketio.on("results")
    def return_results_to_correct_user(message):
        logger.info(f'Emitting user results {message}')
        room = message["session_id"]
        user_id = message["user_id"]
        socketio.emit("json", json.dumps(message), room=room)
        socketio.emit("handled_users", json.dumps({"user_id": user_id}))

    @socketio.on("json")
    def dispatch_execute_blocks(message):
        logger.info(f"Dispatching blocks, json={message}")
        sid = request.sid

        try:
            if "user_id" not in message:
                raise UserIDException("User ID is missing")
            elif not isinstance(message["user_id"], str):
                raise UserIDException("User ID is not a string")
            elif message["user_id"] == "":
                raise UserIDException("User ID is empty")
            user_id = message["user_id"]

            # with app.app_context():
            #     if accepting_user_requests.get(user_id) is False:
            #         raise NotAcceptingRequestException("Not accepting requests from user")
            #     accepting_user_requests.set(user_id, False)

            if test_mode:
                execute_blocks(message, sid, f"ws://127.0.0.1:5000")
            else:
                execute_blocks_celery.delay(message, sid, f"ws://127.0.0.1:8080")

        except UserIDException as e:
            logger.error(e, exc_info=True)
            end_connection = {
                "error": "Your UserID is invalid: " + e.message + ". Please refresh the page and try again."}
            socketio.emit("json", json.dumps(end_connection), room=sid)
        except NotAcceptingRequestException as e:
            logger.error(e, exc_info=True)
            end_connection = {"error": "You have another request in progress. Please wait and try again."}
            socketio.emit("json", json.dumps(end_connection), room=sid)

    return socketio, app


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--production-mode", default=False, action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    production_mode = args.production_mode
    # `production_mode` will be True if the --production_mode flag was used, or False otherwise

    if production_mode:
        # Start the server in production mode
        socketio, app = create_app(False)
        socketio.run(app, port=8080)
    else:
        # Start the server in test mode
        socketio, app = create_app(True)
        socketio.run(app)

    # By default, threading is handled by gevent.spawn:
    # https://github.com/miguelgrinberg/Flask-SocketIO/blob/40007fded6228013ce7e408ea1d9628da8b125fa/src/flask_socketio/__init__.py#L700C36-L700C42
    # https://www.gevent.org/api/gevent.baseserver.html#gevent.baseserver.BaseServer
