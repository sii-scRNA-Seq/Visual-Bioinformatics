import json
import logging
from exception.bad_request_exception import BadRequestException
from block.exception.missing_param_exception import MissingParametersException

from block.basic_filtering import BasicFiltering
from block.qc_plots import QCPlots
from block.qc_filtering import QCFiltering
from block.variable_genes import VariableGenes
from block.pca import PCA
from block.integration import Integration
from block.run_umap import RunUMAP

from block.block_interface import adata_text

import scanpy as sc
import socketio

logger = logging.getLogger()

logger.info("Loading raw data...")
raw_data_cache = {
    "pbmc3k": sc.read_h5ad("../datasets/pbmc3k/PBMC3K.h5ad"),
    "pf_dogga": sc.read_h5ad("../datasets/pf_dogga/MCA_PF_DOGGA.h5ad")
}
logger.info("Finished loading raw data")


def load_data(user_data, block):
    """
    Execute the code for a `LoadData` block.

    Extract the value for `dataset` from `block`, or raise an exception if it is missing.
    Create a copy of the AnnData for the dataset, or raise an exception if the dataset does not exist.

    Parameters:

        - `user_data`: The AnnData for which the code should be executed.
        - `block`: A dictionary mapping parameter names to their values, which should be used while executing the code.

    Return:

        - The resulting AnnData after performing the block's behaviour.
        - The dataset chosen by the user.
        - A dictionary containing the results that will be seen by the user.

    Exceptions: `MissingParametersException`, `Exception`.
    """
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
    """
        Find which parameters from `params` are not in `block`, to validate that the block contains the required parameters.

        Parameters:

            - `params`: A list of the required parameters.
            - `block`: A dictionary mapping block parameters to their values.

        Return:

            - A list of parameters which are in `params` but not in `block`
        """
    missing_params = []
    for param in params:
        if param not in block:
            missing_params.append(param)
    return missing_params


def get_socketio_client(socketio_url):
    """
    Constructs and connects a socketio client. We provide a function to mock for E2E tests.

    Parameters:

    - `socketio_url`: The URL of the websocket server.
    """
    client = socketio.SimpleClient()
    client.connect(socketio_url)
    return client


def disconnect_socketio_client(client):
    """
    Disconnects a socketio client. We provide a function to mock for E2E tests.

    Parameters:

    - `client`: A Websocket client
    """
    client.disconnect()


def emit_synchronous(client, channel, message):
    """
    Emits a message and blocks until it is sent. This forces the final result of execute blocks to be sent before worker
    closes. We provide a function to mock for E2E tests.

    Parameters:

    - `client`: A Websocket client
    - `channel`: The channel to send the message on
    - `message`: The message to send
    """
    client.call(channel, message)


def execute_blocks(message, session_id, socketio_url):
    """
    A function called by a WebSocket request.

    Extract `user_id` from `message` and raise an exception if it is missing or invalid.
    Raise an exception if not accepting requests, otherwise, stop accepting requests from this user.
    Extract `blocks` from `message` and raise an exception if it is missing or invalid.
    Process each individual block in `blocks`, emitting a JSON response to the client for each block.
    Begin accepting request from user again.
    Emit a final JSON to the client to end the connection.

    Parameters:

        - `message`: A dictionary containing the user's request.

    Exceptions: `UserIDException`, `NotAcceptingRequestException`, `BadRequestException`, `MissingParametersException`, `Exception`.
    """

    # standard Python
    sio = get_socketio_client(socketio_url)
    logger.info(f"Executing blocks, json={message}")

    try:
        user_id = message["user_id"]
        logger.info(f"user_id={user_id}")
        user_data = None

        if "blocks" not in message:
            raise BadRequestException("Blocks is missing")
        elif not isinstance(message["blocks"], list):
            raise BadRequestException("Blocks is not a list")
        blocks = message["blocks"]
        dataset = None
        executed_blocks = []

        for current_block_info in blocks:
            logger.info(f"Executing block={current_block_info} user={user_id}")

            if "block_id" not in current_block_info:
                raise BadRequestException("Block ID is missing")
            elif current_block_info["block_id"] == "loaddata" and not executed_blocks:
                user_data, dataset, output_message = load_data(user_data, current_block_info)
            elif current_block_info["block_id"] == "basicfiltering" and "loaddata" in executed_blocks:
                block = BasicFiltering()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "qcplots" and "loaddata" in executed_blocks:
                block = QCPlots()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "qcfiltering" and "loaddata" in executed_blocks:
                block = QCFiltering()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "variablegenes" and "loaddata" in executed_blocks:
                block = VariableGenes()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "pca" and "variablegenes" in executed_blocks:
                block = PCA()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "integration" and "pca" in executed_blocks:
                block = Integration()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] == "runumap" and "pca" in executed_blocks:
                block = RunUMAP()
                user_data, output_message = block.run(user_data, dataset, current_block_info)
            elif current_block_info["block_id"] in ["loaddata", "basicfiltering", "qcplots", "qcfiltering",
                                                    "variablegenes", "pca", "integration", "runumap"]:
                raise BadRequestException("Blocks have an invalid order")
            else:
                raise BadRequestException("Block ID does not match expected values")

            executed_blocks.append(current_block_info["block_id"])

            output_message["blockId"] = current_block_info["block_id"]
            output_message["session_id"] = session_id
            output_message["user_id"] = user_id
            sio.emit("results", output_message)
            logger.debug("emitted:" + json.dumps(output_message))

        logger.debug(f"Finished processing blocks for user={user_id}")
        end_connection = {"end_connection": "end_connection"}

    except BadRequestException as e:
        logger.error(e, exc_info=True)
        end_connection = {"error": "Received a bad request: " + e.message + ". Please refresh the page and try again."}
    except MissingParametersException as e:
        logger.error(e, exc_info=True)
        end_connection = {"error": e.message + ". Please refresh the page and try again."}
    except Exception as e:
        logger.error(e, exc_info=True)
        end_connection = {"error": "Unknown error. Please refresh the page and try again."}
    finally:
        end_connection["session_id"] = session_id
        end_connection["user_id"] = user_id
        emit_synchronous(sio, "results", end_connection)

    disconnect_socketio_client(sio)
