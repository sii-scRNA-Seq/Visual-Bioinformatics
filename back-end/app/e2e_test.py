import pytest
from scipy.sparse import csr_matrix
from unittest.mock import call, patch
import anndata
import json
import numpy as np

from back_end import create_app


@pytest.fixture()
def socketio():
    socketio, app = create_app(test_mode=True)
    app.config.update({"TESTING": True})
    yield socketio


@pytest.fixture()
def app():
    socketio, app = create_app(test_mode=True)
    app.config.update({"TESTING": True})
    yield app


@pytest.fixture()
def socketio_client(socketio, app):
    return socketio.test_client(app)


@pytest.fixture()
def app_client(app):
    return app.test_client()


@pytest.fixture()
def runner(app):
    return app.test_cli_runner()


def get_AnnData(qc_filtering=False):
    if qc_filtering is False:
        counts = csr_matrix(np.array([[0, 1, 2],
                                      [0, 3, 4],
                                      [0, 5, 6],
                                      [0, 7, 8],
                                      [0, 0, 0]]))
        adata = anndata.AnnData(counts)
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
        adata.obs["sample"] = ["1"] * adata.n_obs
    else:
        counts = csr_matrix(np.array([[0, 1, 2, 0],
                                      [0, 3, 4, 1],
                                      [0, 5, 6, 2],
                                      [1, 0, 0, 0],
                                      [1, 1, 1, 1],
                                      [0, 1, 1, 98]]))
        adata = anndata.AnnData(counts)
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars - 1)] + ["MT-Gene"]
        adata.obs["sample"] = ["1"] * adata.n_obs
    return adata


# -------------------- TESTS BEGIN HERE -------------------- #


def test_getuserid_WarnsWhenUserIDIsNone(app_client):
    response = app_client.get("/api/getuserid")
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": "Bad Request",
        "description": "Not a valid user_id",
    }
    assert json.loads(response.data) == message


def test_getuserid_WarnsWhenUserIDIsNotAString(app_client):
    response = app_client.get("/api/getuserid", query_string={
        "user_id": []
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": "Bad Request",
        "description": "Not a valid user_id",
    }
    assert json.loads(response.data) == message


@patch("uuid.uuid4")
def test_getuserid_CreatesUserIdWhenUserIDIsEmpty(mock, app_client):
    mock.return_value = "bob"
    response = app_client.get("/api/getuserid", query_string={
        "user_id": ""
    })
    assert response.status_code == 200
    message = {
        "user_id": "bob"
    }
    assert json.loads(response.data) == message


def test_getuserid_ReturnsUserIdWhenUserIDIsSupplied(app_client):
    response = app_client.get("/api/getuserid", query_string={
        "user_id": "bob"
    })
    assert response.status_code == 200
    message = {
        "user_id": "bob"
    }
    assert json.loads(response.data) == message


def test_getdatasetinfo_ReturnsTheCorrectDatasetInfo(app_client):
    response = app_client.get("/api/getdatasetinfo")
    assert response.status_code == 200
    message = {
        "datasets": [
            {
                "key": "pbmc3k",
                "title": "Peripheral Blood Mononuclear Cells",
                "samples": ["1"],
                "integration_obs": [],
                "grouping_obs": ["leiden", "sample"]
            },
            {
                "key": "pf_dogga",
                "title": "Malaria Cell Atlas P. falciparum",
                "samples": ["0", "3", "5", "10a", "10b"],
                "integration_obs": ["day"],
                "grouping_obs": ["leiden", "day"]
            },
            {
                "key": "covid",
                "title": "COVID19 PBMCs",
                "samples": ["C1"],
                "integration_obs": ["Donor"],
                "grouping_obs": ["leiden", "Donor", "Status"]
            }
        ]
    }
    assert json.loads(response.data) == message


def test_executeblocks_WarnsWhenUserIDIsMissing(socketio_client):
    socketio_client.get_received()
    message = {}
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is missing. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["name"] == "results"
    assert len(received[0]["args"]) == 1
    assert received[0]["args"][0] == expected


def test_executeblocks_WarnsWhenUserIDIsNone(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": None,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is not a string. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["name"] == "results"
    assert len(received[0]["args"]) == 1
    assert received[0]["args"][0] == expected


def test_executeblocks_WarnsWhenUserIDIsNotAString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": 42,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is not a string. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["name"] == "results"
    assert len(received[0]["args"]) == 1
    assert received[0]["args"][0] == expected


def test_executeblocks_WarnsWhenUserIDIsEmptyString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "",
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is empty. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["name"] == "results"
    assert len(received[0]["args"]) == 1
    assert received[0]["args"][0] == expected


def test_executeblocks_ChangesValueOfAcceptingUserRequests(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("flask_caching.Cache.set") as mock:
        with patch("flask_caching.Cache.get", lambda x, y: True):
            with patch("block_execution.get_socketio_client", lambda x: socketio_client):
                with patch("block_execution.disconnect_socketio_client", lambda x: True):
                    with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):

                        socketio_client.emit("json", message)
                        received = socketio_client.get_received()

                        assert mock.call_count == 3
                        assert mock.mock_calls == [call("bob", True), call("bob", False), call("bob", True)]

                        json_received = list(filter(lambda x: x["name"] == "json", received))
                        assert len(json_received) == 1
                        result = json.loads(json_received[0]["args"])
                        assert result["end_connection"] == "end_connection"


def test_executeblocks_WarnsWhenNotAcceptingRequestsFromUser(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("flask_caching.Cache.get", lambda x, y: False):
        with patch("block_execution.get_socketio_client", lambda x: socketio_client):
            with patch("block_execution.disconnect_socketio_client", lambda x: True):
                with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                    socketio_client.emit("json", message)
                    received = socketio_client.get_received()

                    json_received = list(filter(lambda x: x["name"] == "results", received))
                    assert len(json_received) == 1
                    result = json.loads(json_received[0]["args"][0])
                    assert result["error"] == "You have another request in progress. Please wait and try again."


def test_executeblocks_AcceptsRequestsAfterEndConnection(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["end_connection"] == "end_connection"

                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["end_connection"] == "end_connection"


def test_executeblocks_WarnsWhenBlocksIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks is missing. Please refresh the page and try again."


def test_executeblocks_WarnsWhenBlocksIsNone(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": None,
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks is not a list. Please refresh the page and try again."


def test_executeblocks_WarnsWhenBlocksIsNotAList(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": 42,
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks is not a list. Please refresh the page and try again."


def test_executeblocks_WarnsWhenBlockIDIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [{}],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Block ID is missing. Please refresh the page and try again."


def test_executeblocks_WarnsWhenLoadDataIsAfterTheFirstBlock(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "loaddata", "dataset": "pbmc3k"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 2
                result = json.loads(json_received[1]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenBasicFilteringIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "basicfiltering", "min_genes": 0, "min_cells": 0}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenQcPlotsIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "qcplots"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenQcFilteringIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 0, "pct_counts_mt": 0}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenVariableGenesIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenPcaIsBeforeVariableGenes(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "pca"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 2
                result = json.loads(json_received[1]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenIntegrationIsBeforePca(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "integration", "observation": "day"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                result = json.loads(json_received[2]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenRunUmapIsBeforePca(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "runumap", "n_neighbors": 0, "n_pcs": 0}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                result = json.loads(json_received[2]["args"])
                assert result["error"] == "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."


def test_executeblocks_WarnsWhenBlockIDDoesNotMatchExpectedValues(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "Expecttheunexpected"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Received a bad request: Block ID does not match expected values. Please refresh the page and try again."


def test_executeblocks_WarnsForOtherErrors(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering", "min_genes": "min_genes", "min_cells": "min_cells"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 2
                result = json.loads(json_received[1]["args"])
                assert result["error"] == "Unknown error. Please refresh the page and try again."


def test_executeblocks_OnlyOneClientReceivesResponse():
    socketio, app = create_app(test_mode=True)
    app.config.update({"TESTING": True})
    client1 = socketio.test_client(app)
    client2 = socketio.test_client(app)

    client1.get_received()
    client2.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("block_execution.get_socketio_client", lambda x: client1, lambda x: client2):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                client1.emit("json", message)

                client1_received = client1.get_received()
                client1_json_received = list(filter(lambda x: x["name"] == "json", client1_received))

                client2_received = client2.get_received()
                client2_json_received = list(filter(lambda x: x["name"] == "json", client2_received))

    assert len(client1_json_received) == 1
    assert len(client2_json_received) == 0
    assert json.loads(client1_json_received[0]["args"])["end_connection"] == "end_connection"


def test_loaddata_WarnsWhenDatasetIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata"}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Missing parameters: [\"dataset\"]. Please refresh the page and try again."


def test_loaddata_AnnDataIsLoadedCorrectlyForPbmc3k(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"}
        ],
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 2
                assert json.loads(json_received[0]["args"])["blockId"] == "loaddata"
                assert json.loads(json_received[0]["args"])["text"] == "Object with: 2,700 cells and 32,738 genes"
                assert json.loads(json_received[1]["args"])["end_connection"] == "end_connection"


def test_loaddata_AnnDataIsLoadedCorrectlyForPfdogga(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"}
        ],
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 2
                assert json.loads(json_received[0]["args"])["blockId"] == "loaddata"
                assert json.loads(json_received[0]["args"])["text"] == "Object with: 10,000 cells and 5,720 genes"
                assert json.loads(json_received[1]["args"])["end_connection"] == "end_connection"


def test_loaddata_WarnsWhenDatasetDoesNotExist(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "DOES NOT EXIST"},
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                print(received)
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 1
                result = json.loads(json_received[0]["args"])
                assert result["error"] == "Unknown error. Please refresh the page and try again."


def test_basicfiltering_EndToEnd(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering", "min_genes": 1, "min_cells": 1}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                assert json.loads(json_received[1]["args"])["blockId"] == "basicfiltering"
                assert json.loads(json_received[1]["args"])["text"] == "Object with: 2,700 cells and 16,634 genes"
                assert json.loads(json_received[2]["args"])["end_connection"] == "end_connection"


def test_qcplots_EndToEnd(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcplots"}
        ],
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                assert json.loads(json_received[1]["args"])["blockId"] == "qcplots"
                assert len(json.loads(json_received[1]["args"])["image_list"]) == 3
                assert json.loads(json_received[1]["args"])["image_list"][0]["image"][:20] == "b'iVBORw0KGgoAAAANSU"
                assert json.loads(json_received[1]["args"])["image_list"][0][
                    "alt_text"] == "A violin plot displaying the distribution of the n_genes_by_counts observation generated by a QC Plots block"
                assert json.loads(json_received[1]["args"])["image_list"][1]["image"][:20] == "b'iVBORw0KGgoAAAANSU"
                assert json.loads(json_received[1]["args"])["image_list"][1][
                    "alt_text"] == "A violin plot displaying the distribution of the total_UMIs observation generated by a QC Plots block"
                assert json.loads(json_received[1]["args"])["image_list"][2]["image"][:20] == "b'iVBORw0KGgoAAAANSU"
                assert json.loads(json_received[1]["args"])["image_list"][2][
                    "alt_text"] == "A violin plot displaying the distribution of the pct_counts_mt observation generated by a QC Plots block"
                assert json.loads(json_received[2]["args"])["end_connection"] == "end_connection"


def test_qcfiltering_EndToEnd():
    with patch("scanpy.read_h5ad", lambda _: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "sample": "1", "min_n_genes_by_counts": 1, "max_n_genes_by_counts": 4, "pct_counts_mt": 95}
        ],
    }
    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                assert json.loads(json_received[1]["args"])["blockId"] == "qcfiltering"
                assert json.loads(json_received[1]["args"])["text"] == "Object with: 0 cells and 32,738 genes"
                assert json.loads(json_received[2]["args"])["end_connection"] == "end_connection"


def test_variablegenes_EndToEnd(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0}
        ],
    }

    with patch("block_execution.get_socketio_client", lambda x: socketio_client):
        with patch("block_execution.disconnect_socketio_client", lambda x: True):
            with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                socketio_client.emit("json", message)
                received = socketio_client.get_received()
                json_received = list(filter(lambda x: x["name"] == "json", received))
                assert len(json_received) == 3
                assert json.loads(json_received[1]["args"])["blockId"] == "variablegenes"
                assert json.loads(json_received[1]["args"])["image"][
                    "alt_text"] == "A scatter plot displaying dispersions of genes generated by an Identify Highly Variable Genes block"
                assert json.loads(json_received[2]["args"])["end_connection"] == "end_connection"


def test_pca_EndToEnd():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.read_h5ad", lambda _: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "pca"}
        ],
    }
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        with patch("block_execution.get_socketio_client", lambda x: socketio_client):
            with patch("block_execution.disconnect_socketio_client", lambda x: True):
                with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                    socketio_client.emit("json", message)
                    received = socketio_client.get_received()
                    json_received = list(filter(lambda x: x["name"] == "json", received))
                    assert len(json_received) == 4
                    assert json.loads(json_received[2]["args"])["blockId"] == "pca"
                    assert json.loads(json_received[2]["args"])["image"][
                        "alt_text"] == "A scatter plot displaying the contribution of each PC to the total variance in the data, generated by a Principal Component Analysis block"
                    assert json.loads(json_received[3]["args"])["end_connection"] == "end_connection"


def test_integration_EndToEnd(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "pca"},
            {"block_id": "integration", "observation": "day"}
        ],
    }

    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"), patch("scanpy.external.pp.harmony_integrate"), patch("anndata.AnnData.obsm"):
        with patch("block_execution.get_socketio_client", lambda x: socketio_client):
            with patch("block_execution.disconnect_socketio_client", lambda x: True):
                with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                    socketio_client.emit("json", message)
                    received = socketio_client.get_received()
                    json_received = list(filter(lambda x: x["name"] == "json", received))
                    assert len(json_received) == 5
                    assert json.loads(json_received[3]["args"])["blockId"] == "integration"
                    assert json.loads(json_received[3]["args"])["text"] == "Object with: 10,000 cells and 5,720 genes"
                    assert json.loads(json_received[4]["args"])["end_connection"] == "end_connection"


def test_runumap_EndToEnd():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.read_h5ad", lambda _: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "pca"},
            {"block_id": "runumap", "n_neighbors": 10, "n_pcs": 40}
        ],
    }

    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        with patch("block_execution.get_socketio_client", lambda x: socketio_client):
            with patch("block_execution.disconnect_socketio_client", lambda x: True):
                with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                    socketio_client.emit("json", message)
                    received = socketio_client.get_received()
                    json_received = list(filter(lambda x: x["name"] == "json", received))
                    assert len(json_received) == 5
                    assert json.loads(json_received[3]["args"])["blockId"] == "runumap"
                    assert json.loads(json_received[3]["args"])["image"][
                        "alt_text"] == "A UMAP of Leiden clusters using the principal components generated by a Run UMAP block"
                    assert json.loads(json_received[4]["args"])["end_connection"] == "end_connection"


def test_plotReducedDim_EndToEnd():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.read_h5ad", lambda _: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "pca"},
            {"block_id": "runumap", "n_neighbors": 10, "n_pcs": 40},
            {"block_id": "plot_reddim", "reduction": "UMAP", "group_by": "leiden"}
        ],
    }

    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        with patch("block_execution.get_socketio_client", lambda x: socketio_client):
            with patch("block_execution.disconnect_socketio_client", lambda x: True):
                with patch("block_execution.emit_synchronous", lambda cl, ch, msg: cl.emit(ch, msg)):
                    socketio_client.emit("json", message)
                    received = socketio_client.get_received()
                    json_received = list(filter(lambda x: x["name"] == "json", received))
                    assert len(json_received) == 6
                    assert json.loads(json_received[4]["args"])["blockId"] == "plot_reddim"
                    assert json.loads(json_received[4]["args"])["image"][
                        "alt_text"] == "A scatter plot of the data, reduced by UMAP."
                    assert json.loads(json_received[5]["args"])["end_connection"] == "end_connection"
