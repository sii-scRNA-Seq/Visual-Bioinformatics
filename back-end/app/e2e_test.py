import pytest
from scipy.sparse import csr_matrix
from unittest.mock import call, patch
import anndata
import json
import numpy as np

from back_end import create_app


@pytest.fixture()
def socketio():
    with patch("scanpy.datasets.pbmc3k", get_AnnData), patch("scanpy.read_h5ad", lambda _: get_AnnData()):
        socketio, app = create_app(test_mode=True)
        yield socketio


@pytest.fixture()
def app():
    with patch("scanpy.datasets.pbmc3k", get_AnnData), patch("scanpy.read_h5ad", lambda _: get_AnnData()):
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
            {"key": "pbmc3k", "title": "Peripheral Blood Mononuclear Cells", "integration_obs": []},
            {"key": "pf_dogga", "title": "Malaria Cell Atlas P. falciparum", "integration_obs": ["day"]}
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
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsNone(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": None,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is not a string. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsNotAString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": 42,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is not a string. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsEmptyString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "",
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Your UserID is invalid: User ID is empty. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_ChangesValueOfAcceptingUserRequests(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("flask_caching.Cache.set") as mock:
        socketio_client.emit("json", message)
        assert mock.call_count == 2
        assert mock.mock_calls == [call("bob", False), call("bob", True)]
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenNotAcceptingRequestsFromUser(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("flask_caching.Cache.get", lambda x, y: False):
        socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "You have another request in progress. Please wait and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_AcceptsRequestsAfterEndConnection(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected

    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_AcceptsRequestsAfterNotAcceptingRequestException(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [],
    }
    with patch("flask_caching.Cache.get", lambda x, y: False):
        socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "You have another request in progress. Please wait and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected

    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_AcceptsRequestsAfterBadRequestException(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks is missing. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected

    message = {
        "user_id": "bob",
        "blocks": [],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_AcceptsRequestsAfterMissingParametersException(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_genes\", \"min_cells\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected

    message = {
        "user_id": "bob",
        "blocks": [],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_AcceptsRequestsAfterOtherErrors(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering", "min_genes": "min_genes", "min_cells": "min_cells"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Unknown error. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected

    message = {
        "user_id": "bob",
        "blocks": [],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"end_connection": "end_connection"})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks is missing. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsNone(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": None,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks is not a list. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsNotAList(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": 42,
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks is not a list. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlockIDIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [{}],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Block ID is missing. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenLoadDataIsAfterTheFirstBlock(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "loaddata", "dataset": "pbmc3k"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_executeblocks_WarnsWhenBasicFilteringIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "basicfiltering", "min_genes": 0, "min_cells": 0}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenQcPlotsIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "qcplots"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenQcFilteringIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 0, "pct_counts_mt": 0}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenVariableGenesIsBeforeLoadData(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenPcaIsBeforeVariableGenes(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "pca"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


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
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 3
    assert received[2]["args"] == expected


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
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Blocks have an invalid order. Please refresh the page and try again."})
    assert len(received) == 3
    assert received[2]["args"] == expected


def test_executeblocks_WarnsWhenBlockIDDoesNotMatchExpectedValues(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "Expecttheunexpected"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Received a bad request: Block ID does not match expected values. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsForOtherErrors(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering", "min_genes": "min_genes", "min_cells": "min_cells"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Unknown error. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_executeblocks_OnlyOneClientReceivesResponse():
    adata = get_AnnData()
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
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
    client1.emit("json", message)
    client1_received = client1.get_received()
    client2_received = client2.get_received()
    assert len(client1_received) == 1
    assert len(client2_received) == 0
    assert client1_received[0]["args"] == json.dumps({"end_connection": "end_connection"})


def test_loaddata_WarnsWhenDatasetIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"dataset\"]. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_loaddata_AnnDataIsLoadedCorrectlyForPbmc3k(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 5 cells and 3 genes"})
    assert len(received) == 2
    assert received[0]["args"] == expected
    assert received[1]["args"] == json.dumps({"end_connection": "end_connection"})


def test_loaddata_AnnDataIsLoadedCorrectlyForPfdogga(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 5 cells and 3 genes"})
    assert len(received) == 2
    assert received[0]["args"] == expected
    assert received[1]["args"] == json.dumps({"end_connection": "end_connection"})


def test_loaddata_WarnsWhenDatasetDoesNotExist(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "DOES NOT EXIST"},
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Unknown error. Please refresh the page and try again."})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_basicfiltering_FilterCellsAndFilterGenesWork(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "basicfiltering", "min_genes": 1, "min_cells": 1}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 4 cells and 2 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcplots_CorrectlyIdentifiesMTGenesForPbmc3k():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcplots"}
        ],
    }
    with patch("pandas.Index.str.startswith") as mock:
        socketio_client.emit("json", message)
        mock.assert_called_once_with('MT-')


def test_qcplots_CorrectlyIdentifiesMTGenesForPfdogga():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"},
            {"block_id": "qcplots"}
        ],
    }
    with patch("pandas.Index.str.contains") as mock:
        socketio_client.emit("json", message)
        mock.assert_called_once_with('MIT')


# TODO: Test case where qc_plots() is called with a non-existent dataset (should raise an Exception)


def test_qcplots_CallsScanpyFunctions():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcplots"}
        ],
    }
    with patch("scanpy.pp.calculate_qc_metrics") as mock1, patch("scanpy.pl.violin") as mock2:
        socketio_client.emit("json", message)
        mock1.assert_called_once()
        mock2.assert_called_once()


def test_qcplots_ReturnsCorrectString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcplots"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    assert len(received) == 3
    assert json.loads(received[1]["args"])["img"][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(received[1]["args"])["alttext"] == "A violin plot displaying quality control metrics generated by a QC Plots block"
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcfiltering_WarnsWhenAllParametersAreMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_n_genes_by_counts\", \"max_n_genes_by_counts\", \"pct_counts_mt\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenMinNGenesByCountsIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "max_n_genes_by_counts": 42, "pct_counts_mt": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_n_genes_by_counts\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenMaxNGenesByCountsIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 42, "pct_counts_mt": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"max_n_genes_by_counts\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenPctCountsMtIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 42, "max_n_genes_by_counts": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"pct_counts_mt\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_CorrectlyIdentifiesMTGenesForPbmc3k():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 0, "pct_counts_mt": 0}
        ],
    }
    with patch("pandas.Index.str.startswith") as mock:
        socketio_client.emit("json", message)
        mock.assert_called_once_with('MT-')


def test_qcfiltering_CorrectlyIdentifiesMTGenesForPfdogga():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pf_dogga"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 0, "pct_counts_mt": 0}
        ],
    }
    with patch("pandas.Index.str.contains") as mock:
        socketio_client.emit("json", message)
        mock.assert_called_once_with('MIT')


# TODO: Test case where qc_filtering() is called with a non-existent dataset (should raise an Exception)


def test_qcfiltering_CallsScanpyFunctions():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 0, "pct_counts_mt": 0}
        ],
    }
    with patch("scanpy.pp.calculate_qc_metrics") as mock:
        socketio_client.emit("json", message)
        mock.assert_called_once()


def test_qcfiltering_NoFilteringWorks():
    with patch("scanpy.datasets.pbmc3k", lambda: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 5, "pct_counts_mt": 100}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 6 cells and 4 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcfiltering_FilterByMinNGenesByCountsWorks():
    with patch("scanpy.datasets.pbmc3k", lambda: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 1, "max_n_genes_by_counts": 5, "pct_counts_mt": 100}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 5 cells and 4 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcfiltering_FilterByMaxNGenesByCountsWorks():
    with patch("scanpy.datasets.pbmc3k", lambda: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 4, "pct_counts_mt": 100}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 5 cells and 4 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcfiltering_FilterByPctCountsMtWorks():
    with patch("scanpy.datasets.pbmc3k", lambda: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 0, "max_n_genes_by_counts": 5, "pct_counts_mt": 95}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 5 cells and 4 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_qcfiltering_FilterByAllParametersWorks():
    with patch("scanpy.datasets.pbmc3k", lambda: get_AnnData(qc_filtering=True)):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "qcfiltering", "min_n_genes_by_counts": 1, "max_n_genes_by_counts": 4, "pct_counts_mt": 95}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"text": "Object with: 3 cells and 4 genes"})
    assert len(received) == 3
    assert received[1]["args"] == expected
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_variablegenes_WarnsWhenAllParametersAreMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes"}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_mean\", \"max_mean\", \"min_disp\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMinMeanIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "max_mean": 42, "min_disp": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_mean\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMaxMeanIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 42, "min_disp": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"max_mean\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMinDispIsMissing(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 42, "max_mean": 42}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"min_disp\"]. Please refresh the page and try again."})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_CallsScanpyFunctions(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0}
        ],
    }
    with patch("scanpy.pp.normalize_total") as mock1, patch("scanpy.pp.log1p") as mock2, patch("scanpy.pp.highly_variable_genes") as mock3, patch("scanpy.pl.highly_variable_genes") as mock4:
        socketio_client.emit("json", message)
        mock1.assert_called_once()
        mock2.assert_called_once()
        mock3.assert_called_once()
        mock4.assert_called_once()


def test_variablegenes_ReturnsCorrectString(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0}
        ],
    }
    socketio_client.emit("json", message)
    received = socketio_client.get_received()
    assert len(received) == 3
    assert json.loads(received[1]["args"])["img"][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(received[1]["args"])["alttext"] == "A scatter plot displaying dispersions of genes generated by an Identify Highly Variable Genes block"
    assert received[2]["args"] == json.dumps({"end_connection": "end_connection"})


def test_pca_CallsScanpyFunctions():
    adata = get_AnnData()
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
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
    with patch("scanpy.pp.calculate_qc_metrics") as mock1, patch("scanpy.pp.scale") as mock2, patch("scanpy.tl.pca") as mock3, patch("scanpy.pl.pca_variance_ratio") as mock4:
        socketio_client.emit("json", message)
        mock1.assert_called_once()
        mock2.assert_called_once()
        mock3.assert_called_once()
        mock4.assert_called_once()


def test_pca_ReturnsCorrectString():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
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
        socketio_client.emit("json", message)
        received = socketio_client.get_received()
    assert len(received) == 4
    assert json.loads(received[2]["args"])["img"][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(received[2]["args"])["alttext"] == "A scatter plot displaying the contribution of each PC to the total variance in the data, generated by a Principal Component Analysis block"
    assert received[3]["args"] == json.dumps({"end_connection": "end_connection"})


def test_integration_ReturnsCorrectString(socketio_client):
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
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        with patch("scanpy.external.pp.harmony_integrate"), patch("anndata.AnnData.obsm"):
            socketio_client.emit("json", message)
            received = socketio_client.get_received()
            expected = json.dumps({"text": "Object with: 5 cells and 3 genes"})
            assert len(received) == 5
            assert received[3]["args"] == expected
            assert received[4]["args"] == json.dumps({"end_connection": "end_connection"})


def test_runumap_WarnsWhenNNeighborsAndNPcsAreMissing():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
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
            {"block_id": "runumap"}
        ],
    }
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        socketio_client.emit("json", message)
        received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"n_neighbors\", \"n_pcs\"]. Please refresh the page and try again."})
    assert len(received) == 4
    assert received[3]["args"] == expected


def test_runumap_WarnsWhenNNeighborsIsMissing():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 42, "max_mean": 42, "min_disp": 42},
            {"block_id": "pca"},
            {"block_id": "runumap", "n_pcs": 42}
        ],
    }
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        socketio_client.emit("json", message)
        received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"n_neighbors\"]. Please refresh the page and try again."})
    assert len(received) == 4
    assert received[3]["args"] == expected


def test_runumap_WarnsWhenNPcsIsMissing():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 42, "max_mean": 42, "min_disp": 42},
            {"block_id": "pca"},
            {"block_id": "runumap", "n_neighbors": 42}
        ],
    }
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        socketio_client.emit("json", message)
        received = socketio_client.get_received()
    expected = json.dumps({"error": "Missing parameters: [\"n_pcs\"]. Please refresh the page and try again."})
    assert len(received) == 4
    assert received[3]["args"] == expected


def test_runumap_CallsScanpyFunctions(socketio_client):
    socketio_client.get_received()
    message = {
        "user_id": "bob",
        "blocks": [
            {"block_id": "loaddata", "dataset": "pbmc3k"},
            {"block_id": "variablegenes", "min_mean": 0, "max_mean": 0, "min_disp": 0},
            {"block_id": "pca"},
            {"block_id": "runumap", "n_neighbors": 0, "n_pcs": 0}
        ],
    }
    with patch("scanpy.pp.highly_variable_genes"), patch("scanpy.pl.highly_variable_genes"):
        with patch("scanpy.pp.neighbors") as mock1, patch("scanpy.tl.umap") as mock2, patch("scanpy.tl.leiden") as mock3, patch("scanpy.pl.umap") as mock4:
            socketio_client.emit("json", message)
            mock1.assert_called_once()
            mock2.assert_called_once()
            mock3.assert_called_once()
            mock4.assert_called_once()


def test_runumap_ReturnsCorrectString():
    adata = get_AnnData(qc_filtering=True)
    adata.obs["total_counts"] = list(range(0, adata.n_obs))
    with patch("scanpy.datasets.pbmc3k", lambda: adata):
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
        socketio_client.emit("json", message)
        received = socketio_client.get_received()
    assert len(received) == 5
    assert json.loads(received[3]["args"])["img"][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(received[3]["args"])["alttext"] == "A UMAP of Leiden clusters using the principal components generated by a Run UMAP block"
    assert received[4]["args"] == json.dumps({"end_connection": "end_connection"})
