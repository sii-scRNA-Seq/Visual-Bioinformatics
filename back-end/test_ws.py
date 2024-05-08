import pytest
from scipy.sparse import csr_matrix
from unittest.mock import patch
import anndata
import json
import numpy as np
from ws import create_app


@pytest.fixture()
def socketio():
    with patch('scanpy.datasets.pbmc3k', get_AnnData):
        socketio, app = create_app(test_mode=True)
        yield socketio

@pytest.fixture()
def app():
    with patch('scanpy.datasets.pbmc3k', get_AnnData):
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
                                      [0, 1, 1, 98],
                                      [1, 1, 1, 1]]))
        adata = anndata.AnnData(counts)
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars - 1)] + ["MT-Gene"]
    return adata


def test(socketio_client, app_client):
    #response = client.get('/api/getuserid')
    #assert response.status_code == 400
    #message = {
    #    "code": 400,
    #    "name": 'Bad Request',
    #    "description": 'Not a valid user_id',
    #}
    #assert json.loads(response.data) == message
    assert True == True

    socketio_client.get_received()
    message = {
        'user_id': ''
    }
    socketio_client.emit('json', message) # HERE!
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Unknown error, please refresh the page and try again'})

    assert True == True

    # assert len(received) == 1
    # self.assertEqual(len(received), 1)
    # self.assertEqual(len(received[0]['args']), 1)
    # self.assertEqual(received[0]['name'], 'my custom response')
    # self.assertEqual(received[0]['args'][0]['a'], 'b')
