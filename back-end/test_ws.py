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


#################### TESTS BEGIN HERE ####################


def test_getuserid_WarnsWhenUserIDIsNone(app_client):
    response = app_client.get('/api/getuserid')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


@patch('uuid.uuid4')
def test_getuserid_CreatesUserIdWhenUserIDIsEmpty(mock, app_client):
    mock.return_value = 'bob'
    response = app_client.get('/api/getuserid', query_string={
        'user_id': ''
    })
    assert response.status_code == 200
    message = {
        'user_id': 'bob'
    }
    assert json.loads(response.data) == message


def test_getuserid_ReturnsUserIdWhenUserIDIsSupplied(app_client):
    response = app_client.get('/api/getuserid', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 200
    message = {
        'user_id': 'bob'
    }
    assert json.loads(response.data) == message


def test_executeblocks_WarnsWhenUserIDIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {}
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsNone(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': None,
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsNotAString(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 42,
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenUserIDIsEmptyString(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': '',
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Your UserID is invalid, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsNone(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': None,
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlocksIsNotAList(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': 42,
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlockIDIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [{}],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlockIDDoesNotMatchExpectedValues(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [{
            'block_id': 'Expecttheunexpected',
        }],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenBlockIsMissingParameters(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata',},
            {'block_id': 'basicfiltering'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Unknown error, please refresh the page and try again'})
    assert len(received) == 2
    assert received[1]["args"] == expected


### TEMPLATE ###
def test(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'abc',
        'blocks': []
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()[0]["args"]
    expected = json.dumps({'end_connection': 'end_connection'})
    assert received == expected
