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


#################### TESTS BEGIN HERE #################### (# noqa: E266)


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


def test_executeblocks_WarnsWhenLoadDataIsAfterTheFirstBlock(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'loaddata'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_executeblocks_WarnsWhenBasicFilteringIsBeforeLoadData(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'basicfiltering', 'min_genes': 0, 'min_cells': 0}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenQcPlotsIsBeforeLoadData(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'qcplots'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenQcFilteringIsBeforeLoadData(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'qcfiltering', 'min_n_genes_by_counts': 0, 'max_n_genes_by_counts': 0, 'pct_counts_mt': 0}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenVariableGenesIsBeforeLoadData(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'variablegenes', 'min_mean': 0, 'max_mean': 0, 'min_disp': 0}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_executeblocks_WarnsWhenPcaIsBeforeVariableGenes(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'pca'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_executeblocks_WarnsWhenRunUmapIsBeforePca(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 0, 'max_mean': 0, 'min_disp': 0},
            {'block_id': 'runumap', 'n_neighbors': 0, 'n_pcs': 0}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 3
    assert received[2]["args"] == expected


def test_executeblocks_WarnsWhenBlockIDDoesNotMatchExpectedValues(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'Expecttheunexpected'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Received a bad request, please refresh the page and try again'})
    assert len(received) == 1
    assert received[0]["args"] == expected


def test_basicfiltering_WarnsWhenMinGenesAndMinCellsAreMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'basicfiltering'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_genes\', \'min_cells\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_basicfiltering_WarnsWhenMinGenesIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'basicfiltering', 'min_cells': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_genes\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_basicfiltering_WarnsWhenMinCellsIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'basicfiltering', 'min_genes': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_cells\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenAllParametersAreMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'qcfiltering'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_n_genes_by_counts\', \'max_n_genes_by_counts\', \'pct_counts_mt\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenMinNGenesByCountsIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'qcfiltering', 'max_n_genes_by_counts': 42, 'pct_counts_mt': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_n_genes_by_counts\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenMaxNGenesByCountsIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'qcfiltering', 'min_n_genes_by_counts': 42, 'pct_counts_mt': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'max_n_genes_by_counts\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_qcfiltering_WarnsWhenPctCountsMtIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'qcfiltering', 'min_n_genes_by_counts': 42, 'max_n_genes_by_counts': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'pct_counts_mt\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenAllParametersAreMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes'}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_mean\', \'max_mean\', \'min_disp\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMinMeanIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'max_mean': 42, 'min_disp': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_mean\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMaxMeanIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 42, 'min_disp': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'max_mean\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_variablegenes_WarnsWhenMinDispIsMissing(socketio_client, app_client):
    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 42, 'max_mean': 42}
        ],
    }
    socketio_client.emit('json', message)
    received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'min_disp\']'})
    assert len(received) == 2
    assert received[1]["args"] == expected


def test_runumap_WarnsWhenNNeighborsAndNPcsAreMissing():
    adata = get_AnnData(qc_filtering=True)
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    with patch('scanpy.datasets.pbmc3k', lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 0.0125, 'max_mean': 3, 'min_disp': 0.5},
            {'block_id': 'pca'},
            {'block_id': 'runumap'}
        ],
    }
    with patch('scanpy.pp.highly_variable_genes'), patch('scanpy.pl.highly_variable_genes'):
        socketio_client.emit('json', message)
        received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'n_neighbors\', \'n_pcs\']'})
    assert len(received) == 4
    assert received[3]["args"] == expected


def test_runumap_WarnsWhenNNeighborsIsMissing(socketio_client, app_client):
    adata = get_AnnData(qc_filtering=True)
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    with patch('scanpy.datasets.pbmc3k', lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 42, 'max_mean': 42, 'min_disp': 42},
            {'block_id': 'pca'},
            {'block_id': 'runumap', 'n_pcs': 42}
        ],
    }
    with patch('scanpy.pp.highly_variable_genes'), patch('scanpy.pl.highly_variable_genes'):
        socketio_client.emit('json', message)
        received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'n_neighbors\']'})
    assert len(received) == 4
    assert received[3]["args"] == expected


def test_runumap_WarnsWhenNPcsIsMissing(socketio_client, app_client):
    adata = get_AnnData(qc_filtering=True)
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    with patch('scanpy.datasets.pbmc3k', lambda: adata):
        socketio, app = create_app(test_mode=True)
        app.config.update({"TESTING": True})
        socketio_client = socketio.test_client(app)

    socketio_client.get_received()
    message = {
        'user_id': 'bob',
        'blocks': [
            {'block_id': 'loaddata'},
            {'block_id': 'variablegenes', 'min_mean': 42, 'max_mean': 42, 'min_disp': 42},
            {'block_id': 'pca'},
            {'block_id': 'runumap', 'n_neighbors': 42}
        ],
    }
    with patch('scanpy.pp.highly_variable_genes'), patch('scanpy.pl.highly_variable_genes'):
        socketio_client.emit('json', message)
        received = socketio_client.get_received()
    expected = json.dumps({'error': 'Missing parameters: [\'n_pcs\']'})
    assert len(received) == 4
    assert received[3]["args"] == expected


# def test(socketio_client, app_client):
#     socketio_client.get_received()
#     message = {
#         'user_id': 'abc',
#         'blocks': []
#     }
#     socketio_client.emit('json', message)
#     received = socketio_client.get_received()[0]["args"]
#     expected = json.dumps({'end_connection': 'end_connection'})
#     assert received == expected
