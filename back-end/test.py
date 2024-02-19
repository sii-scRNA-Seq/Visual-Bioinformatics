import pytest
from back_end import create_app
from unittest.mock import patch
import anndata
import numpy as np
import json
from scipy.sparse import csr_matrix


@pytest.fixture()
def app():
    app = create_app(test_mode=True)
    app.config.update({
        "TESTING": True,
    })
    yield app


@pytest.fixture()
def client(app):
    return app.test_client()


@pytest.fixture()
def runner(app):
    return app.test_cli_runner()


def get_AnnData():
    counts = csr_matrix(np.array([[0, 1, 2],
                                  [0, 3, 4],
                                  [0, 5, 6],
                                  [0, 7, 8],
                                  [0, 0, 0]]))
    adata = anndata.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    # adata.var_names[-1] = "MT-Gene", give one cell really high count in that column.
    return adata


def test_getuserid_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/getuserid')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


@patch('uuid.uuid4')
def test_getuserid_CreatesUserIdWhenUserIdIsEmpty(mock, client):
    mock.return_value = 'bob'
    response = client.get('/getuserid', query_string={
        'user_id': ''
    })
    assert response.status_code == 200
    message = {
        'user_id': 'bob'
    }
    assert json.loads(response.data) == message


def test_getuserid_ReturnsGivenUserId(client):
    response = client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 200
    message = {
        'user_id': 'bob'
    }
    assert json.loads(response.data) == message


def test_loaddata_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/loaddata')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_loaddata_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/loaddata', query_string={
        'user_id': ''
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_loaddata_AnnDataIsLoadedCorrectly(mock, client):
    mock.return_value = get_AnnData()
    response = client.get('/loaddata', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 200
    message = {
        'text': "AnnData object with n_obs × n_vars = 5 × 3"
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/basicfiltering')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/basicfiltering', query_string={
        'user_id': ''
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/basicfiltering', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenMinGenesIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('basicfiltering', query_string={
        'user_id': 'bob',
        'min_cells': 0
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_genes\']',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenMinCellsIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('basicfiltering', query_string={
        'user_id': 'bob',
        'min_genes': 0
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_cells\']',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenMinGenesAndMinCellsAreMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('basicfiltering', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_genes\', \'min_cells\']',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('/basicfiltering', query_string={
        'user_id': 'bob',
        'min_genes': 1,
        'min_cells': 1
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": ('The blocks you have executed are not a valid order. Please check the order and try again.'),
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterGenesWorks(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob'
    })
    response = client.get('/basicfiltering', query_string={
        'user_id': 'bob',
        'min_genes': 0,
        'min_cells': 1
    })
    assert response.status_code == 200
    message = {
        'text': ("AnnData object with n_obs × n_vars = 5 × 2\n    obs: 'n_genes'\n    var: 'n_cells'")
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterCellsWorks(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob'
    })
    response = client.get('/basicfiltering', query_string={
        'user_id': 'bob',
        'min_genes': 1,
        'min_cells': 0
    })
    assert response.status_code == 200
    message = {
        'text': ("AnnData object with n_obs × n_vars = 4 × 3\n    obs: 'n_genes'\n    var: 'n_cells'")
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterGenesAndCellsWork(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob'
    })
    response = client.get('/basicfiltering', query_string={
        'user_id': 'bob',
        'min_genes': 1,
        'min_cells': 1
    })
    assert response.status_code == 200
    message = {
        'text': ("AnnData object with n_obs × n_vars = 4 × 2\n    obs: 'n_genes'\n    var: 'n_cells'")
    }
    assert json.loads(response.data) == message


# @patch('scanpy.read_10x_mtx')
# def test_basicfiltering_UsesExistingDataForSameParameters(mock_loaddata, client):
#     mock_loaddata.return_value = get_AnnData()
#     client.get('/loaddata', query_string={
#         'user_id': 'bob'
#     })
#     with patch('scanpy.pp.filter_cells') as mock_fc, patch('scanpy.pp.filter_genes') as mock_fg:
#         client.get('/basicfiltering', query_string={
#             'user_id': 'bob',
#             'min_genes': 1,
#             'min_cells': 1
#         })
#         client.get('/basicfiltering', query_string={
#             'user_id': 'bob',
#             'min_genes': 1,
#             'min_cells': 1
#         })
#         mock_fc.assert_called_once()
#         mock_fg.assert_called_once()


def test_qcplots_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/qcplots')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcplots_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/qcplots', query_string={
        'user_id': ''
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcplots_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/qcplots', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcplots_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('/qcplots', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": ('The blocks you have executed are not a valid order. Please check the order and try again.'),
    }
    assert json.loads(response.data) == message


# @patch('scanpy.read_10x_mtx')
# def test_qcplots_ReturnsCorrectString(mock, client):
#     mock.return_value = get_AnnData()
#     client.get('/loaddata', query_string={
#         'user_id': 'bob'
#     })
#     response = client.get('/qcplots', query_string={
#         'user_id': 'bob',
#     })
#     assert response.status_code == 200
#     message = {
#         'text': ("AnnData object with n_obs × n_vars = 5 × 3\n    obs: 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'\n    var: 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'")
#     }
#     assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/qcfiltering')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/qcfiltering', query_string={
        'user_id': ''
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenNGenesByCountsIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('qcfiltering', query_string={
        'user_id': 'bob',
        'pct_counts_mt': 0
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'n_genes_by_counts\']',
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenPctCountsMtIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('qcfiltering', query_string={
        'user_id': 'bob',
        'n_genes_by_counts': 0
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'pct_counts_mt\']',
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenNGenesByCountsAndPctCountsMtAreMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('qcfiltering', query_string={
        'user_id': 'bob'
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'n_genes_by_counts\', \'pct_counts_mt\']',
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob'
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'n_genes_by_counts': 1,
        'pct_counts_mt': 1
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": ('The blocks you have executed are not a valid order. Please check the order and try again.'),
    }
    assert json.loads(response.data) == message


# @patch('scanpy.read_10x_mtx')
# @patch('scanpy.pp.calculate_qc_metrics')
# def test_qcfiltering_NGenesByCountsAndPctCountsMtWork(mocker, mock_loaddata, client):
#     mock_loaddata.return_value = get_AnnData()
#     client.get('/loaddata', query_string={
#         'user_id': 'bob'
#     })
#     #with patch('scanpy.pp.calculate_qc_metrics', wraps=sc.pp.calculate_qc_metrics(new_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)) as mock_qc:
#     spy = mocker.spy('scanpy.pp.calculate_qc_metrics')
#     response = client.get('/qcfiltering', query_string={
#         'user_id': 'bob',
#         'n_genes_by_counts': 1,
#         'pct_counts_mt': 1
#     })
#     assert response.status_code == 200
#     spy.assert_called_once()
