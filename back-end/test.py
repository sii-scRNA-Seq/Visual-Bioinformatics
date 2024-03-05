import pytest
from scipy.sparse import csr_matrix
from unittest.mock import patch
import anndata
import json
import numpy as np
from back_end import create_app


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
        'text': "Object with: 5 cells and 3 genes"
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
    response = client.get('/basicfiltering', query_string={
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
    response = client.get('/basicfiltering', query_string={
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
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
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
        'text': "Object with: 5 cells and 2 genes"
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
        'text': "Object with: 4 cells and 3 genes"
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
        'text': "Object with: 4 cells and 2 genes"
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
        'user_id': '',
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
        'user_id': 'bob',
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
        'user_id': 'bob',
    })
    response = client.get('/qcplots', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_qcplots_CallsScanpyFunctions(mock_loaddata, client):
    adata = get_AnnData()
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    mock_loaddata.return_value = adata

    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    with patch('scanpy.pp.calculate_qc_metrics') as mock1, patch('scanpy.pl.violin') as mock2:
        client.get('/qcplots', query_string={
            'user_id': 'bob',
        })
        mock1.assert_called_once()
        mock2.assert_called_once()


@patch('scanpy.read_10x_mtx')
def test_qcplots_ReturnsCorrectString(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcplots', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 200
    assert json.loads(response.data)['img'][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(response.data)['alttext'] == 'A violin plot displaying quality control metrics generated by a QC Plots block'


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
        'user_id': '',
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
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenGenesByCountsIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'pct_counts_mt': 0
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": "Missing parameters: ['min_n_genes_by_counts', 'max_n_genes_by_counts']",
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenPctCountsMtIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'min_n_genes_by_counts': 0,
        'max_n_genes_by_counts': 0,
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": "Missing parameters: ['pct_counts_mt']",
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenNGenesByCountsAndPctCountsMtAreMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": "Missing parameters: ['min_n_genes_by_counts', 'max_n_genes_by_counts', 'pct_counts_mt']",
    }
    assert json.loads(response.data) == message


def test_qcfiltering_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'min_n_genes_by_counts': 1,
        'max_n_genes_by_counts': 10,
        'pct_counts_mt': 1,
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_qcfiltering_maxNGenesByCountsWorks(mock_loaddata, client):
    mock_loaddata.return_value = get_AnnData(qc_filtering=True)
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'min_n_genes_by_counts': 0,
        'max_n_genes_by_counts': 4,
        'pct_counts_mt': 100,
    })
    assert response.status_code == 200
    message = {
        'text': "Object with: 4 cells and 4 genes"
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_qcfiltering_PctCountsMtWorks(mock_loaddata, client):
    mock_loaddata.return_value = get_AnnData(qc_filtering=True)
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'min_n_genes_by_counts': 0,
        'max_n_genes_by_counts': 5,
        'pct_counts_mt': 95,
    })
    assert response.status_code == 200
    message = {
        'text': "Object with: 4 cells and 4 genes"
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_qcfiltering_maxNGenesByCountsAndPctCountsMtWork(mock_loaddata, client):
    mock_loaddata.return_value = get_AnnData(qc_filtering=True)
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/qcfiltering', query_string={
        'user_id': 'bob',
        'min_n_genes_by_counts': 1,
        'max_n_genes_by_counts': 4,
        'pct_counts_mt': 95,
    })
    assert response.status_code == 200
    message = {
        'text': "Object with: 3 cells and 4 genes"
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/variablegenes')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/variablegenes', query_string={
        'user_id': '',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/variablegenes', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenMinMeanAndMaxMeanAreMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/variablegenes', query_string={
        'user_id': 'bob',
        'min_disp': 0,
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_mean\', \'max_mean\']',
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenMinDispIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/variablegenes', query_string={
        'user_id': 'bob',
        'min_mean': 0,
        'max_mean': 0,
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_disp\']',
    }
    assert json.loads(response.data) == message


def test_variablegenes_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/variablegenes', query_string={
        'user_id': 'bob',
        'min_mean': 0,
        'max_mean': 0,
        'min_disp': 0,
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_variablegenes_CallsScanpyFunctions(mock_loaddata, client):
    mock_loaddata.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    with patch('scanpy.pp.normalize_total') as mock1, patch('scanpy.pp.log1p') as mock2, patch('scanpy.pp.highly_variable_genes') as mock3, patch('scanpy.pl.highly_variable_genes') as mock4:
        client.get('/variablegenes', query_string={
            'user_id': 'bob',
            'min_mean': 0,
            'max_mean': 0,
            'min_disp': 0,
        })
        mock1.assert_called_once()
        mock2.assert_called_once()
        mock3.assert_called_once()
        mock4.assert_called_once()


@patch('scanpy.read_10x_mtx')
def test_variablegenes_ReturnsCorrectString(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/variablegenes', query_string={
        'user_id': 'bob',
        'min_mean': 0,
        'max_mean': 0,
        'min_disp': 0,
    })
    assert response.status_code == 200
    assert json.loads(response.data)['img'][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(response.data)['alttext'] == 'A scatter plot displaying dispersions of genes generated by an Identify Highly Variable Genes block'


def test_pca_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/pca')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_pca_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/pca', query_string={
        'user_id': '',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_pca_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/pca', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_pca_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/pca', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_pca_CallsScanpyFunctions(mock_loaddata, client):
    adata = get_AnnData(qc_filtering=True)
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    mock_loaddata.return_value = adata

    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    with patch('scanpy.pp.calculate_qc_metrics') as mock1, patch('scanpy.pp.regress_out') as mock2, patch('scanpy.pp.scale') as mock3, patch('scanpy.tl.pca') as mock4, patch('scanpy.pl.pca_variance_ratio') as mock5:
        client.get('/pca', query_string={
            'user_id': 'bob',
        })
        mock1.assert_called_once()
        mock2.assert_called_once()
        mock3.assert_called_once()
        mock4.assert_called_once()
        mock5.assert_called_once()


@patch('scanpy.read_10x_mtx')
def test_pca_ReturnsCorrectString(mock, client):
    adata = get_AnnData(qc_filtering=True)
    adata.obs['total_counts'] = list(range(0, adata.n_obs))
    mock.return_value = adata

    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/pca', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 200
    assert json.loads(response.data)['img'][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(response.data)['alttext'] == 'A scatter plot displaying the contribution of each PC to the total variance in the data, generated by a Principle Component Analysis block'


def test_runumap_WarnsUserWhenUserIdIsNone(client):
    response = client.get('/runumap')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenUserIdIsEmpty(client):
    response = client.get('/runumap', query_string={
        'user_id': '',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenUserIdIsNotInCache(client):
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Not a valid user_id',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenNNeighborsIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
        'n_pcs': 0,
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'n_neighbors\']',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenNPcsIsMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
        'n_neighbors': 0,
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'n_pcs\']',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenNNeighborsAndNPcsAreMissing(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
    })
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'n_neighbors\', \'n_pcs\']',
    }
    assert json.loads(response.data) == message


def test_runumap_WarnsUserWhenNoDataIsInUserCache(client):
    client.get('/getuserid', query_string={
        'user_id': 'bob',
    })
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
        'n_neighbors': 0,
        'n_pcs': 0,
    })
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the blocks and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_runumap_CallsScanpyFunctions(mock_loaddata, client):
    mock_loaddata.return_value = get_AnnData()
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    with patch('scanpy.pp.neighbors') as mock1, patch('scanpy.tl.umap') as mock2, patch('scanpy.tl.leiden') as mock3, patch('scanpy.pl.umap') as mock4:
        client.get('/runumap', query_string={
            'user_id': 'bob',
            'n_neighbors': 0,
            'n_pcs': 0,
        })
        mock1.assert_called_once()
        mock2.assert_called_once()
        mock3.assert_called_once()
        mock4.assert_called_once()


@patch('scanpy.read_10x_mtx')
def test_runumap_ReturnsCorrectString(mock, client):
    mock.return_value = get_AnnData(qc_filtering=True)
    client.get('/loaddata', query_string={
        'user_id': 'bob',
    })
    response = client.get('/runumap', query_string={
        'user_id': 'bob',
        'n_neighbors': 10,
        'n_pcs': 40,
    })
    assert response.status_code == 200
    assert json.loads(response.data)['img'][:20] == "b'iVBORw0KGgoAAAANSU"
    assert json.loads(response.data)['alttext'] == 'A UMAP of Leiden clusters using the principle components generated by a Run UMAP block'
