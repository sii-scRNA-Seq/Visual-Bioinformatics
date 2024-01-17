import pytest
from back_end import create_app
from unittest.mock import patch
import anndata
import numpy as np
import json
from scipy.sparse import csr_matrix


@pytest.fixture()
def app():
    app = create_app()
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
    counts = csr_matrix(np.array([[0,1,2],[0,3,4],[0,5,6],[0,7,8],[0,0,0]]))
    adata = anndata.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    return adata


@patch('scanpy.read_10x_mtx')
def test_loaddata_AnnDataIsLoadedCorrectly(mock, client):
    mock.return_value = get_AnnData()
    response = client.get('/loaddata')
    assert response.status_code == 200
    message = {
            'text': "AnnData object with n_obs × n_vars = 5 × 3"
    }
    assert json.loads(response.data) == message


def test_basicfiltering_WarnsUserWhenLoadDataHasNotBeenCalled(client):
    response = client.get('/basicfiltering?min_genes=1&min_cells=1')
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the order and try again.',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterGenesWorks(mock, client):
    ## TODO: Refactor tests to use cache rather than calling loaddata
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('/basicfiltering?min_genes=0&min_cells=1')
    assert response.status_code == 200
    message = {
            'text': "AnnData object with n_obs × n_vars = 5 × 2\n    obs: 'n_genes'\n    var: 'n_cells'"
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterCellsWorks(mock, client):
    ## TODO: Refactor tests to use cache rather than calling loaddata
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('/basicfiltering?min_genes=1&min_cells=0')
    assert response.status_code == 200
    message = {
            'text': "AnnData object with n_obs × n_vars = 4 × 3\n    obs: 'n_genes'\n    var: 'n_cells'"
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_FilterGenesAndCellsWork(mock, client):
    ## TODO: Refactor tests to use cache rather than calling loaddata
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('/basicfiltering?min_genes=1&min_cells=1')
    assert response.status_code == 200
    message = {
            'text': "AnnData object with n_obs × n_vars = 4 × 2\n    obs: 'n_genes'\n    var: 'n_cells'"
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_WarnsUserWhenMinGenesIsMissing(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('basicfiltering?min_cells=0')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_genes\']',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_WarnsUserWhenMinCellsIsMissing(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('basicfiltering?min_genes=0')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_cells\']',
    }
    assert json.loads(response.data) == message


@patch('scanpy.read_10x_mtx')
def test_basicfiltering_WarnsUserWhenMinGenesAndMinCellsAreMissing(mock, client):
    mock.return_value = get_AnnData()
    client.get('/loaddata')
    response = client.get('basicfiltering')
    assert response.status_code == 400
    message = {
        "code": 400,
        "name": 'Bad Request',
        "description": 'Missing parameters: [\'min_genes\', \'min_cells\']',
    }
    assert json.loads(response.data) == message
