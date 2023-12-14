import pytest
from back_end import create_app
from unittest.mock import patch
import anndata
import json

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


@patch('scanpy.read_10x_mtx')
def test_loaddata(mock, client):
    mock.return_value = anndata.AnnData()
    response = client.get('/loaddata/')
    assert response.status_code == 200
    message = {
            'text': "AnnData object with n_obs × n_vars = 0 × 0",
            'other': ''
    }
    assert json.loads(response.data) == message


def test_basicfiltering(client):
    response = client.get('/basicfiltering/1/1/')
    assert response.status_code == 406
    message = {
        "code": 406,
        "name": 'Not Acceptable',
        "description": 'The blocks you have executed are not a valid order. Please check the order and try again.',
    }
    assert json.loads(response.data) == message
