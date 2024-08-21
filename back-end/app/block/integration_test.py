import pytest
from block.integration import Integration
from block.exception.missing_param_exception import MissingParametersException
from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
from unittest.mock import patch


@pytest.fixture()
def example_adata():
    counts = csr_matrix(np.array([[0, 1, 2],
                                  [0, 3, 4],
                                  [0, 5, 6],
                                  [0, 7, 8],
                                  [0, 0, 0]]))
    example_adata = AnnData(counts)
    example_adata.obs_names = [f"Cell_{i:d}" for i in range(example_adata.n_obs)]
    example_adata.var_names = [f"Gene_{i:d}" for i in range(example_adata.n_vars)]

    yield example_adata


def test_param_validation_requires_observation():
    block = Integration()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({})
    assert e_info.value.message == 'Missing parameters: ["observation"]'


def test_run_raisesExceptionIfDatasetDoesNotExist(example_adata):
    block = Integration()
    input = {
        "observation": ""
    }

    with pytest.raises(Exception) as e_info:
        block.run(example_adata, "", input)
    assert str(e_info.value) == "Selected dataset does not exist."


def test_run_raisesExceptionIfObservationDoesNotExistForDataset(example_adata):
    block = Integration()
    input = {
        "observation": ""
    }

    with pytest.raises(Exception) as e_info:
        block.run(example_adata, "pbmc3k", input)
    assert str(e_info.value) == "Selected observation does not exist."


def test_run_callsScanpyFunctions(example_adata):
    block = Integration()
    input = {
        "observation": "day"
    }

    with patch("scanpy.external.pp.harmony_integrate") as mock, patch("anndata.AnnData.obsm"):
        block.run(example_adata, "pf_dogga", input)
    mock.assert_called_once()


def test_run_hasCorrectReturnValues(example_adata):
    block = Integration()
    input = {
        "observation": "day"
    }

    with patch("scanpy.external.pp.harmony_integrate"), patch("anndata.AnnData.obsm"):
        result_adata, result_message = block.run(example_adata.copy(), "pf_dogga", input)
    assert result_adata.n_obs == 5
    assert result_adata.n_vars == 3
    assert result_message["text"] == "Object with: 5 cells and 3 genes"
