import pytest
from block.integration import Integration
from block.exception.missing_param_exception import MissingParametersException
from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
import scanpy as sc


@pytest.fixture()
def example_adata():
    # Must be floats for PCA to work.
    # Harmony requires 30 cells at least
    counts = np.array([[5.0, 1.0, 2.0, 5.0, 3.0],
                       [0.0, 3.0, 4.0, 4.0, 3.0],
                       [2.0, 5.0, 0.0, 4.0, 3.0],
                       [0.0, 7.0, 8.0, 4.0, 3.0],
                       [0.0, 2.0, 0.0, 4.0, 3.0],
                       [1.0, 2.0, 3.0, 4.0, 3.0],
                       [5.0, 1.0, 2.0, 7.0, 1.0]])

    # Make this have the same cells * 6
    counts = np.vstack((counts, counts, counts, counts, counts, counts))
    counts = csr_matrix(counts)

    example_adata = AnnData(counts)
    example_adata.obs_names = [f"Cell_{i:d}" for i in range(example_adata.n_obs)]
    example_adata.var_names = [f"Gene_{i:d}" for i in range(example_adata.n_vars)]

    example_adata.obs["day"] = ["a", "b"] * 21
    sc.tl.pca(example_adata, n_comps=3, copy=False)

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


def test_run_hasCorrectReturnValues(example_adata):
    block = Integration()
    input = {
        "observation": "day"
    }

    result_adata, result_message = block.run(example_adata.copy(), "pf_dogga", input)

    assert result_adata.n_obs == 42
    assert result_adata.n_vars == 5
    assert "X_pca_harmony" in result_adata.obsm
    assert np.array_equal(result_adata.obsm["X_pca"], result_adata.obsm["X_pca_harmony"])
    assert result_message["text"] == "Object with: 42 cells and 5 genes"
