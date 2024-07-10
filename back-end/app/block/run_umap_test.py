import pytest
from block.run_umap import RunUMAP
from block.exception.missing_param_exception import MissingParametersException
from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
from unittest.mock import patch

import scanpy


@pytest.fixture()
def example_adata():
    counts = csr_matrix(np.array([[0, 1, 2, 0],
                                  [0, 3, 4, 1],
                                  [0, 5, 6, 2],
                                  [1, 0, 0, 0],
                                  [1, 1, 1, 1],
                                  [0, 1, 1, 98]]))
    example_adata = AnnData(counts)
    example_adata.obs_names = [f"Cell_{i:d}" for i in range(example_adata.n_obs)]
    example_adata.var_names = [f"Gene_{i:d}" for i in range(example_adata.n_vars - 1)] + ["MT-Gene"]
    example_adata.obs["total_counts"] = list(range(0, example_adata.n_obs))

    yield example_adata


def test_param_validation_requires_n_neighbors():
    block = RunUMAP()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "n_pcs": 42
        })
    assert e_info.value.message == 'Missing parameters: ["n_neighbors"]'


def test_param_validation_requires_n_pcs():
    block = RunUMAP()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "n_neighbors": 42
        })
    assert e_info.value.message == 'Missing parameters: ["n_pcs"]'


def test_param_validation_requires_n_neighbors_and_n_pcs():
    block = RunUMAP()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({})
    assert e_info.value.message == 'Missing parameters: ["n_neighbors", "n_pcs"]'


def test_run_callsScanpyFunctions(example_adata):
    block = RunUMAP()
    input = {
        "n_neighbors": 42,
        "n_pcs": 42
    }

    with patch("scanpy.pp.neighbors", wraps=scanpy.pp.neighbors) as mock1, patch("scanpy.tl.umap", wraps=scanpy.tl.umap) as mock2, patch("scanpy.tl.leiden", wraps=scanpy.tl.leiden) as mock3, patch("scanpy.pl.umap", wraps=scanpy.pl.umap) as mock4:
        block.run(example_adata, input)
    mock1.assert_called_once()
    mock2.assert_called_once()
    mock3.assert_called_once()
    mock4.assert_called_once()