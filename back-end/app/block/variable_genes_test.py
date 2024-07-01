import pytest
from block.variable_genes import VariableGenes
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


def test_param_validation_requires_min_mean():
    block = VariableGenes()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "max_mean": 42,
            "min_disp": 42
        })
    assert e_info.value.message == 'Missing parameters: ["min_mean"]'


def test_param_validation_requires_max_mean():
    block = VariableGenes()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_mean": 42,
            "min_disp": 42
        })
    assert e_info.value.message == 'Missing parameters: ["max_mean"]'


def test_param_validation_requires_min_disp():
    block = VariableGenes()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_mean": 42,
            "max_mean": 42
        })
    assert e_info.value.message == 'Missing parameters: ["min_disp"]'


def test_param_validation_requires_all_parameters():
    block = VariableGenes()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({})
    assert e_info.value.message == 'Missing parameters: ["min_mean", "max_mean", "min_disp"]'


def test_run_callsScanpyFunctions(example_adata):
    block = VariableGenes()
    input = {
        "min_mean": 42,
        "max_mean": 42,
        "min_disp": 42
    }

    with patch("scanpy.pp.normalize_total", wraps=scanpy.pp.normalize_total) as mock1, patch("scanpy.pp.log1p", wraps=scanpy.pp.log1p) as mock2, patch("scanpy.pp.highly_variable_genes", wraps=scanpy.pp.highly_variable_genes) as mock3, patch("scanpy.pl.highly_variable_genes", wraps=scanpy.pl.highly_variable_genes) as mock4:
        block.run(example_adata, input)
    mock1.assert_called_once()
    mock2.assert_called_once()
    mock3.assert_called_once()
    mock4.assert_called_once()
