import pytest
from block.basic_filtering import BasicFiltering
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


def test_param_validation_requires_min_genes():
    block = BasicFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_cells": "25"
        })
    assert e_info.value.message == 'Missing parameters: ["min_genes"]'


def test_param_validation_requires_min_cells():
    block = BasicFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_genes": "25"
        })
    assert e_info.value.message == 'Missing parameters: ["min_cells"]'


def test_param_validation_requires_min_cells_and_min_genes():
    block = BasicFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
        })
    assert e_info.value.message == 'Missing parameters: ["min_genes", "min_cells"]'


def test_run_callsScanpyFunctions(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "25",
        "min_cells": "15"
    }

    with patch("scanpy.pp.filter_cells") as mock1, patch("scanpy.pp.filter_genes") as mock2:
        block.run(example_adata, input)
        mock1.assert_called_once()
        mock2.assert_called_once()


def test_run_filtersNothingWhenValuesSetToZero(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "0",
        "min_cells": "0"
    }

    result_adata, result_message = block.run(example_adata.copy(), input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == example_adata.n_obs


def test_run_filtersCellsWhenValuesSetToNormalValues(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "0",
        "min_cells": "1"
    }

    result_adata, result_message = block.run(example_adata.copy(), input)
    assert result_adata.n_vars == 2
    assert result_adata.n_obs == example_adata.n_obs


def test_run_filtersGenesWhenValuesSetToNormalValues(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "2",
        "min_cells": "0"
    }

    result_adata, result_message = block.run(example_adata.copy(), input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == 4


def test_run_filterBothWhenValuesSetToNormalValues(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "2",
        "min_cells": "2"
    }

    result_adata, result_message = block.run(example_adata.copy(), input)
    assert result_adata.n_vars == 2
    assert result_adata.n_obs == 4


def test_run_hasCorrectReturnValues(example_adata):
    block = BasicFiltering()
    input = {
        "min_genes": "1",
        "min_cells": "1"
    }

    result_adata, result_message = block.run(example_adata.copy(), input)
    assert result_adata.n_obs == 4
    assert result_adata.n_vars == 2
    assert result_message["text"] == "Object with: 4 cells and 2 genes"
