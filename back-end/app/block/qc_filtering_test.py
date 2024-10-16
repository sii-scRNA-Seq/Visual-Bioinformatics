import pytest
from block.qc_filtering import QCFiltering
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
    example_adata.obs["sample"] = ["1"] * example_adata.n_obs

    yield example_adata


def test_param_validation_requires_sample():
    block = QCFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "min_n_genes_by_counts": 42,
            "max_n_genes_by_counts": 42,
            "pct_counts_mt": 42
        })
    assert e_info.value.message == 'Missing parameters: ["sample"]'


def test_param_validation_requires_min_n_genes_by_counts():
    block = QCFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "sample": "sample",
            "max_n_genes_by_counts": 42,
            "pct_counts_mt": 42
        })
    assert e_info.value.message == 'Missing parameters: ["min_n_genes_by_counts"]'


def test_param_validation_requires_max_n_genes_by_counts():
    block = QCFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "sample": "sample",
            "min_n_genes_by_counts": 42,
            "pct_counts_mt": 42
        })
    assert e_info.value.message == 'Missing parameters: ["max_n_genes_by_counts"]'


def test_param_validation_requires_pct_counts_mt():
    block = QCFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
            "sample": "sample",
            "min_n_genes_by_counts": 42,
            "max_n_genes_by_counts": 42
        })
    assert e_info.value.message == 'Missing parameters: ["pct_counts_mt"]'


def test_param_validation_requires_all_parameters():
    block = QCFiltering()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({})
    assert e_info.value.message == 'Missing parameters: ["sample", "min_n_genes_by_counts", "max_n_genes_by_counts", "pct_counts_mt"]'


def test_run_correctlyIdentifiesMTGenesForPbmc3k(example_adata):
    block = QCFiltering()
    input = {
        "sample": "sample",
        "min_n_genes_by_counts": 42,
        "max_n_genes_by_counts": 42,
        "pct_counts_mt": 42
    }

    with patch("pandas.Index.str.startswith") as mock, pytest.raises(Exception):
        block.run(example_adata, "pbmc3k", input)
    mock.assert_called_once_with("MT-")


def test_run_correctlyIdentifiesMTGenesForPf_dogga(example_adata):
    block = QCFiltering()
    input = {
        "sample": "sample",
        "min_n_genes_by_counts": 42,
        "max_n_genes_by_counts": 42,
        "pct_counts_mt": 42
    }

    with patch("pandas.Index.str.contains") as mock, pytest.raises(Exception):
        block.run(example_adata, "pf_dogga", input)
    mock.assert_called_once_with("MIT")


def test_run_raisesExceptionIfDatasetDoesNotExist(example_adata):
    block = QCFiltering()
    input = {
        "sample": "sample",
        "min_n_genes_by_counts": 42,
        "max_n_genes_by_counts": 42,
        "pct_counts_mt": 42
    }

    with pytest.raises(Exception) as e_info:
        block.run(example_adata, "", input)
    assert str(e_info.value) == "Selected dataset does not exist."


def test_run_raisesExceptionIfSampleDoesNotExist(example_adata):
    block = QCFiltering()
    input = {
        "sample": "ThisSampleDoesNotExist",
        "min_n_genes_by_counts": 42,
        "max_n_genes_by_counts": 42,
        "pct_counts_mt": 42
    }

    with pytest.raises(Exception) as e_info:
        block.run(example_adata, "pbmc3k", input)
    assert str(e_info.value) == "Selected sample does not exist."


def test_run_callsScanpyFunctions(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 42,
        "max_n_genes_by_counts": 42,
        "pct_counts_mt": 42
    }

    with patch("scanpy.pp.calculate_qc_metrics", wraps=scanpy.pp.calculate_qc_metrics) as mock:
        block.run(example_adata, "pbmc3k", input)
    mock.assert_called_once()


def test_run_filtersOnlyCellsFromGivenSample(example_adata):
    block = QCFiltering()
    adata = example_adata.copy()
    adata.obs["sample"] = (["0"] * (adata.n_obs - 1)) + ["1"]
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 1000,
        "max_n_genes_by_counts": 5,
        "pct_counts_mt": 100
    }

    result_adata, result_message = block.run(adata, "pbmc3k", input)
    assert result_adata.n_vars == adata.n_vars
    assert result_adata.n_obs == adata.n_obs - 1
    assert result_adata.obs["sample"].to_list().count("0") == result_adata.n_obs
    assert result_adata.obs["sample"].to_list().count("1") == 0


def test_run_filtersNothingWhenValuesSetToBoundaries(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 0,
        "max_n_genes_by_counts": 5,
        "pct_counts_mt": 100
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == example_adata.n_obs


def test_run_filters_min_n_genes_by_counts_WhenValuesSetToNormalValues(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 1,
        "max_n_genes_by_counts": 5,
        "pct_counts_mt": 100
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == 5


def test_run_filters_max_n_genes_by_counts_WhenValuesSetToNormalValues(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 0,
        "max_n_genes_by_counts": 4,
        "pct_counts_mt": 100
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == 5


def test_run_filters_pct_counts_mt_WhenValuesSetToNormalValues(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 0,
        "max_n_genes_by_counts": 5,
        "pct_counts_mt": 95
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == 5


def test_run_filterAllWhenValuesSetToNormalValues(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 1,
        "max_n_genes_by_counts": 4,
        "pct_counts_mt": 95
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_vars == example_adata.n_vars
    assert result_adata.n_obs == 3


def test_run_hasCorrectReturnValues(example_adata):
    block = QCFiltering()
    input = {
        "sample": "1",
        "min_n_genes_by_counts": 1,
        "max_n_genes_by_counts": 4,
        "pct_counts_mt": 95
    }

    result_adata, result_message = block.run(example_adata.copy(), "pbmc3k", input)
    assert result_adata.n_obs == 3
    assert result_adata.n_vars == 4
    assert result_message["text"] == "Object with: 3 cells and 4 genes"
