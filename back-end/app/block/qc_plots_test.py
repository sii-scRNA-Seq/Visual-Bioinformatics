import pytest
from block.qc_plots import QCPlots
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


def test_run_correctlyIdentifiesMTGenesForPbmc3k(example_adata):
    block = QCPlots()
    input = {}

    with patch("pandas.Index.str.startswith") as mock, pytest.raises(Exception):
        block.run(example_adata, "pbmc3k", input)
    mock.assert_called_once_with('MT-')


def test_run_correctlyIdentifiesMTGenesForPf_dogga(example_adata):
    block = QCPlots()
    input = {}

    with patch("pandas.Index.str.contains") as mock, pytest.raises(Exception):
        block.run(example_adata, "pf_dogga", input)
    mock.assert_called_once_with('MIT')


def test_run_raisesExceptionIfDatasetDoesNotExist(example_adata):
    block = QCPlots()
    input = {}

    with pytest.raises(Exception) as e_info:
        block.run(example_adata, "", input)
    assert str(e_info.value) == "Selected dataset does not exist."


def test_run_callsScanpyFunctions(example_adata):
    block = QCPlots()
    input = {}

    with patch("scanpy.pp.calculate_qc_metrics", wraps=scanpy.pp.calculate_qc_metrics) as mock1, patch("scanpy.pl.violin", wraps=scanpy.pl.violin) as mock2:
        block.run(example_adata, "pbmc3k", input)
    mock1.assert_called_once()
    mock2.assert_called_once()