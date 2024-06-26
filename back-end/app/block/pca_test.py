import pytest
from block.pca import PCA
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


def test_run_callsScanpyFunctions(example_adata):
    block = PCA()
    input = {}

    with patch("scanpy.pp.calculate_qc_metrics", wraps=scanpy.pp.calculate_qc_metrics) as mock1, patch("scanpy.pp.scale", wraps=scanpy.pp.scale) as mock2, patch("scanpy.tl.pca", wraps=scanpy.tl.pca) as mock3, patch("scanpy.pl.pca_variance_ratio", wraps=scanpy.pl.pca_variance_ratio) as mock4:
        block.run(example_adata, input)
    mock1.assert_called_once()
    mock2.assert_called_once()
    mock3.assert_called_once()
    mock4.assert_called_once()
