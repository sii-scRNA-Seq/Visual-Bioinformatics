from unittest.mock import patch

import pytest
import scanpy

from block.exception.invalid_param_exception import InvalidParametersException
from block.exception.missing_param_exception import MissingParametersException
from block.plot_reddim import PlotDimensionReduction


@pytest.fixture()
def example_adata():
    example_adata = scanpy.datasets.pbmc3k()
    example_adata = scanpy.pp.subsample(example_adata, n_obs=500, random_state=0, copy=True)  # Subset 500 cells
    example_adata = example_adata[:, 0:200]  # Subset 200 genes
    print(example_adata)
    yield example_adata


def test_param_validation_requires_reduction():
    block = PlotDimensionReduction()
    with pytest.raises(MissingParametersException) as e_info:
        block.validate_parameters({
        })
    assert e_info.value.message == 'Missing parameters: ["reduction"]'


def test_run_throwsErrorIfReductionNotSupported(example_adata):
    block = PlotDimensionReduction()

    input = {
        "reduction": "AMAZING_MAP"
    }

    with pytest.raises(InvalidParametersException) as e_info:
        block.run(example_adata, "pbmc3k", input)

    assert e_info.value.message == 'AMAZING_MAP not supported'


def test_run_TSNE_callsScanpyFunctions(example_adata):
    block = PlotDimensionReduction()
    input = {
        "reduction": 'TSNE',
    }

    example_adata = example_adata.copy()
    example_adata = scanpy.tl.tsne(example_adata, copy=True)

    with patch("scanpy.tl.tsne", wraps=scanpy.tl.tsne) as mock1, \
            patch("scanpy.pl.tsne", wraps=scanpy.pl.tsne) as mock2:
        block.run(example_adata, "pbmc3k", input)

    mock1.assert_not_called()
    mock2.assert_called_once()
    mock2.assert_called_with(example_adata, color=None, show=False)


def test_run_TSNE_generatesTSNEifNotAlreadyGenerated(example_adata):
    block = PlotDimensionReduction()
    input = {
        "reduction": 'TSNE',
    }

    example_adata = example_adata.copy()

    with patch("scanpy.tl.tsne", wraps=scanpy.tl.tsne) as mock1, \
            patch("scanpy.pl.tsne", wraps=scanpy.pl.tsne) as mock2:
        block.run(example_adata, "pbmc3k", input)

    mock1.assert_called_once()
    mock2.assert_called_once()
    mock2.assert_called_with(example_adata, color=None, show=False)


def test_run_PCA_callsScanpyFunctions(example_adata):
    block = PlotDimensionReduction()
    input = {
        "reduction": 'PCA',
    }

    example_adata = example_adata.copy()
    scanpy.tl.pca(example_adata)

    with patch("scanpy.pl.pca", wraps=scanpy.pl.pca) as mock1:
        block.run(example_adata, "pbmc3k", input)

    mock1.assert_called_once()
    mock1.assert_called_with(example_adata, color=None, show=False)


def test_run_PCA_throwsErrorIfPCANotRun(example_adata):
    block = PlotDimensionReduction()

    input = {
        "reduction": "PCA"
    }

    with pytest.raises(InvalidParametersException) as e_info:
        block.run(example_adata, "pbmc3k", input)

    assert e_info.value.message == 'PCA not run yet'


def test_run_UMAP_callsScanpyFunctions(example_adata):
    block = PlotDimensionReduction()
    input = {
        "reduction": 'UMAP',
    }

    example_adata = example_adata.copy()
    scanpy.tl.pca(example_adata)
    scanpy.pp.neighbors(example_adata)
    scanpy.tl.umap(example_adata)

    with patch("scanpy.pl.umap", wraps=scanpy.pl.umap) as mock1:
        block.run(example_adata, "pbmc3k", input)

    mock1.assert_called_once()
    mock1.assert_called_with(example_adata, color=None, show=False)


def test_run_UMAP_throwsErrorIfUMAPNotRun(example_adata):
    block = PlotDimensionReduction()

    input = {
        "reduction": "UMAP"
    }

    with pytest.raises(InvalidParametersException) as e_info:
        block.run(example_adata, "pbmc3k", input)

    assert e_info.value.message == 'UMAP not run yet'


def test_run_usesLeidenColorsIfPresent(example_adata):
    block = PlotDimensionReduction()
    input = {
        "reduction": 'PCA',
    }

    example_adata = example_adata.copy()
    scanpy.tl.pca(example_adata, copy=False)
    scanpy.pp.neighbors(example_adata, copy=False)
    scanpy.tl.leiden(example_adata, copy=False)

    with patch("scanpy.pl.pca", wraps=scanpy.pl.pca) as mock1:
        block.run(example_adata, "pbmc3k", input)

    mock1.assert_called_once()
    mock1.assert_called_with(example_adata, color="leiden", show=False)
